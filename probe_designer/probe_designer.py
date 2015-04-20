#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function
from itertools import groupby
from collections import defaultdict
import os
import traceback

from future.builtins import str
from future.builtins import zip
from future.builtins import map
from past.utils import old_div
from docopt import docopt
import pandas as pd
import arrow
from pathlib import Path
import dataset

import get_seq
import biosearch_designer
import blaster


doc = """
Usage: probe_designer.py TARGETS [-o=OUTPUT] [-m=MIN_PROBES] [-d=DEBUG] [-t=TIMEOUT] [-i=INPUT] [-or=ORGANISM]

Arguments:
    TARGETS Comma-separated list of genes to design probes towards

Options:
    -o --output Whether to write each passed probe to a csv file.
    -m --min_probes Minimum # of Probes to accept for a gene
    -d --debug Whether to print debug output
    -t --timeout How long(seconds) to wait for a response from NCBI.  Default=120 seconds.
    -i --input File containing comma-seperated list of TARGETS.  Overwrites targets
    -or --organism Organism to Target(Default: mouse)
"""

class Probe(object):
    def __init__(self, name, organism):
        self.name = name
        self.organism = organism

def maximize_masking(probes, max_probes=24):
    passed = defaultdict(dict)
    for k, v in groupby(probes, key=lambda x: x['Name']):
        for k1, v1 in groupby(v, key=lambda x: x['Masking']):
            v1 = list(v1)
            n_probes = sum([len(probe_set['Probes']) for probe_set in v1])
            if n_probes >= max_probes:
                passed[k][k1] = v1
    p2 = {}
    for gene, m_dist in passed.items():
        highest_key = sorted(list(m_dist.keys()), reverse=True)[0]
        p2[gene] = m_dist[highest_key]
    return p2

def batch_blast_probes(cds_org, debug, max_probes, min_probes, organism, probes, timeout, local=False):
    g_set = defaultdict(dict)
    db = dataset.connect("sqlite:///db/probes.db")
    p_table = db[organism]
    failed_probes = []
    mass_df = []
    for p_set in probes:
        gene = p_set[0]['Name']
        masking = p_set[0]['Masking']
        if any(p_table.find(Name=gene, Masking=masking)):
            print("Getting probeset from DB for {} at masking {}".format(gene, masking))
            g_set[gene][masking] = pd.DataFrame(list(p_table.find(Name=gene, Masking=masking)))
        merged_pset = []
        for cds_region in p_set:
            reg = cds_region['Probes']
            reg['CDS Region #'] = cds_region['CDS Region #']
            reg['Name'] = gene
            reg['Masking'] = masking
            merged_pset.append(reg)
        probe_df = pd.concat(merged_pset,ignore_index=True)
        if len(probe_df) < min_probes:
            continue
        name_tuple = list(
            zip(probe_df['Name'], probe_df['CDS Region #'].map(str), probe_df['Probe Position*'].map(str)))
        probe_df['Probe Name'] = list(map('-'.join, name_tuple))
        mass_df.append(probe_df)
    if not mass_df:
        return g_set
    df = pd.concat(mass_df, ignore_index=True, join='inner')
    #Choose first probe of every repeat.  Because multiple maskings are allowed, same probe can be chosen many times
    pre = len(df)
    unique_index = [df.index[df['Probe Name']==name][0] for name in df['Probe Name'].unique()]
    from random import shuffle
    shuffle(unique_index)
    df2 = df.ix[unique_index] #Do this in a random order to reduce # of small blast runs
    print(pre, len(df))
    #Try splitting probelist into smaller chunks probe sets
    chunk_size = 2000
    b_res = {}
    print("Blasting: {} total probes".format(len(df)))
    for i in range(0, len(df), chunk_size):
        sub_df = df2.iloc[i:i+chunk_size]
        blaster.blast_probes('', sub_df, timeout=1800, debug=False, organism=cds_org, local=local)
    for gene, sub_df in df.groupby(['Name']):
        b_res = blaster.check_db(sub_df, organism=cds_org)
        if len(b_res) == 0:
            raise Exception("Error:  Nothing found in blast DB for query.  Something is wrong.")
        for masking, v in sub_df.groupby(['Masking']):
            t_res = {probe:b_res.get(probe, [""]) for probe in v['Probe Name']}
            try:
                passed = blaster.filter_probes_based_on_blast(gene, t_res, v, max_probes=max_probes,
                                                                  min_probes=min_probes, debug=debug,
                                                                  max_false_hits=old_div(min_probes, 2))
            except blaster.BlastError as e:
                print(traceback.format_exc())
                print("Blast Failed for {}".format(gene))
                failed_probes.append(gene)
                continue
            print(len(passed), gene, masking)
            g_set[gene][masking] = passed
            # for n, line in passed.T.to_dict().items():
            #     p_table.insert(line)
            p_table.insert_many(passed.T.to_dict().values())
            print("Probeset {}, {} added to table".format(gene, masking))
    return g_set

def filter_probes(gene, sub_df, b_res, max_probes, min_probes, debug, max_false_hits):
    t_res = {probe:b_res.get(probe, [""]) for probe in v['Probe Name']}
    try:
        passed = blaster.filter_probes_based_on_blast(gene, t_res, v, max_probes=max_probes,
                                                            min_probes=min_probes, debug=debug,
                                                            max_false_hits=old_div(min_probes, 2))
    except blaster.BlastError as e:
        print("Blast Failed for {}".format(gene))
    return (gene, passed)


def blast_probes(cds_org, debug, max_probes, min_probes, organism, probes, timeout):
    g_set = defaultdict(dict)
    db = dataset.connect("sqlite:///db/probes.db")
    p_table = db[organism]
    failed_probes = []
    for p_set in probes:
        gene = p_set[0]['Name']
        masking = p_set[0]['Masking']
        if any(p_table.find(Name=gene, Masking=masking)):
            print("Getting probeset from DB for {} at masking {}".format(gene, masking))
            g_set[gene][masking] = pd.DataFrame(list(p_table.find(Name=gene, Masking=masking)))
        merged_pset = []
        for cds_region in p_set:
            reg = cds_region['Probes']
            reg['CDS Region #'] = cds_region['CDS Region #']
            reg['Name'] = gene
            merged_pset.append(reg)
        probe_df = pd.concat(merged_pset,ignore_index=True)
        if len(probe_df) < min_probes:
            continue
        name_tuple = list(
            zip(probe_df['Name'], probe_df['CDS Region #'].map(str), probe_df['Probe Position*'].map(str)))
        probe_df['Probe Name'] = list(map('-'.join, name_tuple))
        try:
            #First Blast
            b_res = blaster.blast_probes(gene, probe_df, timeout=timeout, debug=debug, organism=cds_org)
            #Find probeset based on blast
            passed = blaster.filter_probes_based_on_blast(gene, b_res, probe_df, max_probes=max_probes,
                                                          min_probes=min_probes, debug=debug,
                                                          max_false_hits=old_div(min_probes, 2))
        except blaster.BlastError as e:
            print(traceback.format_exc())
            print("Blast Failed for {}".format(gene))
            failed_probes.append(gene)
            continue
        g_set[gene][p_set[0]['Masking']] = passed
        passed['Masking'] = p_set[0]['Masking']
        for n, line in passed.T.to_dict().items():
            p_table.insert(line)
        print("Probeset added to table")
    # print("Failed on: {}".format(", ".join(failed_probes)))
    return g_set

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join([complement[c] for c in seq])[::-1]

def main(target_genes, max_probes=24, min_probes=24, timeout=120, debug=False, parallel=True, organism="mouse", probe_design='biosearch'):
    b_org = organism.lower()
    if "mouse" == b_org:
        cds_org = '"Mus musculus"[porgn:__txid10090]'
    elif "human" == b_org:
        cds_org = '"Homo sapiens"[porgn:__txid9606]'
    elif "rabbit" == b_org:
        cds_org = '"Oryctolagus cuniculus"[porgn:__txid9986]'
    else:
        raise IOError("Error organism {} not recognized.  Please add to program or reformat".format(organism))
    assert (isinstance(target_genes, list))
    passed_genes = {}
    for gene in target_genes:
        try:
            cds = get_seq.CDS(gene, cds_org).run()
            if len(''.join(cds['CDS List'])) > (20 * min_probes):
                rc_cds = [reverse_complement(c)
                                   for c in cds['CDS List']]
                cds['CDS List'] = rc_cds
                passed_genes[gene] = cds
            else:
                print("CDS too short of %s at %int" % (gene, len(''.join(cds['CDS List']))))
        except get_seq.CDSError as e:
            print(traceback.format_exc())
            print(e)
    if debug:
        write_folder = Path("debug").joinpath("cds")
        try:
            write_folder.mkdir(parents=True)
        except:
            pass
        write_path = write_folder.joinpath("{}.fasta".format(arrow.now()))
        print("Writing to file")
        with open(str(write_path), 'w') as f:
            print("Writing to file2")
            for gene, cds in passed_genes.items():
                f.write(">{} from {}\n".format(gene, organism))
                f.write("{}\n".format(cds))
        print("Writing to file3")
    if probe_design == 'biosearch':
        b_designer = biosearch_designer.Biosearch(organism=b_org, variants=True)
        probes = b_designer.design(passed_genes, min_probes, )
        b_designer.close()
    elif probe_design == 'oligoarray':
        import oligoarray_designer as od
        oda = od.Oligoarray()
        pd = dict(passed_genes.items()[:15])
        probes = oda.design(passed_genes, ml=35, Ml=35, max_tm=103, x_hyb=68)
        probes = [[probe] for probe in probes]
        print("DONE OLIGOAARAY")
    else:
        raise Exception("{} designer not recognized".format(probe_design))
    if debug:
        pass
        # write_folder = Path("debug").joinpath("bsearch")
        # try:
        #     write_folder.mkdir(parents=True)
        # except
        #
        #     pass
        # write_path = write_folder.joinpath("{}.csv".format(arrow.now()))
        # pd.DataFrame(probes).to_csv(str(write_path))
    g_set = batch_blast_probes(cds_org, debug, max_probes, min_probes, organism, probes, timeout)
    return g_set


def get_probes(gene, organism, masking):
    db = dataset.connect("sqlite:///db/probes.db")
    p_table = db[organism]
    return list(p_table.find(Name=gene, Masking=masking))


def get_optimal_probes(gene, organism, min_probes=12, max_probes=24):
    stored_probes = None
    for masking in (5, 4, 3):
        probes = list(get_probes(gene, organism, masking))
        if len(probes) == max_probes:
            return probes
        if len(probes) < min_probes:
            continue
        if stored_pset is not None:
            if len(stored_probes) - 2 > len(probes):
                stored_probes = probes
        stored_probes = probes
    if stored_probes is not None:
        return stored_probes
    raise Exception("No good probeset for {}".format(gene))

def del_probes(target_genes, organism='mouse'):
    db = dataset.connect("sqlite:///db/probes.db")
    p_table = db[organism]
    for gene in target_genes:
        p_table.delete(Name=gene)


# genes = """nrn1,s100a10,rbpms,fbxo2,nefl,tppp3,eg435376,sncg,rbpms2,rlbp1l2,slc17a6,chrnb3,bc089491,cpne9,tubb3,opn4,
# nppb,adcyap1,1700011l03Rik,prph,ctxn3,cd24a,prkcq,cdkn1c,gm687,trhde,igf1,rbms3,gm527,ndp,cyp1b1,rlbp1,
# kcnj13,rdh5,rpe65,rgr,efemp1,bbs2,bbs4"""
# genes = """grin1, grin2a, atp2b2, crmp1, bbc3, adora1, agrn, apoe, c19orf20, calb1, clock,
# dlg4, kifc2, mapk11, mib2, OXR1, pick1, pctk3, pyroxd1, aifm3, bad, bcl2, bok, app,
#         caskin2, grin2c, homer1, homer3, kcnip3, nlgn2, nrgn, nrxn3, amh, c4orf48, rest,
#         met, pten, gad1, efemp1, txnip, hcrtr1, penk, pnoc, crh"""

genes = """
Tbr1
,Rasgrf2
,Pvrl3
,Cux2
,Rorb
,Plcxd2
,Thsd7a
,Kcnk2
,Cplx3
,Sulf2
,Foxp2
,Syt6
,Rprm
,Nr4a2
,Synpr
,Lhx5
,Reln
,Rgs8
,Cartpt
,Tle1
,Tle3
,Mdga1
,Kcnip2
,Rorb
,Cyp39a1
,Lhx2
,Unc5d
,Gpr6
,Mef2c
,Dtx4
,Cux1
,Cux2
,Bhlhe22
,Plxnd1
,Pou3f1
,Marcksl1
,Pou3f3
,Pou3f2
,Nefh
,Cntn6
,Foxo1
,Opn3
,Lix1
,Syt9
,S100a10
,Oma1
,Ldb2
,Crim1
,Pcp4
,Rac3
,Bcl11b
,Crym
,Otx1
,Fezf2
,Igfbp4
,Sox5
,Dkk3
,Tle4
,Sema3e
,Nr4a3
,Lxn
,Foxp2
,Ppp1r1b
,Id2
,Slitrk1
,Tbr1
,Lmo4
,Lmo3
,Ctgf
,Nr4a2
,Wfs1
,DCN
,Htr2c
,Grp
,Gpr101
,Col5a1
,Gpc3
,Prss12
,Ndst4
,Calb1
,Matn2
,Rph3a
,Loxl1
,Plagl1
,Coch
,Itga7
,Lyt
,Gprin3
,Adora2a
,Cd4
,Rgs9
,Drd1
,Serpina9
,Ngfr
,Ntrk1
,Rarb
,Nts
,Pcdh11x
,Slc6a3
,Magell2
,Iyd,
Lyd,
Lm1,
Tiam2,
Robo1,
Cadps2
"""

genes = """
Albhk5
,Ash1
,Axin2
,Bmi1
,bmp4
,bmpr1a
,Brachyury
,Cdh1
,Col5a2
,ctcf
,dazl
,Dnmt1
,dnmt3a
,Dnmt3b
,Dnmt3L
,Dppa2
,Dppa3
,Dppa4
,EED
,ehmt2
,Esrrb
,Ezh2
,Fbxo15
,Fgf4
,Fgfr2
,FoxO1
,Foxa2
,G9A
,Gsk3b
,hdac1
,hdac2
,Hes1
,Id1
,Id2
,Jarid2
,kdm3a
,kdm3b
,kdm4c
,kdm4d
,kdm6a
,klf2
,klf4
,klf5
,Lifr
,Lin28a
,Lin28b
,LSD1
,mael
,METTL14
,METTL3
,Myc
,Mycn
,Nanog
,Nr0b1
,Oct4
,pax6
,Pecam1
,Prdm14
,Rest
,Sall4
,Sdha
,setdb1
,Smad1
,Smad5
,Smad4
,Socs3
,Sox2
,Sp1
,Suz12
,sycp3
,Tbx3
,Tcl1
,Tet1
,tet2
,tet3
,Trim28
,Thy1
,Utf1
,Xist
,Zfp42
"""

debug = True
timeout = 60


if __name__ == '__main__':
    args = docopt(doc, )
    target_genes = [x.strip() for x in args['TARGETS'].strip(',').split(",")]
    if args['--min_probes'] is not None:
        min_probes = int(args['--min_probes'])
    else:
        min_probes = 24
    if args['--debug'] is not None:
        debug = True
    else:
        debug = False
    if args['--timeout'] is not None:
        timeout = int(args['--timeout'])
    else:
        timeout = 120
    if args['--input'] is not None:
        target_genes = [x.strip() for x in args['--input'].strip(',').split(",")]
    if args['--organism'] is not None:
        organism = args['--organism']
    else:
        organism = 'mouse'
    probes = main(target_genes, min_probes, timeout, debug, organism)
    if args['--output'] is not False:
        for gene, probes in probes['Passed'].items():
            try:
                os.mkdir('passed_probes')
            except:
                pass
            out_path = os.path.join('passed_probes', gene + '.csv')
            probes.to_csv(out_path)
    else:
        print(probes['Passed'])
