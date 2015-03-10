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


def flatten_probes(p2):
    # Flatten out probe lists
    flat_probes = {}
    for gene, cds_psets in p2.items():
        for probeset in cds_psets:
            for k, v in probeset.items():
                if any([k == 'Probes', k == 'Target Seq']):
                    continue
                probeset['Probes'][k] = v
        p_list = [x['Probes'] for x in cds_psets]
        probes = pd.concat(p_list, ignore_index=True)
        flat_probes[gene] = probes

    # Give each probe a name
    for gene, probes in flat_probes.items():
        name_tuple = list(zip(probes['Name'], probes['CDS Region #'].map(str), probes['Probe Position*'].map(str)))
        probes['Probe Name'] = list(map('-'.join, name_tuple))
    return flat_probes


def batch_blast_probes(cds_org, debug, max_probes, min_probes, organism, probes, timeout):
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
            continue
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
    df = pd.concat(mass_df, ignore_index=True)
    print("Blasting:", len(df))
    b_res = blaster.blast_probes('', df, timeout=1800, debug=False, organism=cds_org)
    for k,v in df.groupby(['Name', 'Masking']):
        gene, masking = k
        t_res = {probe:b_res[probe] for probe in v['Probe Name']}
        try:
            passed = blaster.filter_probes_based_on_blast(gene, t_res, v, max_probes=max_probes,
                                                              min_probes=min_probes, debug=debug,
                                                              max_false_hits=old_div(min_probes, 2))
        except blaster.BlastError as e:
            print(traceback.format_exc())
            print("Blast Failed for {}".format(gene))
            failed_probes.append(gene)
            continue
        g_set[gene][k] = passed
        for n, line in passed.T.to_dict().items():
            p_table.insert(line)
        print("Probeset {}, {} added to table".format(*k))
    return g_set

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
            continue
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


def main(target_genes, max_probes=24, min_probes=24, timeout=120, debug=False, parallel=True, organism="mouse"):
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
    b_designer = biosearch_designer.Biosearch(organism=b_org, variants=True)
    probes = b_designer.design(passed_genes, min_probes, )
    b_designer.close()
    if debug:
        write_folder = Path("debug").joinpath("bsearch")
        try:
            write_folder.mkdir(parents=True)
        except:
            pass
        write_path = write_folder.joinpath("{}.csv".format(arrow.now()))
        pd.DataFrame(probes).to_csv(str(write_path))
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


# genes = """nrn1,s100a10,rbpms,fbxo2,nefl,tppp3,eg435376,sncg,rbpms2,rlbp1l2,slc17a6,chrnb3,bc089491,cpne9,tubb3,opn4,
# nppb,adcyap1,1700011l03Rik,prph,ctxn3,cd24a,prkcq,cdkn1c,gm687,trhde,igf1,rbms3,gm527,ndp,cyp1b1,rlbp1,
# kcnj13,rdh5,rpe65,rgr,efemp1,bbs2,bbs4"""
# genes = """grin1, grin2a, atp2b2, crmp1, bbc3, adora1, agrn, apoe, c19orf20, calb1, clock,
# dlg4, kifc2, mapk11, mib2, OXR1, pick1, pctk3, pyroxd1, aifm3, bad, bcl2, bok, app,
#         caskin2, grin2c, homer1, homer3, kcnip3, nlgn2, nrgn, nrxn3, amh, c4orf48, rest,
#         met, pten, gad1, efemp1, txnip, hcrtr1, penk, pnoc, crh"""

# genes = """nkx6-1
#         ,pdx1
#         ,mafb
#         ,pax6
#         ,ptf1a
#         ,sox9
#         ,ins
#         ,gcg
#         """
genes = """ins, pdx1, mafa, mafb, nkx6.1, nkx2.2, isl1, pax6, pax4, nd1,
pcsk1, sur1, ldha, mvk, hk1, hk2,
bip, atf5, atf6b, prdx4, ddit3, ero1l, gcg, brn, arx,
grehlin, sox9, hnf6, hes1, ngn3, drd2, rit2, ets1, plagl1, bhlhe40, hmgb3, asxl3,
esrrb, slca1, slca2, gck, hspa5, slc16a1, slc16a4, sst, ppy, neurog3,"""
""""Notes
- assuming esrrb1 is esrrb
- assurming glut1 is SLC2A1
- assuming glut2 is SLC2A2
- assuming hk4 is GCK
- assuming bip is HSPA5
- assuming mct1 is SLC16A1
- assuming mct4 is SLC16A4
- assuming kir6.2 is KCNJ11
- assuming glucagon is gcg
- assuming somatostatin is SST
- assuming pancreatic polypetide is PPY
- assuming ngn3 is NEUROG3
- I don't know what BRN is
"""
# genes = """
# ,Pcdha1,	Pcdhb1,	Pcdhga1,
# ,Pcdha2,	Pcdhb2,	Pcdhga2
# ,Pcdha3,	Pcdhb3,	Pcdhga3
# ,Pcdha4,	Pcdhb4,	Pcdhga4
# ,Pcdha5,	Pcdhb5,	Pcdhga5
# ,Pcdha6,	Pcdhb6,	Pcdhga6
# ,Pcdha7,	Pcdhb7,	Pcdhga7
# ,Pcdha8,	Pcdhb8,	Pcdhga8
# ,Pcdha9,	Pcdhb9,	Pcdhga9
# ,Pcdha10,	Pcdhb10,	Pcdhga10
# ,Pcdha11,	Pcdhb11,	Pcdhga11
# ,Pcdha12,	Pcdhb12,	Pcdhga12
# ,Pcdhac1,	Pcdhb13,	Pcdhgb1
# ,Pcdhac2,	Pcdhb14,	Pcdhgb2
# ,Pcdhb15,	Pcdhgb4,
# ,Pcdhb16,	Pcdhgb5,
# ,Pcdhb17,	Pcdhgb6,
# ,Pcdhb18,	Pcdhgb7,
# ,Pcdhb19,	Pcdhgb8,
# ,Pcdhb20,	Pcdhgc3,
# ,Pcdhb21,	Pcdhgc4,
# ,Pcdhb22,	Pcdhgc5,
# """
# genes = """Wfs1, DCN, Htr2c, Grp, Gpr101, Col5a1, Gpc3, Prss12, Ndst4, Calb1, Matn2, Rph3a, Loxl1, Plagl1, Coch, Itga7,
#             Iyd, Pvalb, Slc6a1, vamp1, Map4k3, Amigo1, Amigo2, Col15a1, Ccdc3, Lct, Trhr"""

# genes = """CCL2
# ,ACT2
# ,CCL21
# ,CCL28
# ,CXCL13
# ,CCR5
# ,Il1A
# ,IL4
# ,IL10
# ,IL17A
# ,IL18
# ,TNFSF13B
# ,CTLA4
# ,CD40LG
# ,CD5
# ,CD40
# ,CD1D
# ,B2M
# ,CD74
# ,CD72"""
genes = """
Tbr1,
Rasgrf2,
Pvrl3,
Cux2,
Rorb,
Plcxd2,
Thsd7a,
Kcnk2,
Cplx3,
Sulf2,
Foxp2,
Syt6,
Rprm,
Nr4a2,
Synpr,
Pcp4,
Gad1,
Pvalb,
Sst,
Htr3a,
Vip,
Reln,
Cck,
Npy,
Lhx6,
Calb2,
Pde1a,
Lphn2,
Kcnip2,
Rgs10,
Nov,
Cpne5,
Slc5a7,
Crh,
Pax6,
Cxcl14,
Gda,
Sema3e,
"""
genes = """Esrrb
,Nanog
,Tcl1
,Sox2
,Tet1
,Tbx3
,Dppa3
,Prdm14
,Dnmt3b
,Dnmt1
,Sdha
,Kdm6a
,Fbxo15
,Lin28a
,Alkbh5
,Mettl3
"""
#genes = "Dnmt1,Dnmt3b"
target_genes = [x.strip() for x in genes.split(",")]
debug = True
timeout = 60
probes = main(target_genes, min_probes=12, max_probes=24,
              timeout=timeout, debug=debug, organism='mouse')
dir_name = "passed_probes_iPSC4"
for gene in target_genes:
    try:
        os.mkdir(dir_name)
    except:
        pass
    passed = False
    for masking in [5,4,3]:
        probes = get_probes(gene, "mouse", masking)
        if len(probes) > 12:
            out_path = os.path.join(dir_name, gene + '_Mask={}_N={}.csv'.format(masking, len(probes)))
            print(gene, masking, len(probes))
            pd.DataFrame(probes).to_csv(out_path)
            passed = True
    if passed == False:
        print("Couldn't Find probes for {}".format(gene))

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
