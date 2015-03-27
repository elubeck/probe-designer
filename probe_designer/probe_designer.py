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
    #Choose first probe of every repeat.  Because multiple maskings are allowed, same probe can be chosen many times
    pre = len(df)
    unique_index = [df.index[df['Probe Name']==name][0] for name in df['Probe Name'].unique()]
    from random import shuffle
    shuffle(unique_index)
    df = df.ix[unique_index] #Do this in a random order to reduce # of small blast runs
    print(pre, len(df))
    #Try splitting probelist into smaller chunks probe sets
    chunk_size = 2000
    b_res = {}
    print("Blasting: {} total probes".format(len(df)))
    for i in range(0, len(df), chunk_size):
        sub_df = df.iloc[i:i+chunk_size]
        print("Blasting: {} to {}".format(i, i+chunk_size-1))
        blaster.blast_probes('', sub_df, timeout=1800, debug=False, organism=cds_org)
        # sub_bres = []
        # b_res.update(sub_bres)
    print("Getting blast results from DB")
    #b_res = blaster.check_db(df, organism=cds_org)
    #b_res = blaster.blast_probes('', df, timeout=1800, debug=False, organism=cds_org)
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
            g_set[gene][masking] = passed
            for n, line in passed.T.to_dict().items():
                p_table.insert(line)
            print("Probeset {}, {} added to table".format(gene, masking))
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
,Drd1a
,Serpina9
,Ngfr
,Ntrk1
,Rarb
,Nts
,Pcdh11x
,Slc6a3
,Magell2
,Lyd,
Lm1,
Tiam2,
Robo1,
Cadps2
"""

genes = [u'Smug1',
 u'Npr1',
 u'Axin2',
 u'Kcnh3',
 u'Apbb2',
 u'Glp1r',
 u'Ppfibp1',
 u'Ano2',
 u'Foxa1',
 u'Maob',
 u'Wif1',
 u'Rnf144b',
 u'Snca',
 u'Tdo2',
 u'Slc17a8',
 u'Neto2',
 u'Pde5a',
 u'Fezf1',
 u'Cdh8',
 u'Anxa2',
 u'Esr2',
 u'Anxa4',
 u'Grid2',
 u'Sun2',
 u'Itgb8',
 u'Nenf',
 u'Kcng4',
 u'Syt6',
 u'A230065H16Rik',
 u'Rgs4',
 u'Kcnmb4',
 u'Gpx3',
 u'Stk24',
 u'Dap',
 u'Prss35',
 u'Dach1',
 u'Lrmp',
 u'Nov',
 u'Clca1',
 u'Fam40b',
 u'Iyd',
 u'Nfia',
 u'Ccnd2',
 u'Gabrr1',
 u'Thbs2',
 u'Tph2',
 u'Serpina3k',
 u'Dlk1',
 u'Rrm2',
 u'9630033F20Rik',
 u'Adam19',
 u'Cartpt',
 u'Lhx1',
 u'Slc26a10',
 u'Crim1',
 u'3110047P20Rik',
 u'F5',
 u'Slc6a2',
 u'Ptpru',
 u'Zdhhc2',
 u'Arl10',
 u'Serpina9',
 u'Kcnq4',
 u'Prox1',
 u'Kitl',
 u'Agrp',
 u'Prokr2',
 u'Arc',
 u'2310042E22Rik',
 u'Col23a1',
 u'Gfra1',
 u'Pstpip1',
 u'Sipa1l2',
 u'Cacng5',
 u'Trim16',
 u'Fam163b',
 u'Kcnj5',
 u'Zfp114',
 u'Ism1',
 u'Slc41a3',
 u'LOC433088',
 u'Npnt',
 u'Pcdhac1',
 u'Calb1',
 u'Prkg2',
 u'Slc18a3',
 u'Cdh11',
 u'Lgals1',
 u'Scn5a',
 u'Sema7a',
 u'Camk1d',
 u'Foxp2',
 u'Kcnc3',
 u'Doc2b',
 u'Pacrg',
 u'B3galt2',
 u'Gpr155',
 u'Fabp7',
 u'Irs4',
 u'Nts',
 u'Mlec',
 u'Gabrq',
 u'Rai14',
 u'Ankrd34b',
 u'Cplx3',
 u'Crym',
 u'Nhlh2',
 u'Baiap3',
 u'Pln',
 u'Fnbp1l',
 u'Rreb1',
 u'Fgd5',
 u'Spint2',
 u'Tmed3',
 u'Fam65b',
 u'Isoc1',
 u'Slc10a4',
 u'Layn',
 u'Nqo2',
 u'Cdh23',
 u'Cdh24',
 u'Ly6h',
 u'Lxn',
 u'Bok',
 u'Cyp39a1',
 u'Hs6st3',
 u'Ltk',
 u'Itgb5',
 u'Panx2',
 u'Scai',
 u'Chst8',
 u'Ngef',
 u'BC034076',
 u'Hspb1',
 u'Cables2',
 u'Bhlhe41',
 u'Lypd6b',
 u'Gng12',
 u'Col9a1',
 u'Clptm1l',
 u'Me2',
 u'Slc16a2',
 u'Smoc1',
 u'Sv2b',
 u'Rftn1',
 u'Mki67',
 u'Tox',
 u'Vav2',
 u'Entpd3',
 u'Ntng1',
 u'Acvr1c',
 u'Ddit4l',
 u'Slc17a6',
 u'Slc17a7',
 u'AF529169',
 u'Tcn2',
 u'Ctgf',
 u'Ust',
 u'Fam129a',
 u'Gpc5',
 u'Tmem215',
 u'Trpv4',
 u'Chrna5',
 u'Gchfr',
 u'Etv1',
 u'Htr5b',
 u'Clic5',
 u'Lats2',
 u'Ntn1',
 u'Chrnb3',
 u'Prss12',
 u'Lrp1b',
 u'Zfhx3',
 u'Ntsr1',
 u'LOC433228',
 u'Stac',
 u'Gpr133',
 u'Moxd1',
 u'Ghsr',
 u'Chrna6',
 u'Sh3kbp1',
 u'Tspan18',
 u'Lipa',
 u'Loxl2',
 u'Gucy2c',
 u'Bcl2l11',
 u'LOC381076',
 u'Rcn1',
 u'Glra3',
 u'Tspan11',
 u'Glra1',
 u'1810041L15Rik',
 u'Rorb',
 u'Sulf1',
 u'Esyt3',
 u'Clca2',
 u'Sln',
 u'Chrna7',
 u'Slit3',
 u'Arhgap25',
 u'Chrna3',
 u'Cacna1g',
 u'Fam5c',
 u'Pitpnc1',
 u'Cda',
 u'Fam102b',
 u'Kcng3',
 u'Tac2',
 u'Vsnl1',
 u'Kcnj14',
 u'Cmklr1',
 u'Tanc1',
 u'Gpr101',
 u'A930038C07Rik',
 u'Adcyap1',
 u'Tpd52l1',
 u'Ptpn3',
 u'A830036E02Rik',
 u'Amigo2',
 u'Ifit3',
 u'Pcdh11x',
 u'9430028L06Rik',
 u'Ttc6',
 u'Allc',
 u'Iqcj',
 u'Bace1',
 u'Itga7',
 u'Gnb4',
 u'Prkcd',
 u'Slc22a3',
 u'Id2',
 u'Fam46a',
 u'Zfp618',
 u'Ddc',
 u'Icam5',
 u'Gsg1l',
 u'Acaa2',
 u'Map3k15',
 u'Th',
 u'LOC432928',
 u'Dusp5',
 u'Abcd2',
 u'Nrp1',
 u'Rgs16',
 u'Fjx1',
 u'Enox2',
 u'Phactr2',
 u'Fign',
 u'Atp10a',
 u'Limk1',
 u'Trim36',
 u'Slc8a2',
 u'Kit',
 u'Cdh9',
 u'Ccdc3',
 u'Npy',
 u'Rab3b',
 u'Pomc',
 u'Gpr50',
 u'Myo16',
 u'Nxph4',
 u'Plxdc2',
 u'Chrna2',
 u'Gphn',
 u'Prcp',
 u'1600021P15Rik',
 u'1700019N12Rik',
 u'Ace',
 u'Lypd1',
 u'Prlr',
 u'Gabrg1',
 u'Ndst4',
 u'Camk4',
 u'Prkcg',
 u'Ntrk1',
 u'Rps5',
 u'Lpl',
 u'Cyld',
 u'Prkcq',
 u'Fam19a2',
 u'C230071H18Rik',
 u'Gls',
 u'Rasl11b',
 u'Calcb',
 u'Lhfpl2',
 u'Coch',
 u'Rab37',
 u'Acan',
 u'Sema5b',
 u'Sntb1',
 u'Zfpm2',
 u'1110007C09Rik',
 u'Qrfpr',
 u'Tiam2',
 u'Ramp3',
 u'Grm1',
 u'Fam69c',
 u'Gda',
 u'Syt10',
 u'Pax6',
 u'Prune2',
 u'Satb2',
 u'Nrgn',
 u'Npy1r',
 u'Vipr2',
 u'Dlx1',
 u'Dpp6']

#target_genes = [x.strip() for x in genes.split(",")]
target_genes = genes
debug = True
timeout = 60
#probes = main(target_genes, min_probes=16, max_probes=64,
#              timeout=timeout, debug=debug, organism='mouse')
dir_name = "passed_probes_wholeBrain_ABA_fine_structure_top_5"
n_probes = 0
for gene in target_genes:
    try:
        os.mkdir(dir_name)
    except:
        pass
    passed = False
    for masking in [5,4,3]:
        probes = get_probes(gene, "mouse", masking)
        if len(probes) > 16:
            out_path = os.path.join(dir_name, gene + '_Mask={}_N={}.csv'.format(masking, len(probes)))
            print(gene, masking, len(probes))
            pd.DataFrame(probes).to_csv(out_path)
            passed = True
    if passed:
        n_probes += 1
    if passed == False:
        print("Couldn't Find probes for {}".format(gene))
print(n_probes)

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
