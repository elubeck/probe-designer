#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division, print_function
from future.builtins import str
from future.builtins import zip
from future.builtins import map
from past.utils import old_div
from itertools import groupby
from collections import defaultdict
from docopt import docopt
import os
import pandas as pd
import get_seq
import biosearch_designer
import blaster
import arrow
from pathlib import Path
import dataset

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


def blast_probes(cds_org, debug, max_probes, min_probes, organism, probes, timeout):
    g_set = defaultdict(dict)
    db = dataset.connect("sqlite:///db/probes.db")
    p_table = db[organism]
    for gene, masked_probes in groupby(probes, key=lambda x: x['Name']):
        for p_set in sorted(masked_probes, key=lambda x: x['Masking'], reverse=True):
            probe_df = p_set['Probes']
            probe_df['Name'] = gene
            probe_df['CDS Region #'] = p_set['CDS Region #']
            if len(probe_df) < min_probes:
                continue
            name_tuple = list(
                zip(probe_df['Name'], probe_df['CDS Region #'].map(str), probe_df['Probe Position*'].map(str)))
            probe_df['Probe Name'] = list(map('-'.join, name_tuple))
            try:
                b_res = blaster.blast_probes(gene, probe_df, timeout=timeout, debug=debug, organism=cds_org)
            except:
                print("Blast Failed for {}".format(gene))
                continue
            passed = blaster.filter_probes_based_on_blast(gene, b_res, probe_df, max_probes=max_probes,
                                                          min_probes=min_probes, debug=debug,
                                                          max_false_hits=old_div(min_probes, 2))
            g_set[gene][p_set['Masking']] = passed
            passed['Masking'] = p_set['Masking']
            for n, line in passed.T.to_dict().items():
                p_table.insert(line)
            print("Probeset added to table")
    return g_set


def main(target_genes, max_probes=24, min_probes=24, timeout=120, debug=False, parallel=True, organism="mouse"):
    if "mouse" == organism.lower():
        cds_org = '"Mus musculus"[porgn:__txid10090]'
        b_org = organism.lower()
    elif "human" == organism.lower():
        cds_org = '"Homo sapiens"[porgn:__txid9606]'
        b_org = organism.lower()
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
    probes = b_designer.design(passed_genes, min_probes,)
    b_designer.close()
    if debug:
        write_folder = Path("debug").joinpath("bsearch")
        try:
            write_folder.mkdir(parents=True)
        except:
            pass
        write_path = write_folder.joinpath("{}.csv".format(arrow.now()))
        pd.DataFrame(probes).to_csv(str(write_path))
    return blast_probes(cds_org, debug, max_probes, min_probes, organism, probes, timeout)


def get_probes(gene, organism, masking):
    db = dataset.connect("sqlite:///db/probes.db")
    p_table = db[organism]
    return p_table.find(Name=gene, masking=masking)

def get_optimal_probes(gene, organism, min_probes=12, max_probes=24):
    stored_probes= None
    for masking in (5,4,3):
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
# genes = """nrn1,s100a10,rbpms,fbxo2,nefl,tppp3,eg435376,sncg,rbpms2,rlbp1l2,slc17a6,chrnb3,bc089491,cpne9,tubb3,opn4,
# nppb,adcyap1,1700011l03Rik,prph,ctxn3,cd24a,prkcq,cdkn1c,gm687,trhde,igf1,rbms3,gm527,ndp,cyp1b1,rlbp1,
# kcnj13,rdh5,rpe65,rgr,efemp1,bbs2,bbs4"""
# genes = """grin1, grin2a, atp2b2, crmp1, bbc3, adora1, agrn, apoe, c19orf20, calb1, clock,
#         dlg4, kifc2, mapk11, mib2, OXR1, pick1, pctk3, pyroxd1, aifm3, bad, bcl2, bok, app,
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
genes = """ins, pdx1, mafa, mafb, nkx6.1, nkx2.2, islet1, pax6, pax4, nd1,
pcsk1, sur1, ldha, mvk, hk1, hk2,
bip, atf5, atf6b, prdx4, ddit3, ero1l, gcg, brn, arx,
grehlin, sox9, hnf6, hes1, ngn3, drd2, rit2, ets1, plagl1, bhlhe40, hmgb3, asxl3,
esrrb, slca1, slca2, gck, hspa5, slc16a1, slc16a4, sst, ppy, neurog3"""
"""Notes
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

target_genes = [x.strip() for x in genes.split(",")]
min_probes = 24
debug = True
timeout = 180
probes = main(target_genes, min_probes=12, max_probes=24,
              timeout=timeout, debug=debug, organism='human')
for gene, probes in probes['Passed'].items():
    try:
        os.mkdir('passed_probes3p3')
    except:
        pass
    out_path = os.path.join('passed_probes3p3', gene + '.csv')
    probes.to_csv(out_path)

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
