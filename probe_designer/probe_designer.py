#!/usr/bin/env python
# -*- coding: utf-8 -*-
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


def maximize_masking(probes, max_probes=24):
    passed = defaultdict(dict)
    for k, v in groupby(probes, key=lambda x: x['Name']):
        for k1, v1 in groupby(v, key=lambda x: x['Masking']):
            v1 = list(v1)
            n_probes = sum([len(probe_set['Probes']) for probe_set in v1])
            if n_probes >= max_probes:
                passed[k][k1] = v1
    p2 = {}
    for gene, m_dist in passed.iteritems():
        highest_key = sorted(m_dist.keys(), reverse=True)[0]
        p2[gene] = m_dist[highest_key]
    return p2


def flatten_probes(p2):
    # Flatten out probe lists
    flat_probes = {}
    for gene, cds_psets in p2.iteritems():
        for probeset in cds_psets:
            for k, v in probeset.iteritems():
                if any([k == 'Probes', k == 'Target Seq']):
                    continue
                probeset['Probes'][k] = v
        p_list = map(lambda x: x['Probes'], cds_psets)
        probes = pd.concat(p_list, ignore_index=True)
        flat_probes[gene] = probes

    # Give each probe a name
    for gene, probes in flat_probes.iteritems():
        name_tuple = zip(probes['Name'], probes['CDS Region #'].map(str), probes['Probe Position*'].map(str))
        probes['Probe Name'] = map('-'.join, name_tuple)
    return flat_probes


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
        except:
            print("No probes could be designed towards CDS of %s" % gene)
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
            for gene, cds in passed_genes.iteritems():
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
    # good_probes = maximize_masking(probes, max_probes)
    # flat_probes = {gene:probes for gene, probes in flatten_probes(p_set).iteritems()
    # if len(probes)>min_probes}
    blasted_probes = {}
    passed_probes = {}
    g_set = defaultdict(dict)
    print("Breakpoint")
    import dataset
    db = dataset.connect("sqlite:///db/probes.db")
    p_table = db[organism]
    for gene, masked_probes in groupby(probes, key=lambda x: x['Name']):
        for p_set in sorted(masked_probes, key=lambda x: x['Masking'], reverse=True):
            probe_df = p_set['Probes']
            probe_df['Name'] = gene
            probe_df['CDS Region #'] = p_set['CDS Region #']
            if len(probe_df) < min_probes:
                continue
            name_tuple = zip(probe_df['Name'], probe_df['CDS Region #'].map(str), probe_df['Probe Position*'].map(str))
            probe_df['Probe Name'] = map('-'.join, name_tuple)
            try:
                b_res = blaster.blast_probes(gene, probe_df, timeout=timeout, debug=debug, organism=cds_org)
            except:
                print("Blast Failed for {}".format(gene))
                continue
            passed = blaster.filter_probes_based_on_blast(gene, b_res, probe_df, max_probes=max_probes,
                                                          min_probes=min_probes, debug=debug,
                                                          max_false_hits=min_probes / 2)
            g_set[gene][p_set['Masking']] = passed
            passed['Masking'] = p_set['Masking']
            for n, line in passed.T.to_dict().iteritems():
                p_table.insert(line)
            print("probeset added to table")
            print("Done Lvl")
    return g_set
    # return {'Blast': blasted_probes,
    #         'Passed': passed_probes}

    #
    # for gene, probe_df in flat_probes.iteritems():
    # try:
    # blasted_probes[gene] = blaster.blast_probes(gene, probe_df, timeout=timeout, debug=debug, organism=cds_org)
    # passed_probes[gene] = blaster.filter_probes_based_on_blast(gene, blasted_probes[gene], probe_df,
    # max_false_hits=8, debug=debug, max_probes=24,
    # min_probes=16)
    # except:
    # print("Failed to blast probes for %s" %gene)


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
genes = """grin1, grin2a, atp2b2, crmp1, bbc3, adora1, agrn, apoe, c19orf20, calb1, clock,
        dlg4, kifc2, mapk11, mib2, OXR1, pick1, pctk3, pyroxd1, aifm3, bad, bcl2, bok, app,
        caskin2, grin2c, homer1, homer3, kcnip3, nlgn2, nrgn, nrxn3, amh, c4orf48, rest,
        met, pten"""

genes = """nkx6-1
        ,pdx1
        ,mafb
        ,pax6
        ,ptf1a
        ,sox9
        ,ins
        ,gcg
        """
target_genes = [x.strip() for x in genes.strip(',').split(",")]
min_probes = 24
debug = True
timeout = 180
probes = main(target_genes, min_probes=12, max_probes=24,
              timeout=timeout, debug=debug, organism='human')
for gene, probes in probes['Passed'].iteritems():
    try:
        os.mkdir('passed_probes2p1')
    except:
        pass
    out_path = os.path.join('passed_probes2p1', gene + '.csv')
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
        for gene, probes in probes['Passed'].iteritems():
            try:
                os.mkdir('passed_probes')
            except:
                pass
            out_path = os.path.join('passed_probes', gene + '.csv')
            probes.to_csv(out_path)
    else:
        print(probes['Passed'])
