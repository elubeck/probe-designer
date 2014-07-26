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




def maximize_masking(probes):
    passed = defaultdict(dict)
    for k, v in groupby(probes, key=lambda x: x['Name']):
        for k1, v1 in groupby(v, key=lambda x: x['Masking']):
            v1 = list(v1)
            n_probes = sum([len(probe_set['Probes']) for probe_set in v1])
            if n_probes >= 24:
                passed[k][k1] = v1
    p2 = {}
    for gene, m_dist in passed.iteritems():
        highest_key = sorted(m_dist.keys(), reverse=True)[0]
        p2[gene] = m_dist[highest_key]
    return p2

def flatten_probes(p2):
    #Flatten out probe lists
    flat_probes = {}
    for gene, cds_psets in p2.iteritems():
        for probeset in cds_psets:    
            for k,v in probeset.iteritems():
                if any([k == 'Probes', k == 'Target Seq']):
                    continue
                probeset['Probes'][k] = v
        p_list = map(lambda x: x['Probes'], cds_psets)
        probes = pd.concat(p_list, ignore_index=True)
        flat_probes[gene] = probes

    #Give each probe a name
    for gene, probes in flat_probes.iteritems():
        name_tuple = zip(probes['Name'], probes['CDS Region #'].map(str), probes['Probe Position*'].map(str))
        probes['Probe Name'] = map('-'.join, name_tuple)
    return flat_probes


def main(target_genes, min_probes=24, timeout=120, debug=False):
    assert(isinstance(target_genes, list))
    passed_genes = {}
    for gene in target_genes:
        cds = get_seq.CDS(gene).run()
        if len(''.join(cds['CDS List'])) > (20*min_probes):
            passed_genes[gene] = cds
    biosearch = biosearch_designer.Biosearch()
    probes = biosearch.design(passed_genes, min_probes)
    biosearch.close()
    good_probes = maximize_masking(probes)
    flat_probes = flatten_probes(good_probes)
    blasted_probes = {}
    passed_probes = {}
    for gene, probe_df in flat_probes.iteritems():
        blasted_probes[gene] = blaster.blast_probes(gene, probe_df, timeout=timeout, debug=debug)
        passed_probes[gene] = blaster.filter_probes_based_on_blast(gene, blasted_probes[gene], probe_df,
        max_false_hits=8, debug=debug)
    return {'Blast': blasted_probes,
            'Passed': passed_probes}



doc = """
Usage: probe_designer.py TARGETS [-o=OUTPUT] [-m=MIN_PROBES] [-d=DEBUG] [-t=TIMEOUT]

Arguments:
    TARGETS Comma-separated list of genes to design probes towards

Options:
    -o --output Whether to write each passed probe to a csv file.
    -m --min_probes Minimum # of Probes to accept for a gene
    -d --debug Whether to print debug output
    -t --timeout How long(seconds) to wait for a response from NCBI.  Default=120 seconds.
"""

if __name__ == '__main__':
    args = docopt(doc,)
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
        timeout=int(args['--timeout'])
    else:
        timeout=120
    probes = main(target_genes, min_probes, timeout, debug)
    if args['--output'] is not False:
        for gene, probes in probes['Passed'].iteritems():
            try:
                os.mkdir('passed_probes')
            except:
                pass
            out_path = os.path.join('passed_probes', gene+'.csv')
            probes.to_csv(out_path)
    else:
        print(probes['Passed'])
