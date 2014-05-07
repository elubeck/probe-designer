#!/usr/bin/env python
# -*- coding: utf-8 -*-
import get_seq
import biosearch_designer
from itertools import groupby
import pandas as pd
import blaster

min_probes = 24


def maximize_masking(probes):
    passed = []
    for k, v in groupby(probes, key=lambda x: x['Name']):
        for k1, v1 in groupby(v, key=lambda x: x['Masking']):
            v1 = list(v1)
            n_probes = sum([len(probe_set['Probes']) for probe_set in v1])
            if n_probes >= 24:
                passed.append((k, k1, v1))
    p2 = {}
    for k, v in groupby(passed, key=lambda x: x[0]):
        p = [i for mask, i 
             in sorted(groupby(v, key=lambda x: x[1]), reverse=True)]
        p2[k] = list(p[0])[0][2]
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


def main(target_genes):
    assert(isinstance(target_genes, list))
    passed_genes = {}
    for gene in target_genes:
        handle = get_seq.get_mRNA_cdss(gene)
        cds = get_seq.merge_cds(handle)
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
        blasted_probes[gene] = blaster.blast_probes(gene, probe_df)
        passed_probes[gene] = blaster.filter_probes_based_on_blast(gene, blasted_probes[gene], probe_df, max_false_hits=8)
    return {'Blast': blasted_probes,
            'Passed': passed_probes}

            
if __name__ == '__main__':
    probes = main(['PNPLA5', 'Adamts8'])
