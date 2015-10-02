# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from builtins import *

import dataset

from progressbar import ProgressBar
from probe_designer.intron_designer import ProbeFilter
import probe_designer.utils.align_probes
import probe_designer.mRNA_designer
import probe_designer.probe_refiner

high_priority = ['Alkbh5',
                 'Ash1l',
                 'Axin2',
                 'Bmi1',
                 'bmp4',
                 'bmpr1a',
                 'T',
                 'Cdh1',
                 'Col5a2',
                 'ctcf',
                 'dazl',
                 'Dnmt1',
                 'dnmt3a',
                 'Dnmt3b',
                 'Dnmt3L',
                 'Dppa2',
                 'Dppa3',
                 'Dppa4',
                 'EED',
                 'ehmt2',
                 'Esrrb',
                 'Ezh2',
                 'Fbxo15',
                 'Fgf4',
                 'FGF5',
                 'Fgfr2',
                 'Foxa2',
                 'FoxO1',
                 'Ehmt2',
                 'Gsk3b',
                 'hdac1',
                 'hdac2',
                 'Hes1',
                 'Id1',
                 'Id2',
                 'Jarid2',
                 'kdm3a',
                 'kdm3b',
                 'kdm4c',
                 'kdm4d',
                 'kdm6a',
                 'klf2',
                 'klf4',
                 'klf5',
                 'Lifr',
                 'Lin28a',
                 'Lin28b',
                 'Kdm1a',
                 'mael',
                 'METTL14',
                 'METTL3',
                 'Myc',
                 'Mycn',
                 'Nanog',
                 'Nr0b1',
                 'pax6',
                 'Pecam1',
                 'Prdm14',
                 'Rest',
                 'Sall4',
                 'Sdha',
                 'setdb1',
                 'Smad1',
                 'Smad4',
                 'Smad5',
                 'Socs3',
                 'Sox2',
                 'Sp1',
                 'Suz12',
                 'sycp3',
                 'Tbx3',
                 'Tcl1',
                 'Tet1',
                 'tet2',
                 'tet3',
                 'Thy1',
                 'Trim28',
                 'Utf1',
                 'Xist',
                 'Zfp42',
                 'Mecp2',
                 'MBD3', ]
low_priority = [
    'Hck',
    'Ifitm1',
    'Il6st',
    'Kit',
    'Lefty1',
    'Lefty2',
    'Lifr',
    'Nodal',
    'Nog',
    'Numb',
    'Pten',
    'Sfrp2',
    'Tdgf1',
    'Fgf4',
    'Fgf5',
    'Gdf3',
    'Lefty1',
    'Lefty2',
    'Nodal',
    'Brix1',
    'Cd9',
    'Diap2',
    'Ifitm2',
    'Igf2bp2',
    'Lin28a',
    'Podxl',
    'Rest',
    'Sema3a',
    'Tert',
    'Il6st',
    'Lifr',
    'Socs3',
    'Stat3',
    'Fgfr1',
    'Fgfr2',
    'Fgfr3',
    'Fgfr4',
    'Ptch1',
    'Ptchd2',
    'Smo',
    'Gli1',
    'Gli2',
    'Gli3',
    'Sufu',
    'Ncstn',
    'Notch1',
    'Notch2',
    'Notch3',
    'Notch4',
    'Psen1',
    'Psen2',
    'Psenen',
    'Rbpjl',
    'Hes1',
    'Hey1',
    'Acvr1',
    'Acvr1b',
    'Acvr1c',
    'Acvr2a',
    'Acvr2b',
    'Acvrl1',
    'Amhr2',
    'Bmpr1a',
    'Bmpr1b',
    'Bmpr2',
    'Eng',
    'Ltbp1',
    'Ltbp2',
    'Ltbp3',
    'Ltbp4',
    'Rgma',
    'Tgfbr1',
    'Tgfbr2',
    'Tgfbr3',
    'Tgfbrap1',
    'Crebbp',
    'E2f5',
    'Ep300',
    'Rbl1',
    'Rbl2',
    'Smad1',
    'Smad2',
    'Smad3',
    'Smad4',
    'Smad5',
    'Smad6',
    'Smad7',
    'Smad9',
    'Sp1',
    'Zeb2',
    'Fzd1',
    'Fzd2',
    'Fzd3',
    'Fzd4',
    'Fzd5',
    'Fzd6',
    'Fzd7',
    'Fzd8',
    'Fzd9',
    'Lrp5',
    'Lrp6',
    'Vangl2',
    'Bcl9',
    'Bcl9l',
    'Ctnnb1',
    'Lef1',
    'Nfat5',
    'Nfatc1',
    'Nfatc2',
    'Nfatc3',
    'Nfatc4',
    'Pygo2',
    'Tcf7l1',
    'Tcf7',
    'Tcf7l2',
    'Mecp2',
    'MBD3',
]

gene_list = {
    gene[0].upper() + gene[1:].lower()
    for pset in [high_priority, low_priority] for gene in pset
}

dbs = ['es_mrna_probe_cds_only', 'es_mrna_probe_full_tx']
cds_only = [True, False]
for db, cds_only in zip(, cds_only):
    # Design All Genes
    gene_lookup = probe_designer.mRNA_designer.batch_design2(gene_list,
                                              cds_only=cds_only,
                                              db_name=db)

probe_filterer = ProbeFilter()


def filter_wrapper(gene):
    res = [r["seq"] for r in probes.find(target=gene)]
    f_probes = []
    if len(res) > 24:
        iters = 0
        f_probes = []
        while iters < 5:
            max_off = 1 + iters * 2000
            f_probes = probe_filterer.run(res, gene,
                                          n_probes=48,
                                          max_off_target=max_off)
            if len(f_probes) >= 24: break
            iters += 1
        f_probes = [{'target': gene, 'seq': probe} for probe in f_probes]
    return f_probes


flat_probes = []
all_probes = {}
for db in :
    db_name = "sqlite:///db/{}_filtered.db".format(db)
    probe_db_filtered = dataset.connect(db_name)
    filtered_probe_table = probe_db_filtered['mouse']
    db_loc = "sqlite:///db/{}.db".format(db)
    intron_db = dataset.connect(db_loc)
    probes = intron_db['unfiltered_probes']
    genes = [g['target'] for g in probes.distinct('target')
             if g['target'] in gene_list]
    used_probes = {
        p['target']
        for p in filtered_probe_table.distinct('target')
    }
    remaining_genes = set(genes).difference(used_probes)
    p_bar = ProgressBar(maxval=len(remaining_genes))
    from multiprocessing import Pool, cpu_count
    p = Pool(processes=cpu_count())
    for n, f_probes in enumerate(p.imap(filter_wrapper, remaining_genes)):
        if any(f_probes):
            filtered_probe_table.insert_many(f_probes)
            flat_probes.append(f_probes)
        p_bar.update(n)
    p_bar.finish()

    # Get Genes with >=24 probes
    tall_probes = dict()
    for gene_d in filtered_probe_table.distinct('target'):
        gene = gene_d['target']
        probes = list(
            {p['seq']
             for p in filtered_probe_table.find(target=gene)})
        if len(probes) >= 24:
            tall_probes[gene] = probes
    all_probes[db] = tall_probes

##### Pick between CDS or full tx
from itertools import groupby
p = [(k, v, db_n)
     for db_n, db in all_probes.iteritems() for k, v in db.iteritems()]

zak_picks = [
    'Dazl', 'Dppa2', 'Fgf4', 'Fgf5', 'Kdm4d', 'Klf2', 'Klf4', 'Lin28a',
    'Lin28b', 'Mecp2', 'Mettl14', 'Sall4', 'Socs3', 'Sox2', 'Xist', 'Zfp42',
    'Foxa2', 'Hdac2'
]

all_p2 = {}
for gene, vals in groupby(sorted(p), lambda x: x[0]):
    p_dict = {val[2]: val[1] for val in vals}
    if gene in zak_picks:
        all_p2[gene] = p_dict['dbs_mrna_probe_full_tx']
    else:
        try:
            all_p2[gene] = p_dict['es_mrna_probe_cds_only']
        except:
            continue

print(len(all_p2))

# all_p2 = {}
# for gene, vals in groupby(sorted(p), lambda x: x[0]):
#     p_dict = {val[2]: val[1] for val in vals}
#     if len(p_dict) == 1:
#         all_p2[gene] = p_dict.values()[0]
#     else:
#         if len(p_dict.values()[0]) == len(p_dict.values()[1]):
#             all_p2[gene] = p_dict['es_mrna_probe_cds_only']
#         else:
#             all_p2[gene] = max(p_dict.values(), key=lambda x: len(x))

# Search for redundant nested sequences
p_set2 = probe_designer.probe_refiner.probe_set_refiner(all_p2)

p_set2 = {k: v for k, v in p_set2.items() if len(v) >= 24}

# Blast Everything to check that all probes together don't cause big problems
flat_probes = [
    {"gene": gene,
     "seq": probe,
     "number": n,
     "name": "{}-{}".format(gene, n)}
    for gene, probes in p_set2.iteritems() for n, probe in enumerate(probes)
]

import blaster2
fasta_str = "\n".join(">{name}\n{seq}".format(**probe)
                      for probe in flat_probes)
blast_hits2 = blaster2.local_blast_query(fasta_str,
                                         db='gencode_tracks_reversed')
blast_res2 = blaster2.parse_hits(blast_hits2, match_thresh=20, strand=1)

from collections import defaultdict
# For every gene hit in blast results see what probes hit it
hit_to_probe = defaultdict(list)
for probe, hits in blast_res2.iteritems():
    for hit in hits:
        hit_to_probe[hit].append(probe)

# Get probes that hit unintended target
bad_hits = [(hit, [probe for probe in probes
                   if probe.split('-')[0].lower() not in hit.lower()])
            for hit, probes in hit_to_probe.iteritems()]

bad_hits = sorted([(k, hits) for k, hits in bad_hits if any(hits)],
                  key=lambda x: len(x[1]),
                  reverse=True)

# Drop probes if more than thresh probes hit x-target or copy number is greater than c_thresh
thresh = 5
c_thresh = 50
drop_list = []
for hit, probes in bad_hits:
    ensembl_name = hit.split(',')[0]
    if len(probes) >= thresh or probe_filterer.get_copynum(
        [ensembl_name]) > c_thresh:
        for probe in probes:
            drop_list.append(probe)

drop_probes = set(drop_list)

# Make final list of good probes
flat_probes_filtered = [probe for probe in flat_probes
                        if probe['name'] not in drop_probes]

from itertools import groupby
# Count remaining probes
grouped_probes = [(group, list(probes))
                  for group, probes in groupby(flat_probes_filtered,
                                               key=lambda x: x['gene'])]
passed = {g for g, probes in grouped_probes if len(probes) >= 24}
grouped_probes2 = [(g, p) for g, p in grouped_probes if g in passed]
p_set2 = {g: {p['seq'] for p in probes} for g, probes in grouped_probes2}

import csv
with open("temp/zak_mrna_8-19-15_2.csv", 'wb') as f_out:
    c = csv.writer(f_out)
    c.writerow(['Gene', "# Probes"])
    for k, v in p_set2.iteritems():
        c.writerow([k, len(v)])

n_oligos = len([g for k, v in p_set2.items() for g in v])

with open("temp/zak_8-24-15.csv", "wb") as f_out:
    c = csv.writer(f_out)
    for gene, probes in p_set2.items():
        for n, probe in enumerate(probes):
            c.writerow([gene, n, probe])

######### BLAT ALL PROBES ###############

# try:
#     isinstance(p_set2, dict)
# except:
    # p_set2 = defaultdict(list)
    # with open("temp/zak_8-24-15.csv", "r") as f_out:
    #     c = csv.reader(f_out)
    #     for line in c:
    #         gene, n, probe = line
    #         p_set2[gene].append(probe)

import probe_designer.utils.align_probes
probe_designer.utils.align_probes.blat_pset(p_set2, 'zak_8-18-15_91.bed')
