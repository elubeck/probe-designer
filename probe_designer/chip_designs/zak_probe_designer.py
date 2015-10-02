from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from builtins import *
import dataset
import intron_designer2
from progressbar import ProgressBar

import probe_designer.mRNA_designer

high_priority = ['Albhk5',
                 'Ash1',
                 'Axin2',
                 'Bmi1',
                 'bmp4',
                 'bmpr1a',
                 'Brachyury',
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
                 'G9A',
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
                 'LSD1',
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

zak_genes = {
    gene.lower()
    for pset in [high_priority, low_priority] for gene in pset
}

# probes = mRNA_designer.batch_design(zak_genes)

####### Filter New and old Probes  #############
dbs = ["sqlite:///db/mrna_probes.db"]
filterer = intron_designer2.ProbeFilter()
probe_db_filtered = dataset.connect("sqlite:///db/mRNA_probes_8-4-15.db")
filtered_probe_table = probe_db_filtered['mouse']


def filter_wrapper(gene):
    res = [r["seq"] for r in probes.find(target=gene)]
    f_probes = []
    if len(res) >= 20:
        f_probes = [{'target': gene,
                     'seq': p} for p in filterer.run(res, gene,
                                                     n_probes=48)]
    return f_probes

# Commented out to save runtime
for db in dbs:
    intron_db = dataset.connect(db)
    probes = intron_db['unfiltered_probes']
    genes = [g['target'] for g in probes.distinct('target')]
    used_probes = {
        p['target']
        for p in filtered_probe_table.distinct('target')
    }
    remaining_genes = set(genes).difference(used_probes)
    p_bar = ProgressBar(maxval=len(remaining_genes))
    from multiprocessing import Pool, cpu_count
    p = Pool(processes=cpu_count())
    from itertools import imap
    for n, f_probes in enumerate(p.imap(filter_wrapper, remaining_genes)):
        if any(f_probes):
            filtered_probe_table.insert_many(f_probes)
        p_bar.update(n)
    p_bar.finish()

# Get Genes with >=24 probes
all_probes = dict()
for gene_d in filtered_probe_table.distinct('target'):
    gene = gene_d['target']
    probes = [p['seq'] for p in filtered_probe_table.find(target=gene)]
    if len(probes) >= 24:
        all_probes[gene] = probes

# Search for redundant nested sequences
p_set2 = probe_designer.mRNA_designer.probe_set_refiner(all_probes)

p_set2 = {k: v for k, v in p_set2.items() if len(v) >= 24}

import blaster2

# Blast Everything to check that all probes together don't cause big problems
flat_probes = [
    {"gene": gene,
     "seq": probe,
     "number": n,
     "name": "{}-{}".format(gene, n)}
    for gene, probes in p_set2.iteritems() for n, probe in enumerate(probes)
]
# Commented out to save runtime
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
                   if probe.split('-')[0] not in hit.lower()])
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
    if len(probes) >= thresh or filterer.get_copynum(
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
passed = [g for g, probes in grouped_probes if len(probes) >= 40]

filterer = intron_designer2.ProbeFilter()

######### BLAT ALL PROBES ###############
import dataset
import probe_designer.align_probes
from progressbar import ProgressBar
from multiprocessing import Pool, cpu_count
probe_alignments = dataset.connect("sqlite:///db/mrna_probes_aligned.db")
mouse_alignments = probe_alignments['mouse']


def blat_wrapper(d):
    try:
        gene, passed = d
        gene = gene[0].upper() + gene[1:]
        if not passed:
            return {}
        chr = probe_designer.align_probes.find_chromosome(gene)
        chr_positions = probe_designer.align_probes.blat_probes(passed, chr)
        for pos in chr_positions:
            pos.update({'target': gene})
        return chr_positions
    except:
        return {}


n_genes = len(p_set2)
p_bar = ProgressBar(maxval=n_genes)
p = Pool(processes=cpu_count())
from itertools import imap
for n, chr_positions in enumerate(imap(blat_wrapper, p_set2.items())):
    mouse_alignments.insert_many(chr_positions)
    p_bar.update(n)

with open('temp/zak_mrna_8-4-15.bed', 'wb') as f_out:
    bed = csv.writer(f_out, delimiter='\t')
    for row in mouse_alignments:
        bed.writerow([row['tName'], row['tStart'], row['tEnd'], row['strand'],
                      row['target'] + '-' + row['qName']])
