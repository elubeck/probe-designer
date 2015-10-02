from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from builtins import *

import dataset

from progressbar import ProgressBar
from intron_designer2 import ProbeFilter
import probe_designer.mRNA_designer


# Get all genes
# db = dataset.connect("sqlite:///db/refGene.db")
# mouse = db['mouse']
# genes = {gene['name2'] for gene in mouse.distinct('name2')}
import probe_designer.probe_refiner

genes = ['Pde1a', 'Adgrl2', 'Kcnip2', 'Rgs10', 'Nov', 'Cpne5', 'Slc5a7', 'Crh',
         'Pax6', 'Gda', 'Sema3e', 'Gfap', 'Mfge8', 'Gja1', 'Slco1c1', 'Hexb',
         'Mrc1', 'Lyve1', 'Lyz2', 'Itpr2', 'Rhob', 'Omg', 'Klk6', 'Cldn5',
         'Acta2', 'Sox2', 'Pax6', 'Dcx', 'Th', 'Nes', 'Tbr1', 'Wfs1', 'Dcn',
         'Htr2c', 'Grp', 'Gpr101', 'Col5a1', 'Gpc3', 'Prss12', 'Ndst4',
         'Calb1', 'Matn2', 'Rph3a', 'Loxl1', 'Plagl1', 'Coch', 'Itga7', 'Lyd',
         'Amigo2']

c_genes = ['Pde1a', 'Adgrl2', 'Kcnip2', 'Rgs10', 'Nov', 'Cpne5', 'Slc5a7',
           'Crh', 'Pax6', 'Gda', 'Sema3e', 'Gfap', 'Mfge8', 'Gja1', 'Slco1c1',
           'Hexb', 'Mrc1', 'Lyve1', 'Lyz2', 'Itpr2', 'Rhob', 'Omg', 'Klk6',
           'Cldn5', 'Acta2', 'Sox2', 'Pax6', 'Dcx', 'Th', 'Nes', 'Tbr1',
           'Loxl1', 'Calb1', 'Col5a1', 'Lyd', 'Amigo2']

dbs = ["brain_mrna_probes_cds_only", "brain_mrna_probes_full_tx"]
cds_only = [True, False]
for db, cds_only in zip(dbs, cds_only):
    # Design All Genes
    gene_lookup = probe_designer.mRNA_designer.batch_design2(genes,
                                              cds_only=cds_only,
                                              db_name=db)

probe_filterer = ProbeFilter(copy_num='brain',
                             db='gencode_tracks_reversed_introns+mRNA')


def filter_wrapper(gene):
    res = [r["seq"] for r in probes.find(target=gene)]
    f_probes = []
    if len(res) > 24:
        iters = 0
        f_probes = []
        while iters < 30:
            max_off = 1 + iters * 2000
            f_probes = probe_filterer.run(res, gene,
                                          n_probes=24,
                                          max_off_target=max_off)
            if len(f_probes) >= 24: break
            iters += 1
        f_probes = [{'target': gene, 'seq': probe} for probe in f_probes]
    return f_probes


print("FILTERING")
flat_probes = []
all_probes = {}
for db in dbs:
    db_name = "sqlite:///db/{}_filtered.db".format(db)
    probe_db_filtered = dataset.connect(db_name)
    filtered_probe_table = probe_db_filtered['mouse']
    db_loc = "sqlite:///db/{}.db".format(db)
    intron_db = dataset.connect(db_loc)
    probes = intron_db['unfiltered_probes']
    genes = [g['target'] for g in probes.distinct('target')
             if g['target'] in genes]
    used_probes = {
        p['target']
        for p in filtered_probe_table.distinct('target')
    }
    remaining_genes = set(genes).difference(used_probes)
    p_bar = ProgressBar(maxval=len(remaining_genes))
    from itertools import imap
    # p = Pool(processes=cpu_count())
    for n, f_probes in enumerate(imap(filter_wrapper, remaining_genes)):
        if any(f_probes):
            filtered_probe_table.insert_many(f_probes)
            flat_probes.append(f_probes)
        p_bar.update(n)
    p_bar.finish()

    # Get Genes with >=24 probes
    tall_probes = dict()
    for gene_d in filtered_probe_table.distinct('target'):
        gene = gene_d['target']
        if gene in genes:
            probes = list(
                {p['seq']
                 for p in filtered_probe_table.find(target=gene)})
            if len(probes) >= 24:
                tall_probes[gene] = probes
    all_probes[db] = tall_probes

from itertools import groupby

p = [(k, v, db_n)
     for db_n, db in all_probes.iteritems() for k, v in db.iteritems()]

all_p2 = {}
for gene, vals in groupby(sorted(p), lambda x: x[0]):
    p_dict = {val[2]: val[1] for val in vals}
    if len(p_dict) == 1:
        all_p2[gene] = p_dict.values()[0]
    else:
        if len(p_dict.values()[0]) == len(p_dict.values()[1]):
            all_p2[gene] = p_dict["brain_mrna_probes_cds_only"]
        else:
            all_p2[gene] = max(p_dict.values(), key=lambda x: len(x))

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
blast_hits2 = blaster2.local_blast_query(
    fasta_str,
    db='gencode_tracks_reversed_introns+mRNA')
blast_res2 = blaster2.parse_hits(blast_hits2, match_thresh=18, strand=1)

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

# CHANGED TO p_set22 from p_set2
p_set22 = {g: {p['seq'] for p in probes} for g, probes in grouped_probes2}

n_oligos = len([g for k, v in p_set2.items() for g in v])

#TEMP CODE, DELETE BLOCK


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement[c] for c in seq.upper())[::-1]


from random import sample
import csv
sub_set = {
    k: sample(probes, 24)
    for k, probes in p_set2.iteritems() if k in c_genes
}

used_bridges = []
with open('temp/9-14-15/barcode_c_primers.csv', 'r') as f_in:
    for n, line in enumerate(csv.reader(f_in)):
        if n == 0: continue
        used_bridges.append(line[1])
used_bridges += [reverse_complement(bridge) for bridge in used_bridges]

with open('temp/validated_bridges-1k_old.txt', 'r') as f_in:
    bridges = [line.strip('\n') for line in f_in]

available_bridges = set(bridges).difference(used_bridges)
to_use = sample(available_bridges, len(sub_set))

gene_bridge = {}
with open("temp/brain_probes_9-27-15.csv", "wb") as f_out:
    c = csv.writer(f_out)
    for bridge, (gene, probes) in zip(to_use, sub_set.items()):
        gene_bridge[gene] = bridge
        for n, probe in enumerate(sample(probes, 24)):
            c.writerow([gene, n, "{} TATA {}".format(
                reverse_complement(bridge), probe)])

from itertools import product
barcodes = list(product(range(5), range(5)))
gene_barcodes = dict(zip(gene_bridge.keys(), barcodes))

initiators = ["gCATTCTTTCTTgAggAgggCAgCAAACgggAAgAg",
              "AgCTCAgTCCATCCTCgTAAATCCTCATCAATCATC",
              "AAAgTCTAATCCgTCCCTgCCTCTATATCTCCACTC",
              "CACATTTACAgACCTCAACCTACCTCCAACTCTCAC",
              "CACTTCATATCACTCACTCCCAATCTCTATCTACCC"]

bridge_template = "{} TATA {}"

with open("temp/brain_barcodes_9-27-15.csv", 'wb') as f_out:
    fo = csv.writer(f_out)
    fo.writerow(['Gene', 'Hyb 1', 'Hyb 2'])
    for bridge, barcode in gene_barcodes.iteritems():
        fo.writerow([bridge, barcode[0], barcode[1]])

with open("temp/brain_bridge_9-27-15.csv", "wb") as f_out:
    c = csv.writer(f_out)
    for gene, bridge in gene_bridge.items():
        for n, init_n in enumerate(gene_barcodes[gene]):
            c.writerow(["{}-{}".format(gene, n + 1),
                        "{} TATA {}".format(initiators[init_n], bridge)])

#TEMP CODE, DELETE BLOCK

######### BLAT ALL PROBES ###############

import probe_designer.utils.align_probes
probe_designer.utils.align_probes.blat_pset(p_set2, 'brain_8-18-15_5.bed')
