"""
This file was used to design probes for 10k genes. 
Both the old database from the previous intron design project and a new database designing more liberal intron probes were used.  These databases were generated on whole genome runs using the intron_designer2.design_introns() function. 
All subsequent steps are documented in this file.
"""

from __future__ import print_function, division
import dataset
import probe_designer.intron_designer2
from collections import defaultdict
from itertools import groupby
import probe_designer.blaster2
from progressbar import ProgressBar
import probe_designer.mRNA_designer
import probe_designer.blaster2
import csv
import json
from probe_designer.get_seq import reverse_complement
from collections import Counter

filterer = probe_designer.intron_designer2.ProbeFilter()

####### Filter Probes  #############
dbs = ["sqlite:///db/intron_probes_10k_2.db"]

# Open Filtered Probe Database
filterer = probe_designer.intron_designer2.ProbeFilter()
intron_db_filtered = dataset.connect(
    "sqlite:///db/intron_probes_filtered_10k_2.db")
filtered_probe_table = intron_db_filtered['mouse']

for db in dbs:
    intron_db = dataset.connect(db)
    probes = intron_db['mouse']
    genes = [g['Name'] for g in probes.distinct('Name')]
    used_probes = {
        p['target']
        for p in filtered_probe_table.distinct('target')
    }
    remaining_genes = set(genes).difference(used_probes)
    p_bar = ProgressBar(maxval=len(remaining_genes))
    for n, gene in enumerate(remaining_genes):
        if gene in used_probes: continue
        res = [r["Probe (5'-> 3')"] for r in probes.find(Name=gene)]
        if len(res) >= 20:
            f_probes = [{'target': gene,
                         'seq': p} for p in filterer.run(res, gene,
                                                         n_probes=25)]
            filtered_probe_table.insert_many(f_probes)
        p_bar.update(n)
    p_bar.finish()

# Get Genes with >=24 probes
all_probes = {}
for gene_d in filtered_probe_table.distinct('target'):
    gene = gene_d['target']
    probes = [p['seq'] for p in filtered_probe_table.find(target=gene)]
    if len(probes) >= 24:
        all_probes[gene] = probes