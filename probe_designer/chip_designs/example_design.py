import os
import sys

from probe_designer.mRNA_designer import RNARetriever2, design_step
from probe_designer.probe_refiner import ProbeFilter, probe_set_refiner

retriever = RNARetriever2()
filterer = ProbeFilter(db='gencode_tracks_reversed_introns+mRNA', copy_num='brain')

genes = ['Actg1', 'Actr3', 'Adnp', 'atf4', 'atg5', 'B2m', 'Bcl2', 'BDNF',
         'Casp8', 'Cdc42', 'Cebpa', 'Cebpb', 'chat', 'Cited2', 'csf1', 'CSf1r',
         'dbh', 'Ddc', 'ddit4', 'ddx4', 'Dhx9', 'drd2', 'Fcgr2b', 'Gba',
         'Hhex', 'Hif1a', 'Huwe1', 'Il1a', 'Il1b', 'Il6', 'Insl3', 'Jak3',
         'Jarid2', 'lamp2', 'lc3b', 'Mapk8', 'marco', 'Nfatc1', 'Nfkbia',
         'Numb', 'Ogt', 'Pag1', 'pdzl', 'Pknox1', 'pnmt', 'Ppia', 'Prkcb',
         'psmc4', 'rab1a', 'Rela', 'S100a4', 'S100a6', 'Sfpi1', 'Skil',
         'Smarca2', 'smurf1', 'Stat1', 'Th', 'Tial1', 'TNFa', 'Trim30a', 'ubc',
         'usp14', 'usp18']
passed = {}
for gene in genes:
    try:
        name = gene[0].upper() + gene[1:].lower()
        # Make best set of probes
        name, probes1, seq = design_step(name, cds_only=True, length=35, spacing=0,
                                        gc_target=0.55, gc_min=0.35, gc_max=0.75)
        # Filter crappy probes out
        probes = filterer.run(probes1, name, match_thresh=18, n_probes=48,
                                max_off_target=50, off_target_hits=6)
        passed[gene] = probes
        print(gene, len(probes1), len(probes),  len(''.join(seq)))
        if len(probes) < 24:
            name, probes, seq = design_step(name, cds_only=False, length=35, spacing=0,
                                            gc_target=0.55, gc_min=0.35, gc_max=0.75)
            # Filter crappy probes out
            probes2 = filterer.run(probes, name, match_thresh=18, n_probes=48,
                                  max_off_target=50, off_target_hits=6)

            print(gene, len(probes), len(probes2), len(''.join(seq)))
            if len(probes2) > len(passed[gene]):
                passed[gene] = probes
    except AssertionError:
        print("Couldn't find {}".format(name))

import arrow
import csv
with open("tim_genes_{}.csv".format(arrow.now().timestamp), 'w') as f_out:
    cv = csv.writer(f_out)
    for name, genes in passed.items():
        for gene in genes:
            cv.writerow([name, gene])


print("DONE")
