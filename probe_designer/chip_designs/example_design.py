import os
import sys

from probe_designer.mRNA_designer import RNARetriever2, design_step
from probe_designer.probe_refiner import ProbeFilter, probe_set_refiner

retriever = RNARetriever2()
filterer = ProbeFilter(db='gencode_tracks_reversed_introns+mRNA', copy_num='brain')

# # Make best set of probes
# name, probes, sequence = design_step('Pgk1', cds_only=True, length=35,
#                                      spacing=0, gc_target=0.55, gc_min=None, gc_max=None)
# probes = set(probes)
# # Filter crappy probes out
# probes = filterer.run(probes, name, match_thresh=18, n_probes=48, max_off_target=2000)
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
        name, probes, seq = design_step(name, cds_only=True, length=26, spacing=0,
                                        gc_target=0.55, gc_min=0.35, gc_max=0.75)
        probes = filterer.run(probes, name, match_thresh=14, n_probes=48,
                                max_off_target=2000, off_target_hits=6)
        print(gene, len(probes), len(seq))
        passed[gene] = probes
    except AssertionError:
        print("Couldn't find {}".format(name))

print("BOB")
print("DONE")



 # This is useful if you need to refine a set of probes
 p_set2 = probe_set_refiner(all_p2, block_size=20) # Removes redundant probes of length block_size

 # Only keep probes with at least 24 probes per set
 p_set2 = {k: v for k, v in p_set2.items() if len(v) >= 24}

 # blast everything to check that all probes together don't cause big problems
 flat_probes = [
     {"gene": gene,
      "seq": probe,
      "number": n,
      "name": "{}-{}".format(gene, n)}
     for gene, probes in p_set2.items() for n, probe in enumerate(probes)
     ]

 import blaster
 fasta_str = "\n".join(">{name}\n{seq}".format(**probe)
                       for probe in flat_probes)
 blast_hits2 = blaster.local_blast_query(
     fasta_str,
     db='gencode_tracks_reversed_introns+mrna')
 blast_res2 = blaster.parse_hits(blast_hits2, match_thresh=18, strand=1)

 from collections import defaultdict
 # for every gene hit in blast results see what probes hit it
 hit_to_probe = defaultdict(list)
 for probe, hits in blast_res2.items():
     for hit in hits:
         hit_to_probe[hit].append(probe)

 # get probes that hit unintended target
 bad_hits = [(hit, [probe for probe in probes
                    if probe.split('-')[0].lower() not in hit.lower()])
             for hit, probes in hit_to_probe.items()]

 bad_hits = sorted([(k, hits) for k, hits in bad_hits if any(hits)],
                   key=lambda x: len(x[1]),
                   reverse=true)

 # drop probes if more than thresh probes hit x-target or copy number is greater than c_thresh
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

 # make final list of good probes
 flat_probes_filtered = [probe for probe in flat_probes
                         if probe['name'] not in drop_probes]

 from itertools import groupby
 # count remaining probes
 grouped_probes = [(group, list(probes))
                   for group, probes in groupby(flat_probes_filtered,
                                                key=lambda x: x['gene'])]
 passed = {g for g, probes in grouped_probes if len(probes) >= 24}
 grouped_probes2 = [(g, p) for g, p in grouped_probes if g in passed]
 p_set2 = {g: {p['seq'] for p in probes} for g, probes in grouped_probes2}

 n_oligos = len([g for k, v in p_set2.items() for g in v])

 with open("temp/intron_11-15.csv", "wb") as f_out:
     c = csv.writer(f_out)
     for gene, probes in p_set2.items():
         for n, probe in enumerate(probes):
             c.writerow([gene, n, probe])
