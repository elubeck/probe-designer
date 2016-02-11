import os
import sys
os.chdir(os.path.split(os.getcwd())[0])

from progressbar import ProgressBar
from tempfile import NamedTemporaryFile
from collections import defaultdict
import random
import csv
from utils.misc import gc_count, n_probes
from oligoarray_designer import OligoarrayDesigner
from progressbar import ProgressBar
from intron_designer import IntronRetriever
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from probe_refiner import ProbeFilter
import dataset

intron_probes = defaultdict(list)
with open("temp/introns_10-18-15.csv", "r") as f_in:
    for n, line in enumerate(csv.reader(f_in)):
        intron_probes[line[0]].append(line[2])

passed = {k for k, v in intron_probes.items() if len(v) >= 24}

# First get used probes
intron_db = dataset.connect("sqlite:///db/intron_probes_4.db")
probe_db = intron_db['mouse']
used_probes = set([row['Name'] for row in probe_db.distinct("Name")])
passed = passed.union(used_probes)
o = OligoarrayDesigner(
    blast_db='/home/eric/blastdb/old_format/transcribed_mouse2')
chunksize = 12
probe_size = 35
chunks = []
intron_retriever = IntronRetriever(skip_probes=passed, reversed=True)
p_bar = ProgressBar(maxval=intron_retriever.tot_introns).start()
for n, gene in enumerate(intron_retriever):
    p_bar.update(n)
    gene_chunks = []
    for chunk in gene:
        gene_chunks.append(chunk)
        if n_probes(gene_chunks) > 400:
            break
    if chunk.id in used_probes: continue
    if n_probes(gene_chunks) < 50:
        continue
    chunks.append(gene_chunks)
    if len(chunks) == chunksize:
        with NamedTemporaryFile("w") as fasta_input:
            # Break records into 1000nt chunks
            for gene_chunks in chunks:
                # Reduce potential probeset size to drop design time
                while n_probes(gene_chunks, probe_size) > 400:
                    gene_chunks = random.sample(gene_chunks,
                                                len(gene_chunks) - 1)
                SeqIO.write(gene_chunks, fasta_input, 'fasta')  # write all sequences to file
            fasta_input.flush()
            res = o.run(fasta_input.name,
                        max_dist=1000,
                        min_length=probe_size,
                        max_length=probe_size,
                        max_tm=100,
                        cross_hyb_temp=72,
                        min_tm=74,
                        secondary_struct_temp=76,
                        max_oligos=400,
                        timeout=20,
                        blast_timeout=3)
        probe_db.insert_many(res)
        p_bar.update(n)
        chunks = []
p_bar.finish()

hits = (probes['target'] for probes in probe_db.distinct("target"))
p_match = ({probe['seq']
            for probe in probe_db.find(target=hit)} for hit in hits)
long_enough = [p_set for p_set in p_match if len(p_set) >= 24]

probe_filterer = ProbeFilter(db='gencode_tracks_reversed_introns+mRNA')


def filter_wrapper(gene):
    res = {r["seq"] for r in probes.find(target=gene)}
    f_probes = []
    if len(res) >= 24:
        iters = 0
        f_probes = []
        while iters < 5:
            max_off = 1 + iters * 2000
            f_probes = probe_filterer.run(res,
                                          gene,
                                          n_probes=48,
                                          max_off_target=max_off)
            if len(f_probes) >= 24: break
            iters += 1
        f_probes = [{'target': gene, 'seq': probe} for probe in f_probes]
    return f_probes


print("Filtering")
flat_probes = []
all_probes = {}
dbs = ['intron_probes_4']
for db in dbs:
    print(db)
    db_name = "sqlite:///db/{}_filtered.db".format(db)
    probe_db_filtered = dataset.connect(db_name)
    filtered_probe_table = probe_db_filtered['mouse']
    db_loc = "sqlite:///db/{}.db".format(db)
    intron_db = dataset.connect(db_loc)
    probes = intron_db['mouse']
    genes = [g['target'] for g in probes.distinct('target')]
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

    # get genes with >=24 probes
    tall_probes = dict()
    for gene_d in filtered_probe_table.distinct('target'):
        gene = gene_d['target']
        probes = list(
            {p['seq']
             for p in filtered_probe_table.find(target=gene)})
        if len(probes) >= 24:
            tall_probes[gene] = probes
    all_probes[db] = tall_probes

# ##### pick between cds or full tx
# from itertools import groupby
# p = [(k, v, db_n) for db_n, db in all_probes.items() for k, v in db.items()]

# all_p2 = {}
# for gene, vals in groupby(sorted(p), lambda x: x[0]):
#     p_dict = {val[2]: val[1] for val in vals}
#     if len(p_dict) == 1:
#         all_p2[gene] = list(p_dict.values())[0]
#     else:
#         if len(list(p_dict.values())[0]) == len(list(p_dict.values())[1]):
#             all_p2[gene] = p_dict["brain_mrna_probes_cds_only"]
#         else:
#             all_p2[gene] = max(p_dict.values(), key=lambda x: len(x))

# search for redundant nested sequences
import probe_refiner
p_set2 = probe_refiner.probe_set_refiner(all_p2)

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
