import dataset
import csv
from operator import itemgetter
from itertools import groupby
from joblib import Parallel, delayed

db = dataset.connect("sqlite:///refGene.db")
mouse = db['mouse']
genes = mouse.distinct('name2')


def get_range(sequence, min_size=20):
    #Find longest contiguous sequences
    ranges = []
    for k, g in groupby(enumerate(sequence), lambda (i, x): i - x):
        group = map(itemgetter(1), g)
        start = group[0]
        end = group[-1]
        if end - start > min_size:
            ranges.append((start, end))
    return ranges


def get_intron(targets, exons):
    intron_set = []
    # Pull out intron1
    for gene in targets:
        if gene['strand'] == '+':
            start = int(gene['txStart'])
            end = int(gene['txEnd'])
            e_start = map(int, gene['exonStarts'].strip(',').split(","))
            e_end = map(int, gene['exonEnds'].strip(',').split(","))
            # gene['exon_pos'] = {p for start, end in zip(e_start, e_end)
            #                     for p in range(start, end + 1)}
        elif gene['strand'] == '-':
            # Flip things around
            start = int(gene['txEnd'])
            end = int(gene['txStart'])
            e_end = map(int, gene['exonStarts'].strip(',').split(","))[::-1]
            e_start = map(int, gene['exonEnds'].strip(',').split(","))[::-1]
            # gene['exon_pos'] = {p for start, end in zip(e_end, e_start)
            #                     for p in range(start, end + 1)}
        if e_start[0] != start:
            # if gap between start and 1st exon
            intron1 = (start, e_start[0] - 1)
        elif len(e_start) != 1:
            # if multiple exons
            intron1 = (e_end[0] + 1, e_start[1] - 1)
        else:
            # If only 1 exon
            intron1 = (e_end[0], end)
        intron1 = sorted(intron1)  # Need to flip tuple around if '-' strand
        gene['intron1'] = intron1
        gene['intron1_len'] = intron1[1] - intron1[0]
        # Create a list of used positions
        intron_set.append(gene)

    # Make set of points that are only intron1, not any other exons
    # exon_pos = {p for gene in intron_set for p in gene['exon_pos']}
    exon_pos = exons[gene['chrom']]  # Get all the exons for chromosome
    sorted(intron_set, key=lambda x: x['intron1_len'])
    n_range = []
    for intronp1 in intron_set:
        intron = intronp1['intron1']
        intron1_range = set(range(intron[0], intron[1] + 1))
        non_overlapping = sorted(intron1_range.difference(exon_pos))
        final_range = get_range(non_overlapping)
        intron_len = sum([end - start for start, end in final_range])
        n_range.append((intron_len, final_range))
    # Make the longest intron first
    n_range = sorted(n_range, key=lambda x: x[0], reverse=True)
    if n_range[0][0] == 0 and len(e_end) != 1:
        print(gene['name2'], len(e_end))
    return (n_range[0][1], gene['strand'], gene['chrom'])


def wrapper(vals):
    return (vals[0], get_intron(vals[1]))


def get_exon_pos(table):
    from collections import defaultdict
    exon_space = defaultdict(set)
    for line in table:
        chrom = line['chrom']
        e_start = map(int, line['exonStarts'].strip(',').split(","))
        e_end = map(int, line['exonEnds'].strip(',').split(","))
        for start, end in zip(e_start, e_end):
            for i in range(start, end + 1):
                exon_space[chrom].add(i)
    return exon_space

# # THIS CODE GETS COORDINATES
# db2 = dataset.connect("sqlite:///intron1_coords.db")
# table2 = db2['mouse']
# from pathlib import Path
# from Bio import SeqIO
# chrom_root = Path("/Users/Eric/Downloads/chromFaMasked/")
# for chromosome in table2.distinct('chrom'):
#     chrom_code = chromosome['chrom']
#     chrom_path = chrom_root.joinpath("{}.fa.masked".format(chrom_code))
#     chrom_seq = SeqIO.read(str(chrom_path), 'fasta')
#     for row in table2.find(chrom=chrom_code):
#         coords_str = row['coords'].strip("[]")
#         # Intron doesn't have assigned coords
#         if not coords_str:
#             continue
#         for coord_pair in coords_str.split(", ("):
#             start, end = map(int, coord_pair.strip("()").split(", "))
#             seq = chrom_seq[start:end+1]
#             if row['strand'] == "-":
#                 seq = seq.reverse_complement()
#             # Drop masked sequences
#             sub_seqs = [sub_seq for sub_seq in seq.seq.split("N") if len(sub_seq) >= 20]
#             # TODO: Spit out to FASTA FILE.  Design

# Make a CSV file containing the full coordinates to every gene in the genome
# with open("gene_pos.csv", 'w') as f:
#     field_names = ['strand', 'start', 'end', 'chrom']
#     csv_f = csv.DictWriter(f, fieldnames=field_names)
#     csv_f.writeheader()
#     for gene in mouse.distinct('name2'):
#         start = []
#         end = []
#         gene_name = gene['name2']
#         for record in mouse.find(name2=gene_name):
#             start.append(int(record['txStart']))
#             end.append(int(record['txEnd']))
#         r = {
#             'strand': record['strand'],
#             'start': min(start),
#             'end': max(end),
#             'chrom': record['chrom']
#         }
#         csv_f.writerow(r)

ranges3 = {}
for n, gene in enumerate(genes):
    if n % 5000 == 0:
        print(n)
    name = gene['name2']
    hits = mouse.find(name2=name)
    if name in ranges3:
        raise Exception("This shouldn't happen")
    ranges3[name] = get_intron(hits, exons)

# cols = ['name2', 'coords', 'strand', 'chrom']
# flat = [dict(zip(cols, [name] + list(map(str, line))))
#         for name, line in ranges3.iteritems()]
# table2.insert_many(flat)

# from matplotlib import pyplot as plt
# import seaborn as sns
# import numpy as np
# dog = [(sum([end-start for start, end in pair]), len(pair))
#        for pair, strand in ranges3.values()]
# dog = np.array(dog)
# dog = dog[dog[:, 0] != 0, :]
# dog2 = [(sum([end-start for start, end in pair]), len(pair))
#        for pair, strand in ranges2.values()]
# dog2 = np.array(dog2)
# dog2 = dog2[dog2[:, 0] != 0, :]
# plt.hist(dog[:, 0], range=(0, 10000), bins=100)
# plt.show()
