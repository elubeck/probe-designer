from __future__ import print_function, division
import csv
import collections
import signal
import sys
from collections import Counter, defaultdict
from time import sleep
from itertools import groupby
from collections import Counter
import random

import dataset
import pandas as pd
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline as blastn
from future.builtins import object
from tempfile import NamedTemporaryFile

# Get merged embryonic 11.5 encode data
with open('db/encode_counts.csv', 'r') as f:
    counts = {line[0]: float(line[1]) for line in csv.reader(f)}


def local_blast_query(query, db='gencode_tracks_reversed'):
    """
    Blast fasta query via NCBI blastn.
    :param query: Fasta string
    :return: a handle to an xml formatted blast result
    """
    output_handle = NamedTemporaryFile(suffix='.xml')
    with NamedTemporaryFile("w") as f:
        f.write(query)
        f.flush()
        res = blastn(query=f.name,
                     db=db,
                     out=output_handle.name,
                     outfmt=5,
                     word_size=11,
                     num_threads=8)
        res()
    return output_handle


def parse_hits(handle, strand=1, match_thresh=13):
    """
    Parses the results from a blast search.  Returns all blast hits on a given strand
    :param handle: results xml handle from blast search
    :param strand: which strand to search for hits.  By default searches complementary strand (-1).
    :return: dictionary of query:[target...] hits for each query
    """
    #Parse blast hits
    gene_hits = defaultdict(list)
    for record in NCBIXML.parse(handle):
        gene_hits[record.query] = [
        ]  # This was added so probes that couldn't match anything weren't dropped
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                #Check that hit is on opposite strand
                if hsp.frame[1] == strand:
                    block_len = [len(block) for block in hsp.match.split(" ")]
                    if any(True for b in block_len if b >= match_thresh):
                        gene_hits[record.query].append(alignment.hit_def)
    handle.close()
    return gene_hits


def blast_probes(probes, timeout=120, debug=False, ):
    """
    Runs entire blast routine for a given gene returning
    only the complementary hits for a probeset.  Stores result in a db.
    :param debug: Whether to print program output
    :param local: Run local or global blast
    :return:
    """
    res = local_blast_query(fasta_query, )
    gene_hits = parse_hits(res)
    res.close()
    return gene_hits


def probe_to_fasta(probe):
    return ">{Name},{Probe #}\n{Probe (5'-> 3')}\n".format(**probe)


def get_copynum(hits):
    off_target = {name: counts[name] for name in hits if name in counts.keys()}
    false_hits = sum(off_target.values())
    return false_hits


def gc_count(probe):
    return len([1 for c in probe.lower() if c in ['c', 'g']]) / len(probe)


def filter_probes(probe_table, off_target_thresh=7, min_probes=24):
    # First get copy #'s of everything
    c_num = []
    probe_lookup = {}
    for probe in probe_table:
        target = probe['target']
        target_pos = probe['hits_refseq'].split(',').index(target)
        encode_refs = [
            encode_name
            for n, encode_name in enumerate(probe['hits_encode'].split(','))
            if n != target_pos
        ]
        probe['hit_num'] = get_copynum(encode_refs)
        probe_lookup[probe['query']] = probe['sequence']
        c_num.append((probe['query'], probe['hit_num']))

    # First pick probes with no off target hits
    passed_probes = [probe for probe in c_num if probe[1] == 0]
    remaining_probes = set(c_num)
    remaining_probes.difference_update(passed_probes)

    # Sort remaining probes by off target hits
    choices = sorted(list(remaining_probes), key=lambda x: x[1])
    choice_index = 0
    off_target = Counter([])

    # Iterate until at least 24 probes are picked or no more probes are available
    while len(passed_probes) < min_probes and choice_index != len(
        remaining_probes):
        selected = choices[choice_index]
        # Get off target hits for chosen probe
        off_target_hit = [k for k in hits[selected[0]]
                          if target != k.split(',')[1]]
        test_counter = off_target + Counter(off_target_hit)

        # Check if adding this probe makes anything go over off target threshold
        over_thresh = [True for k in test_counter.values()
                       if k >= off_target_thresh]
        if not any(over_thresh):
            passed_probes.append(selected)
            off_target = test_counter
        choice_index += 1

    # If more than min_probes probes chosen optimize on GC
    if len(passed_probes) > min_probes:
        # Get GC counts of every probe
        probe_gc = [(probe[0], gc_count(probe_lookup[probe[0]]))
                    for probe in passed_probes]
        gc_target = 0.55
        gc_range = 0.01
        multiplier = 1
        chosen_gc = []
        # Get closest probe set to 0.55 gc
        while len(chosen_gc) < min_probes:
            gc_min = gc_target - gc_range * multiplier
            gc_max = gc_target + gc_range * multiplier
            chosen_gc = [probe for probe, gc in probe_gc
                         if gc_max >= gc >= gc_min]
            multiplier += 1

        # If too many probes still exists choose a random subset
        if len(chosen_gc) >= min_probes:
            chosen_gc = random.sample(chosen_gc, min_probes)
        passed_probes = chosen_gc
    else:
        passed_probes = [probe for probe, n in passed_probes]
    finished_probes = [probe_lookup[probe] for probe in passed_probes]
    return finished_probes

# intron_db = dataset.connect("sqlite:///db/intron_probes.db")
# blast_db = dataset.connect("sqlite:///db/intron_blast.db")
# blast_table = blast_db['mouse']
# probe_db = intron_db['mouse']
# dog_cat = filter_probes(list(blast_table.find(target='Pgk1')), )
# strand = "+"
# match_thresh = 18
# used_probes = set([row['target'] for row in blast_table.distinct("target", match_len=match_thresh)])
# used_genes = list(probe_db.distinct("Name"))
# from progressbar import ProgressBar
# p_bar = ProgressBar(maxval=len(used_genes)).start()
# for ngene, gene in enumerate(used_genes):
#     gname = gene['Name']
#     if gname in used_probes:
#         continue
#     # Format probes for blast and later lookup
#     probes = list(probe_db.find(Name=gname))
#     probe_lookup = {"{Name},{Probe #}".format(**probe):probe["Probe (5'-> 3')"] for probe in probes}
#     fasta_str = "".join(probe_to_fasta(probe) for probe in probes)
#     res_handle = local_blast_query(fasta_str)
#     if strand == "+":
#         skey = 1
#     elif strand == '-':
#         skey = -1
#     else:
#         raise Exception("BOB")
#     hits = parse_hits(res_handle, strand=skey, match_thresh=match_thresh)
#     db_vals = [{'query':query,'sequence':probe_lookup[query], 'match_len':match_thresh,
#                 'target': query.split(',')[0], 'strand': strand,
#                 'hits_encode':",".join(x.split(",")[0] for x in hit_l),
#                 'hits_refseq':",".join(x.split(",")[1] for x in hit_l),}
#                 for query, hit_l in hits.iteritems()]
#     blast_table.insert_many(db_vals)
#     p_bar.update(ngene)
# p_bar.finish()

# # hit_vals = []
# # for probe_name, matches in hits.iteritems():
# #     gene_name, probe_num = name.split(",")
# #     gencode_id, refseq = map(list, zip(*[match.split(',') for match in matches]))
# #     if not gene_name in refseq:
# #         raise Exception("Query not found in hits")
# #     del gencode_id[refseq.index(gene_name)]
# #     hit_vals.append((probe_name, get_copynum(gencode_id)))
# # hit_vals = sorted(hit_vals, key=lambda x: x[1], reverse=True)

# # from itertools import groupby

# # off_thresh = 8
# # min_probes = 24
# # passed_probes = hits.keys() 
# # while len(passed_probes) >= min_probes:
# #     flat = [(hit, k) for k, hit_l in hits.iteritems() if k in passed_probes
# #             for hit in hit_l]
# #     grouped_hits = [(k, [probe for k1, probe in probes])
# #                     for k, probes in groupby(flat, key=lambda x: x[0])]
# #     grouped_hits = sorted(grouped_hits, key=lambda x: len(x[1]), reverse=True)
# #     bad_probes = grouped_hits[0][1]
# #     probe_rank = [(probe, len(True for k, group in grouped_hits if probe in group))
# #         for probe in bad_probes]
# #     worst_probe = max(probe_rank, key=lambda x: x[1])
