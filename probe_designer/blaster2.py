from __future__ import print_function, division
import csv
import collections
import signal
import sys
from collections import Counter, defaultdict
from time import sleep

import dataset
import pandas as pd
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline as blastn
from future.builtins import object
from tempfile import NamedTemporaryFile

# Get merged embryonic 11.5 encode data
with open('db/encode_counts.csv', 'r') as f:
    counts = {line[0]: float(line[1]) for line in csv.reader(f)}


def local_blast_query(query, db='gencode_tracks.fas'):
    """
    Blast fasta query via NCBI blastn.
    :param query: Fasta string
    :return: a handle to an xml formatted blast result
    """
    output_file = "temp/blast_res.xml"
    with NamedTemporaryFile("w") as f:
        f.write(query)
        f.flush()
        res = blastn(query=f.name, db=db, out=output_file, outfmt=5,
                    word_size=11, num_threads=8)
        res()
    output_handle = open(output_file, 'r')
    return output_handle



def parse_hits(handle, strand=-1, match_thresh=13):
    """
    Parses the results from a blast search.  Returns all blast hits on a given strand
    :param handle: results xml handle from blast search
    :param strand: which strand to search for hits.  By default searches complementary strand (-1).
    :return: dictionary of query:[target...] hits for each query
    """
    #Parse blast hits
    gene_hits = defaultdict(list)
    for record in NCBIXML.parse(handle):
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                #Check that hit is on opposite strand
                if hsp.frame[1] == strand:
                    block_len = [len(block) for block in hsp.match.split(" ")]
                    if any(True for b in block_len if b >= match_thresh):
                        gene_hits[record.query].append(alignment.hit_def)
    return gene_hits

def blast_probes(probes, timeout=120, debug=False,):
    """
    Runs entire blast routine for a given gene returning
    only the complementary hits for a probeset.  Stores result in a db.
    :param debug: Whether to print program output
    :param local: Run local or global blast
    :return:
    """
    res = local_blast_query(fasta_query,)
    gene_hits = parse_hits(res)
    res.close()
    return gene_hits


def flatten(l):
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el

def probe_to_fasta(probe):
    return ">{Name},{Probe #}\n{Probe (5'-> 3')}\n".format(**probe)
                
def get_copynum(hits):
    off_target = {name: counts[name] for name in hits}
    false_hits = sum(off_target.values())
    return false_hits

intron_db = dataset.connect("sqlite:///db/intron_probes.db")
blast_db = dataset.connect("sqlite:///db/intron_blast.db")
blast_table = blast_db['mouse']
probe_db = intron_db['mouse']
used_genes = list(probe_db.distinct("Name"))
probes = list(probe_db.find(Name='Atg5'))
probe_lookup = {"{Name},{Probe #}".format(**probe):probe["Probe (5'-> 3')"] for probe in probes}
fasta_str = "".join(probe_to_fasta(probe) for probe in probes)
res_handle = local_blast_query(fasta_str)
match_thresh = 13
strand = "+"
if strand == "+":
    skey = 1
elif strand == '-':
    skey = -1
else:
    raise Exception("BOB")
hits = parse_hits(res_handle, strand=skey, match_thresh=match_thresh)
db_vals = [{'query':query,'sequence':probe_lookup[query], 'match_len':match_thresh,
            'target': query.split(',')[0], 'strand': strand,
            'hits_encode':",".join(x.split(",")[0] for x in hit_l),
            'hits_refseq':",".join(x.split(",")[1] for x in hit_l),
            }
            for query, hit_l in hits.iteritems()]
blast_table.insert_many(db_vals)

# hit_vals = []
# for probe_name, matches in hits.iteritems():
#     gene_name, probe_num = name.split(",")
#     gencode_id, refseq = map(list, zip(*[match.split(',') for match in matches]))
#     if not gene_name in refseq:
#         raise Exception("Query not found in hits")
#     del gencode_id[refseq.index(gene_name)]
#     hit_vals.append((probe_name, get_copynum(gencode_id)))
# hit_vals = sorted(hit_vals, key=lambda x: x[1], reverse=True)

# from itertools import groupby

# off_thresh = 8
# min_probes = 24
# passed_probes = hits.keys() 
# while len(passed_probes) >= min_probes:
#     flat = [(hit, k) for k, hit_l in hits.iteritems() if k in passed_probes
#             for hit in hit_l]
#     grouped_hits = [(k, [probe for k1, probe in probes])
#                     for k, probes in groupby(flat, key=lambda x: x[0])]
#     grouped_hits = sorted(grouped_hits, key=lambda x: len(x[1]), reverse=True)
#     bad_probes = grouped_hits[0][1]
#     probe_rank = [(probe, len(True for k, group in grouped_hits if probe in group))
#         for probe in bad_probes]
#     worst_probe = max(probe_rank, key=lambda x: x[1])
