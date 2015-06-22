from __future__ import division, print_function

import csv
import gzip
import random
import sys
from collections import Counter
from math import ceil
from itertools import groupby
from intron_locator import IntronGetter, IntronIterator

import dataset
from Bio import SeqIO
from Bio.Seq import Seq
from tempfile import NamedTemporaryFile
from progressbar import ProgressBar
import blaster2
from Bio.SeqRecord import SeqRecord
from pathlib import Path

from oligoarray_designer import Oligoarray, OligoarrayDesigner

csv.field_size_limit(sys.maxsize)  # Prevent field size overflow.


def gc_count(probe):
    return len([1 for c in probe.lower() if c in ['c', 'g']]) / len(probe)


class IntronRetriever(object):
    def get_single_intron(self, gene, chunk_size=1000):
        row = self.intron_getter(gene)
        chrom_code = row['chrom']
        chrom_path = self.chromo_folder.joinpath(
            "{}.fa.masked".format(chrom_code))
        chrom_seq = SeqIO.read(str(chrom_path), 'fasta')
        return self.get_intron(row, chrom_seq, chunk_size)

    def get_intron(self, row, chrom_seq, chunk_size=1000):
        name2 = row['name2']
        if not row['introns']:
            yield None
            return
        for start, end in row['introns']:
            seq = chrom_seq[start:end + 1]
            if row['strand'] == "+":  #TODO: CHECK THAT THIS IS CORRECT
                seq = seq.reverse_complement()
            # Drop masked sequences
            for sub_seq in seq.seq.split("N"):
                if len(sub_seq) >= 20:
                    record = SeqRecord(sub_seq,
                                       name=name2,
                                       id=name2,
                                       description='')
                    for start in range(0, len(record), chunk_size):
                        yield record[start:start + chunk_size]

    def __iter__(self):
        for chromosome in self.intron_getter.table.distinct('chrom'):
            chrom_code = chromosome['chrom']
            chrom_path = self.chromo_folder.joinpath(
                "{}.fa.masked".format(chrom_code))
            chrom_seq = SeqIO.read(str(chrom_path), 'fasta')
            for row in IntronIterator(chrom_code):
                # for record in self.get_intron(row, chrom_seq):
                #     yield record
                gene_records = self.get_intron(row, chrom_seq)
                records = list(gene_records)
                if any(records):
                    yield records

    def __init__(self, organism='mouse'):
        self.organism = organism
        self.chromo_folder = Path("db/chromFaMasked/")
        self.intron_getter = IntronGetter()
        self.tot_introns = self.intron_getter.tot_records


def n_probes(chunk_list, probe_size=35):
    # Get # of probes that could be made from set of gene chunks
    return sum(len(chunk) // probe_size for chunk in chunk_list)


class Designer(object):
    def run_oligoarray(self, chunks, max_probes=150, probe_size=None):
        """
        Chunks is a list of gene fragments lists. [Gene[fragment_1...fragment_m],...Gene_n]
        """
        assert (isinstance(chunks[0], list))
        if not probe_size:
            probe_size = self.probe_size
        flat_probes = []
        with NamedTemporaryFile("w") as fasta_input:
            # Break records into 1000nt chunks
            for gene_chunks in chunks:
                if n_probes(gene_chunks) > max_probes:
                    # Reduce potential probeset size to drop design time
                    while n_probes(gene_chunks, probe_size) > max_probes:
                        gene_chunks = random.sample(gene_chunks,
                                                    len(gene_chunks) - 1)
                SeqIO.write(gene_chunks, fasta_input,
                            'fasta')  # write all sequences to file
            fasta_input.flush()
            res = o.run(fasta_input.name,
                        max_dist=1000,
                        min_length=probe_size,
                        max_length=probe_size,
                        max_tm=100,
                        cross_hyb_temp=72,
                        min_tm=74,
                        secondary_struct_temp=76,
                        max_oligos=100,
                        timeout=15)
            for name, probes in res:
                probes['Percent GC'] = 100 * probes["Probe (5'-> 3')"].map(
                    gc_count)
                flat_probes.append(proves.to_dict().values())
        return flat_probes

    def __iter__(self):
        used_probes = set([row['Name']
                           for row in self.probe_db.distinct("Name")])
        chunksize = 12
        probe_size = self.probe_size
        chunks = []
        intron_retriever = IntronRetriever()
        for n, gene in enumerate(intron_retriever):
            if gene[0].id in used_probes:
                continue
            name = gene[0].id.lower()
            chunks.append(gene)
            if len(chunks) == chunksize:
                probes = self.run_oligoarray(chunks, probe_size=probe_size)
                self.probe_db.insert_many(res)
                chunks = []
            p_bar.update(n)
        p_bar.finish()

    def run(self, probe_size=35):
        pass

    def __init__(self, organism='mouse'):
        intron_db = dataset.connect("sqlite:///db/intron_probes3.db")
        self.probe_db = intron_db[organism]
        self.o = OligoarrayDesigner(
            blast_db='/home/eric/blastdb/old_format/transcribed_mouse2')
        self.probe_size = probe_size




class ProbeFilter(object):
    def run_batch(self):
        hit_val2 = sorted(hit_vals, key=lambda x: x[0].split(',')[0])
        finished_probes = {}
        # Iterate through every probeset
        for target, probe_set in groupby(hit_val2,
                                         key=lambda x: x[0].split(',')[0]):
            pass

    def filter_probes(self, probe_set, target, hits, probe_lookup,
                      max_hits=10000,
                      off_target_thresh=7,
                      probe_num=48, ):
        # Filter out probes with greater than max_hits off target hits
        probe_set = [probe for probe in probe_set if probe[1] < max_hits]
        # First pick probes with no off target hits
        passed_probes = [probe for probe in probe_set if probe[1] == 0]
        remaining_probes = set(probe_set)
        remaining_probes.difference_update(passed_probes)

        # Sort remaining probes by off target hits
        choices = sorted(list(remaining_probes), key=lambda x: x[1])
        choice_index = 0
        off_target = Counter([])

        # Iterate until at least probe_num probes are picked or no more probes are available
        while len(passed_probes) < probe_num and choice_index != len(
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

        # If more than 24 probes chosen optimize on GC
        if len(passed_probes) > probe_num:
            # Get GC counts of every probe
            probe_gc = [(probe[0], gc_count(probe_lookup[probe[0]]))
                        for probe in passed_probes]
            gc_target = 0.55
            gc_range = 0.01
            multiplier = 1
            chosen_gc = []
            # Get closest probe set to 0.55 gc
            while len(chosen_gc) < probe_num:
                gc_min = gc_target - gc_range * multiplier
                gc_max = gc_target + gc_range * multiplier
                chosen_gc = [probe for probe, gc in probe_gc
                             if gc_max >= gc >= gc_min]
                multiplier += 1
            # If too many probes still exists choose a random subset
            if len(chosen_gc) != probe_num:
                chosen_gc = random.sample(chosen_gc, probe_num)
            passed_probes = chosen_gc
        else:
            passed_probes = [probe for probe, n in passed_probes]
        return passed_probes

    def blast2copynum(self, hits, drop_self=True):
        """
        drop_self = True: make sure probes that probes that hit themselves are not counted in blast tally
        """
        hit_vals = []
        for probe_name, matches in hits.iteritems():
            gene_name = probe_name.split(",")[0]
            if any(matches):  # Added incase no matches
                gencode_id, refseq = map(list, zip(*[match.split(',')
                                                     for match in matches]))
            else:
                gencode_id, refseq = [], []
            if drop_self:
                if not gene_name in refseq:
                    continue
                    bad_count += 1
                    # raise Exception("Query not found in hits")
                del gencode_id[refseq.index(gene_name)]
            try:
                hit_vals.append((probe_name, self.get_copynum(gencode_id)))
            except:
                print("Failed @ {}".format(gencode_id))
        hit_vals = sorted(hit_vals, key=lambda x: x[1], reverse=True)
        return hit_vals

    def run_blast(self, probe_lookup, match_thresh, strand=None, db='gencode_tracks_reversed'):
        if not strand:
            strand = self.strand
        fasta_str = "\n".join(">{}\n{}".format(*items)
                              for items in probe_lookup.iteritems())
        res = blaster2.local_blast_query(fasta_str,
                                         db=db)
        hits = blaster2.parse_hits(res,
                                   strand=strand,
                                   match_thresh=match_thresh)
        return hits

    def run(self, flat_probes, gene_name,
            match_thresh=18,
            n_probes=48,
            max_off_target=10000):
        """
        flat_probes is a list of probes.  Sequence only
        """
        probe_lookup = {
            "{},{}".format(gene_name, n): probe
            for n, probe in enumerate(flat_probes)
        }
        res = self.run_blast(probe_lookup, match_thresh, db=self.db)
        hit_vals = self.blast2copynum(res)
        filtered_probes = self.filter_probes(hit_vals, gene_name, res,
                                             probe_lookup,
                                             probe_num=n_probes,
                                             max_hits=max_off_target)
        finished_probes = [probe_lookup[probe] for probe in filtered_probes]
        return finished_probes

    def get_copynum(self, hits):
        off_target = {
            name: self.counts[name]
            for name in hits if name in self.counts.keys()
        }
        false_hits = sum(off_target.values())
        return false_hits

    def __init__(self,
                 db='gencode_tracks_reversed',
                 strand='+',
                 copy_num='embryonic11.5'):
        self.db = db
        if strand == '+':
            self.strand = 1
        else:
            self.strand = -1
        if copy_num == 'embryonic11.5':
            # Get merged embryonic 11.5 encode data
            with open('db/encode_counts.csv', 'r') as f:
                self.counts = {
                    line[0]: float(line[1])
                    for line in csv.reader(f)
                }
        elif copy_num == 'brain':
            with open('db/brain_counts.csv', 'r') as f:
                self.counts = {
                    line[0]: float(line[1])
                    for line in csv.reader(f)
                }
        else:
            raise Exception('need a valid copy_num database')


# intron_db_filtered = dataset.connect("sqlite:///db/intron_probes_filtered.db.bk")
# filtered_probe_table = intron_db_filtered['mouse']
# targets = [row['target'] for row in filtered_probe_table.distinct('target')]
# passed = [t_name for t_name in targets if len(list(filtered_probe_table.find(target=t_name))) >= 24]



def design_introns():
    # First get used probes
    intron_db = dataset.connect("sqlite:///db/intron_probes3.db")
    probe_db = intron_db['mouse']
    used_probes = set([row['Name'] for row in probe_db.distinct("Name")])
    o = OligoarrayDesigner(
        blast_db='/home/eric/blastdb/old_format/transcribed_mouse2')
    chunksize = 12
    probe_size = 35
    chunks = []
    # Check if a good probeset was already designed
    intron_db_filtered = dataset.connect("sqlite:///db/intron_probes_filtered.db.bk")
    filtered_probe_table = intron_db_filtered['mouse']
    intron_retriever = IntronRetriever()
    p_bar = ProgressBar(maxval=intron_retriever.tot_introns).start()
    for n, gene in enumerate(intron_retriever):
        p_bar.update(n)
        gene_chunks = []
        for chunk in gene:
            gene_chunks.append(chunk)
            if n_probes(gene_chunks) > 150: break
        if chunk.id in used_probes: continue
        if len(list(filtered_probe_table.find(target=chunk.id))) >= 24: continue
        # Tagging of introns in db should be added around here
        chunks.append(gene_chunks)
        if len(chunks) == chunksize:
            with NamedTemporaryFile("w") as fasta_input:
                # Break records into 1000nt chunks
                for gene_chunks in chunks:
                    if n_probes(gene_chunks) > 150:
                        # Reduce potential probeset size to drop design time
                        while n_probes(gene_chunks, probe_size) > 150:
                            gene_chunks = random.sample(gene_chunks,
                                                        len(gene_chunks) - 1)
                    SeqIO.write(gene_chunks, fasta_input,
                                'fasta')  # write all sequences to file
                fasta_input.flush()
                res = o.run(fasta_input.name,
                            max_dist=1000,
                            min_length=probe_size,
                            max_length=probe_size,
                            max_tm=100,
                            cross_hyb_temp=72,
                            min_tm=74,
                            secondary_struct_temp=76,
                            max_oligos=100,
                            timeout=15)
            for name, probes in res:
                probes['Percent GC'] = 100 * probes["Probe (5'-> 3')"].map(
                    gc_count)
                probe_db.insert_many(probes.T.to_dict().values())
            p_bar.update(n)
            chunks = []
    p_bar.finish()