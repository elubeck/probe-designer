from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from builtins import *

import copy
import csv
import random
import signal
import sys
from collections import Counter, defaultdict
from itertools import groupby
import json
from tempfile import NamedTemporaryFile

import dataset
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from operator import itemgetter
from pathlib import Path
from progressbar import ProgressBar

import blaster2
import get_seq
from oligoarray_designer import OligoarrayDesigner

csv.field_size_limit(sys.maxsize)  # Prevent field size overflow.


def gc_count(probe):
    """
    Returns the % GC for a string
    """
    return len([1 for c in probe.lower() if c in ['c', 'g']]) / len(probe)


class TimeoutError(Exception):
    pass


class timeout(object):
    def __init__(self, seconds=1, error_message='Timeout'):
        self.seconds = seconds
        self.error_message = error_message

    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)

    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)

    def __exit__(self, type, value, traceback):
        signal.alarm(0)


# TODO: FINISH ME
class mRNARetriever(object):
    """
    Gets mRNAs from mouse genome.
    """

    def get_mRNA(self, target, chunk_size=1000):
        """
        Given a target(refseq name) returns the sequence of probes in increments of chunk_size.  Aligns all variants of a gene to ensure correct common sequence amongst alternatively spliced transcripts.
        """
        hits = []
        for file_name in Path("db/").glob("mouse.*.rna.fna"):
            with open(str(file_name), 'r') as fin:
                for record in SeqIO.parse(fin, 'fasta'):
                    if "({})".format(
                        target.lower()) in record.description.lower():
                        hits.append(record)
        for sub_seq in get_seq.CDS(target).align_seqs(hits):
            if len(sub_seq) >= 20:
                record = SeqRecord(Seq(sub_seq),
                                   name=target,
                                   id=target,
                                   description='')
                for start in range(0, len(record), chunk_size):
                    yield record[start:start + chunk_size]

    def __iter__(self):
        """
        Iterates over entire refseq gene list.
        """
        with open("db/Refseq2Gene.tsv", "r") as f_in:
            ref_names = {
                gene_name
                for rnum, gene_name in csv.reader(f_in,
                                                  delimiter='\t')
            }
            for gene in ref_names:
                for record in self.get_mRNA(gene):
                    yield record

    def __init__(self, organism='mouse'):
        self.organism = organism


def n_probes(chunk_list, probe_size=35):
    # Get # of probes that could be made from set of gene chunks
    return sum(len(chunk) // probe_size for chunk in chunk_list)


class ProbeDesigner(object):
    def run_oligoarray(self, chunks, n_target_probes=150):
        """
        Wrapper code to run oligoarray over a list of sequences.
        """
        with NamedTemporaryFile("w") as fasta_input:
            # Break records into 1000nt chunks
            for gene_chunks in chunks:
                if n_probes(gene_chunks, self.probe_size) > n_target_probes:
                    # Reduce potential probeset size to drop design time
                    while n_probes(gene_chunks,
                                   self.probe_size) > n_target_probes:
                        gene_chunks = random.sample(gene_chunks,
                                                    len(gene_chunks) - 1)
                SeqIO.write(gene_chunks, fasta_input, 'fasta')
            fasta_input.flush()
            res = self.o.run(fasta_input.name,
                             max_dist=self.max_dist,
                             min_length=self.probe_size,
                             max_length=self.probe_size,
                             max_tm=self.max_tm,
                             cross_hyb_temp=self.cross_hyb_temp,
                             min_tm=self.min_tm,
                             secondary_struct_temp=self.secondary_struct_temp,
                             max_oligos=n_target_probes,
                             timeout=15)
        for name, probes in res:
            probes['Percent GC'
                   ] = 100 * probes["Probe (5'-> 3')"].map(gc_count)
            yield probes

    def design(self, target):
        gene = [list(mRNARetriever().get_mRNA(target))]
        for probe in self.run_oligoarray(gene):
            print("TODO: ADD to DB")

    def batch_design(self, target_list, chunksize=3):
        p_bar = ProgressBar(maxval=len(target_list)).start()
        chunks = []
        n = 0
        probe_list = []
        while True:
            if len(chunks) == chunksize or n == len(target_list):
                for name, probes in self.run_oligoarray(chunks):
                    for probe in probes.T.to_dict().values():
                        probe_list.append(probe)
                    # TODO: ADD to DB
                chunks = []
                p_bar.update(n)
                if n == len(target_list):
                    break
            gene_name = target_list[n]
            chunks.append(list(mRNARetriever().get_mRNA(gene_name)))
            n += 1
        p_bar.finish()
        return probe_list

    def __init__(self,
                 db='/home/eric/blastdb/old_format/transcribed_mouse',
                 max_dist=1000,
                 probe_size=35,
                 max_tm=100,
                 cross_hyb_temp=72,
                 min_tm=74,
                 secondary_struct_temp=76, ):
        self.o = OligoarrayDesigner(blast_db=db)
        self.probe_size = probe_size
        self.max_dist = max_dist
        self.cross_hyb_temp = cross_hyb_temp
        self.min_tm = min_tm
        self.max_tm = max_tm
        self.secondary_struct_temp = secondary_struct_temp


def filter_probes_by_blast(probes):
    def probe_to_fasta(probe):
        return ">{Name},{Probe #}\n{Probe (5'-> 3')}\n".format(**probe)

    iv = [probe
          for probe_set in probes for probe in probe_set.T.to_dict().values()]
    probe_lookup = {
        "{Name},{Probe #}".format(**probe): probe["Probe (5'-> 3')"]
        for probe in probes
    }
    fasta_str = "".join(probe_to_fasta(probe) for probe in probes)
    # for n, probe in enumerate(probes.values):
    #     probe_lookup.append(("{},{}".format(gene, n), probe))
    # probe_lookup = dict(probe_lookup )
    # fasta_str = "\n".join(">{}\n{}".format(k, v) for k,v in probe_lookup.iteritems())

    strand = 1
    res = blaster2.local_blast_query(fasta_str, db='gencode_tracks_reversed')
    hits = blaster2.parse_hits(res, strand=strand, match_thresh=18)

    import pdb
    pdb.set_trace()
    hit_vals = []
    for probe_name, matches in hits.iteritems():
        gene_name = probe_name.split(",")[0]
        gencode_id, refseq = map(list, zip(*[match.split(',')
                                             for match in matches]))
        if not gene_name in refseq:
            continue
            bad_count += 1
            # raise Exception("Query not found in hits")
        del gencode_id[refseq.index(gene_name)]
        try:
            hit_vals.append((probe_name, blaster2.get_copynum(gencode_id)))
        except:
            print("Failed @ {}".format(gencode_id))
    hit_vals = sorted(hit_vals, key=lambda x: x[1], reverse=True)

    def gc_count(probe):
        return len([1 for c in probe.lower() if c in ['c', 'g']]) / len(probe)

    # Iterate everyg ene
    from itertools import groupby
    from collections import Counter
    hit_val2 = sorted(hit_vals, key=lambda x: x[0].split(',')[0])
    off_target_thresh = 7
    finished_probes = {}
    # Iterate through every probeset
    for target, probe_set in groupby(hit_val2,
                                     key=lambda x: x[0].split(',')[0]):
        probe_set = list(probe_set)

        # First pick probes with no off target hits
        passed_probes = [probe for probe in probe_set if probe[1] == 0]
        remaining_probes = set(probe_set)
        remaining_probes.difference_update(passed_probes)

        # Sort remaining probes by off target hits
        choices = sorted(list(remaining_probes), key=lambda x: x[1])
        choice_index = 0
        off_target = Counter([])

        # Iterate until at least 24 probes are picked or no more probes are available
        while len(passed_probes) < 24 and choice_index != len(remaining_probes):
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
        if len(passed_probes) > 24:
            # Get GC counts of every probe
            probe_gc = [(probe[0], gc_count(probe_lookup[probe[0]]))
                        for probe in passed_probes]
            gc_target = 0.55
            gc_range = 0.01
            multiplier = 1
            chosen_gc = []
            # Get closest probe set to 0.55 gc
            while len(chosen_gc) < 24:
                gc_min = gc_target - gc_range * multiplier
                gc_max = gc_target + gc_range * multiplier
                chosen_gc = [probe for probe, gc in probe_gc
                             if gc_max >= gc >= gc_min]
                multiplier += 1

            # If too many probes still exists choose a random subset
            if len(chosen_gc) != 24:
                chosen_gc = random.sample(chosen_gc, 24)
            passed_probes = chosen_gc
        else:
            passed_probes = [probe for probe, n in passed_probes]
        finished_probes[target] = [probe_lookup[probe]
                                   for probe in passed_probes]
        return finished_probes


def sub_seq_splitter(seq, size,
                     gc_target=0.55,
                     stride=1,
                     spacing=2,
                     gc_min=None,
                     gc_max=None,
                     o_mode=1):
    """
    Takes sequence and designs best probes of given size closests to gc_target.
    """
    assert (isinstance(seq, unicode))
    # Generate all probes
    # Probes start and ending position includes spacing.
    probes = []
    for seed in range(0, len(seq) - size, stride):
        p_seq = seq[seed:seed + size]
        p_gc = gc_count(p_seq)
        if gc_min and p_gc < gc_min:
            continue
        if gc_max and p_gc > gc_max:
            continue
        probes.append({
            'start': seed - spacing,
            'end': seed + size + spacing,
            'gc': p_gc,
            'seq': p_seq
        })
    probes = sorted(probes, key=lambda x: abs(gc_target - x['gc']))

    # Build list of overlapping probes
    probe_indices = [set(range(probe['start'], probe['end']))
                     for probe in probes]
    if o_mode == 1:
        overlapping = []
        for probe_n, probe_ind in enumerate(probe_indices):
            for i in range(len(probe_indices)):
                if i != probe_n:
                    if any(probe_ind & probe_indices[i]):
                        overlapping.append((probes[probe_n], probes[i]))
        overlapping_2 = {
            probe['seq']: [hit for query, hit in over_probes]
            for probe, over_probes in groupby(overlapping,
                                              key=lambda x: x[0])
        }
    if o_mode == 2:
        seq_lookup = {probe['seq']: probe for probe in probes}
        overlapping = defaultdict(list)
        for probe_n, probe_ind in enumerate(probe_indices):
            for i in range(probe_n, len(probe_indices)):
                if i != probe_n:
                    if any(probe_ind & probe_indices[i]):
                        overlapping[probes[probe_n]['seq']].append(probes[i])
        overlapping_b = copy.deepcopy(overlapping)
        for k, v in overlapping.iteritems():
            pl = seq_lookup[k]
            for i2 in v:
                overlapping_b[i2['seq']].append(pl)
        # Get probes that don't overlap anything!
        for seq in set(seq_lookup.keys()).difference(overlapping_b.keys()):
            overlapping_b[seq] = []
        overlapping_2 = {
            seq: [hit for query, hit in over_probes][0]
            for seq, over_probes in groupby(overlapping_b.items(),
                                            key=lambda x: x[0])
        }

        # Group probes by difference from target gc
    probe_groups = [(k, list(v))
                    for k, v in groupby(probes,
                                        key=lambda x: abs(gc_target - x['gc']))
                    ]
    probe_index = 0
    passed_probes = []
    disallowed_probes = []

    # Iterate by probe GC and choose non-overlapping probes
    while probe_index != len(probe_groups):
        gc, probe_group = probe_groups[probe_index]
        # Filter already overlapping probes and convert list of dict to list of seqs
        # TODO:  This(set vs list comprehension) is a dirty hack that prevents probes that hit two different places from being chosen twice
        probe_group = {
            probe['seq']
            for probe in probe_group if probe['seq'] not in disallowed_probes
        }
        # Iterate until no probes overlap in probe_group
        if len(probe_group) > 1:
            n_iters = 0
            while True:
                om = [(probe_seq, [overlapped
                                   for overlapped in overlapping_2[probe_seq]
                                   if overlapped['seq'] in probe_group])
                      for probe_seq in probe_group]
                om = sorted(om, key=lambda x: len(x[1]), reverse=True)
                # If worst probe doesn't overlap any other probes
                if len(om[0][1]) == 0:
                    break
                # Otherwise throw out worst probe
                else:
                    probe_group = [probe_seq for probe_seq, overind in om[1:]]
                n_iters += 1
            probe_group = [probe_seq for probe_seq, overind in om[1:]]
        for probe_seq in probe_group:
            passed_probes.append(probe_seq)
            for overlapped in overlapping_2[probe_seq]:
                disallowed_probes.append(overlapped['seq'])
        probe_index += 1

    # DEBUG_START
    passed_or = defaultdict(list)
    for probe in passed_probes:
        for probe_m in probes:
            if probe_m['seq'] == probe:
                passed_or[probe].append(probe_m)
    dog_head = sorted([i['start'] for p in passed_or.values() for i in p])
    dog_poop = [i2 - i1 for i2, i1 in zip(dog_head[1:], dog_head[:-1])]
    for v in dog_poop:
        if v < 37:
            raise Exception("Probe design fail")
    # DEBUG_END
    return passed_probes


def probe_set_refiner(pset_i, block_size=18):
    """
    Checks that probes don't hit the same target of block_size.  Drops probes if they do.  Returns a filtered probeset.
    Makes the assumption that all probes from different sets that hit the same target must be dropped.  This assumption can be faulty in the case of non-mRNA targeting sequences such as adapters added to probes.  Prefered behavior in this case would be to drop nothing or only mRNA binding sequences.
    """
    pset = copy.deepcopy(pset_i)
    p_lookup = defaultdict(list)
    flat_pfrags = []
    for gene, probes in pset.iteritems():
        for probe in probes:
            for i in range(0, len(probe) - block_size):
                p_lookup[probe[i:i + block_size]].append(gene)
                flat_pfrags.append(probe[i:i + block_size])
    # Get probes that hit multiple targets
    counts = [(k, v) for k, v in Counter(flat_pfrags).iteritems() if v > 1]
    counts = [(bad_probe, n_hits) for bad_probe, n_hits in counts
              if len(set(p_lookup[bad_probe])) > 1]
    for probe, n_hits in counts:
        for target in set(p_lookup[probe]):
            for n, probe_seq in enumerate(pset[target]):
                if probe in probe_seq:
                    del pset[target][n]
    return pset


def batch_design(genes, max_time=180):
    db = dataset.connect("sqlite:///db/mrna_probes.db")
    p_table = db['unfiltered_probes']
    cds_table = db['cds_table']
    mr = mRNARetriever()
    cds_records = {}
    gene_probes = defaultdict(list)
    pbar = ProgressBar(maxval=len(genes)).start()
    failed = []
    used = {g['target'] for g in p_table.distinct('target')}
    for n, gene in enumerate(genes):
        pbar.update(n)
        if gene in used:
            print("Found {}".format(gene))
            gene_probes[gene] = [probe['seq']
                                 for probe in p_table.find(target=gene)]
            continue
        try:
            with timeout(seconds=max_time):
                try:
                    gene_records = [
                        str(record.seq.reverse_complement())
                        for record in mr.get_mRNA(gene,
                                                  chunk_size=100000)
                    ]
                    rec_table = [
                        {"target": gene,
                         "seq": cds,
                         "number": rec_num}
                        for rec_num, cds in enumerate(gene_records)
                    ]
                    cds_table.insert_many(rec_table)
                except:
                    failed.append(gene)
                cds_records[gene] = gene_records
                for chunk in cds_records[gene]:
                    for probe in sub_seq_splitter(str(chunk), 35,
                                                  gc_min=0.35,
                                                  o_mode=2):
                        gene_probes[gene].append(probe)
                    flat_p = [{'target': gene,
                               'seq': probe} for probe in gene_probes[gene]]
                    p_table.insert_many(flat_p)
        except TimeoutError:
            print("Timedout at {}".format(gene))
    pbar.finish()
    return gene_probes


def design_step(gene, max_time=180):
    rr = RNARetriever2()
    probes = []
    try:
        with timeout(seconds=max_time):
            try:
                gene_records = [
                    str(record.seq)
                    for record in rr.get_single_rna(gene,
                                                    chunk_size=100000)
                ]
            except:
                print("Record lookup fail at {}".format(gene))
                return ("FAILED", [], [])
            for chunk in gene_records:
                for probe in sub_seq_splitter(str(chunk), 35,
                                              gc_min=0.35,
                                              o_mode=2):
                    probes.append(probe)
            return gene, probes, gene_records
    except TimeoutError:
        print("Timedout at {}".format(gene))
        return ("FAILED", [], [])


def batch_design2(genes, max_time=180):
    db = dataset.connect("sqlite:///db/mrna_probes2.db")
    p_table = db['unfiltered_probes']
    cds_table = db['cds_table']
    mr = mRNARetriever()
    gene_probes = defaultdict(list)
    pbar = ProgressBar(maxval=len(genes)).start()
    used = {g['target'] for g in p_table.distinct('target')}
    rr = RNARetriever2()
    gi = (gene for gene in genes if gene not in used)
    gene_probes = {}
    from multiprocessing import Pool, cpu_count
    from functools import partial
    d_wrap = partial(design_step, max_time=max_time)
    p = Pool(cpu_count())
    for n, vals in enumerate(p.imap(d_wrap, gi)):
        gene, probes, gene_records = vals
        if gene == "FAILED": continue
        rec_table = [
            {"target": gene,
             "seq": cds,
             "number": rec_num} for rec_num, cds in enumerate(gene_records)
        ]
        cds_table.insert_many(rec_table)
        flat_p = [{'target': gene, 'seq': probe} for probe in probes]
        p_table.insert_many(flat_p)
        gene_probes[gene] = probes
        pbar.update(n)
    pbar.finish()
    return gene_probes


class SequenceGetter(object):
    """
    Gets sequences from refgenome table
    """

    def get_range(self, sequence, min_size=20):
        #Find longest contiguous sequences
        ranges = []
        for k, g in groupby(enumerate(sequence), lambda (i, x): i - x):
            group = [itemgetter(1)(v) for v in g]
            start = group[0]
            end = group[-1]
            if end - start > min_size:
                ranges.append((start, end))
        return ranges

    def get_positions(self, gene):
        """
        Finds position of exons for gene
        """

        def splitter(string):
            return [int(k) for k in string.strip(',').split(',')]

        if gene['strand'] == '+':
            start = int(gene['txStart'])
            end = int(gene['txEnd'])
            e_start = splitter(gene['exonStarts'])
            e_end = splitter(gene['exonEnds'])
        elif gene['strand'] == '-':
            # Flip things around
            start = int(gene['txEnd'])
            end = int(gene['txStart'])
            e_start = splitter(gene['exonStarts'])[::-1]
            e_end = splitter(gene['exonEnds'])[::-1]
        return (start, end, e_start, e_end)

    def __init__(self, table=None, ):
        if table is None:
            db = dataset.connect("sqlite:///db/refGene.db")
            self.table = db['mouse']
        else:
            self.table = table
        self.tot_records = len(list(self.table.distinct('name2')))


class RNAGetter(SequenceGetter):
    def get_rna(self, gene_isoforms, ):
        """
        Get all rnas from a gene and return in transcribed order.
        Only chooses exons in every isoform.
        """
        rna_set = []
        for gene in gene_isoforms:
            tx_start, tx_end, e_start, e_end = self.get_positions(gene)
            exon_pos = {
                p
                for end, start in zip(e_end, e_start)
                for p in range(*sorted([end, start]))
            }
            rna_set.append(exon_pos)

        # Get Regions in Every Isoform
        conserved_regions = set.intersection(*rna_set)
        rna_regions = []
        for contiguous_frag in self.get_range(conserved_regions):
            rna_regions.append(contiguous_frag)
        # Remove duplicate rnas from spliceoforms and return sorted result
        rnas = sorted(set(rna_regions), key=lambda x: x[0])
        if gene['strand'] == "+":
            return rnas
        elif gene['strand'] == '-':
            return rnas[::-1]

    def run_gene(self, gene):
        """
        Wrapper for finding results for a gene
        """
        hits = list(self.table.find(name2=gene))
        return {
            'chrom': hits[0]['chrom'],
            'strand': hits[0]['strand'],
            'name2': gene,
            'exons': self.get_rna(hits, )
        }

    def __init__(self, table=None, ):
        super(RNAGetter, self).__init__(table, )


class RNAIterator(RNAGetter):
    def __iter__(self):
        for chromosome in self.table.distinct('chrom'):
            chrom_code = chromosome['chrom']
            for row in self.table.distinct('name2', chrom=chrom_code):
                name2 = row['name2']
                yield self.run_gene(name2)

    def __init__(self, ):
        super(self.__class__, self).__init__()


def gc_count(probe):
    return len([1 for c in probe.lower() if c in ['c', 'g']]) / len(probe)


class RNARetriever2(object):
    def get_single_rna(self, gene, chunk_size=1000):
        row = self.rna_getter.run_gene(gene)
        chrom_code = row['chrom']
        chrom_path = self.chromo_folder.joinpath(
            "{}.fa.masked".format(chrom_code))
        chrom_seq = SeqIO.read(str(chrom_path), 'fasta')
        return self.get_rna(row, chrom_seq, chunk_size)

    def get_rna(self, row, chrom_seq, chunk_size=1000):
        name2 = row['name2']
        if not row['exons']:
            yield None
            return
        for start, end in row['exons']:
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
        for chromosome in self.rna_getter.table.distinct('chrom'):
            chrom_code = chromosome['chrom']
            chrom_path = self.chromo_folder.joinpath(
                "{}.fa.masked".format(chrom_code))
            chrom_seq = SeqIO.read(str(chrom_path), 'fasta')
            for row in ChromIntronIterator(chrom_code,
                                           skip_rnas=self.skip_probes):
                gene_records = self.get_rna(row, chrom_seq)
                records = list(gene_records)
                if any(records):
                    yield records

    def __init__(self, organism='mouse', skip_probes=None):
        self.organism = organism
        self.chromo_folder = Path("db/{}/chromFaMasked/".format(organism))
        self.rna_getter = RNAGetter()
        if not skip_probes:
            skip_probes = []
        self.skip_probes = skip_probes
        self.tot_rnas = self.rna_getter.tot_records - len(skip_probes)


if __name__ == '__main__':
    pass
