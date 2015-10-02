from __future__ import division, print_function

import csv
import random
import sys
from tempfile import NamedTemporaryFile
from pathlib import Path

import dataset
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from intron_locator import IntronGetter, ChromIntronIterator
from progressbar import ProgressBar
from oligoarray_designer import OligoarrayDesigner
from utils.misc import  gc_count, n_probes

csv.field_size_limit(sys.maxsize)  # Prevent field size overflow.



class IntronRetriever(object):
    def get_single_intron(self, gene, chunk_size=1000):
        row = self.intron_getter.run_gene(gene)
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
            seq = chrom_seq[start:end]
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
        chrom_i = self.intron_getter.table.distinct('chrom')
        if self.reversed is True:
            chrom_i = reversed(list(chrom_i))
        for chromosome in chrom_i:
            chrom_code = chromosome['chrom']
            chrom_path = self.chromo_folder.joinpath(
                "{}.fa.masked".format(chrom_code))
            chrom_seq = SeqIO.read(str(chrom_path), 'fasta')
            for row in ChromIntronIterator(chrom_code,
                                           skip_introns=self.skip_probes):
                gene_records = self.get_intron(row, chrom_seq)
                records = list(gene_records)
                if any(records):
                    yield records

    def __init__(self, organism='mouse', skip_probes=None, reversed=False):
        self.organism = organism
        self.chromo_folder = Path("db/{}/chromFaMasked/".format(organism))
        self.intron_getter = IntronGetter(in_all=True)
        if skip_probes is None:
            self.skip_probes = []
        else:
            self.skip_probes = skip_probes
        self.tot_introns = self.intron_getter.tot_records - len(skip_probes)
        self.reversed = reversed



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


# intron_db_filtered = dataset.connect("sqlite:///db/intron_probes_filtered.db.bk")
# filtered_probe_table = intron_db_filtered['mouse']
# targets = [row['target'] for row in filtered_probe_table.distinct('target')]
# passed = [t_name for t_name in targets if len(list(filtered_probe_table.find(target=t_name))) >= 24]


def design_introns(reversed=False):
    # First get used probes
    intron_db = dataset.connect("sqlite:///db/intron_probes_10k_3.db")
    probe_db = intron_db['mouse']
    used_probes = set([row['Name'] for row in probe_db.distinct("Name")])
    o = OligoarrayDesigner(
        blast_db='/home/eric/blastdb/old_format/transcribed_mouse2')
    chunksize = 12
    probe_size = 35
    chunks = []
    # Check if a good probeset was already designed
    # intron_db_filtered = dataset.connect(
    #     "sqlite:///db/intron_probes_filtered.db.bk")
    # filtered_probe_table = intron_db_filtered['mouse']
    intron_retriever = IntronRetriever(skip_probes=used_probes,
                                       reversed=reversed)
    p_bar = ProgressBar(maxval=intron_retriever.tot_introns).start()
    for n, gene in enumerate(intron_retriever):
        p_bar.update(n)
        gene_chunks = []
        for chunk in gene:
            gene_chunks.append(chunk)
            if n_probes(gene_chunks) > 200:
                break
        if chunk.id in used_probes: continue
        # if len(list(filtered_probe_table.find(target=chunk.id))) >= 24:
        #     continue
        # Tagging of introns in db should be added around here
        chunks.append(gene_chunks)
        if len(chunks) == chunksize:
            with NamedTemporaryFile("w") as fasta_input:
                # Break records into 1000nt chunks
                for gene_chunks in chunks:
                    # Reduce potential probeset size to drop design time
                    while n_probes(gene_chunks, probe_size) > 200:
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
