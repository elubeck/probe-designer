from __future__ import division, print_function

import csv
import gzip
import random
import sys
from collections import Counter
from math import ceil
from itertools import groupby

import dataset
from Bio import SeqIO
from Bio.Seq import Seq
from tempfile import NamedTemporaryFile
from progressbar import ProgressBar
import blaster2

from oligoarray_designer import Oligoarray, OligoarrayDesigner

csv.field_size_limit(sys.maxsize) # Prevent field size overflow.

def gc_count(probe):
    return len([1 for c in probe.lower() if c in ['c', 'g']]) / len(probe)


zak_genes = [
    'Albhk5', 'Ash1', 'Axin2', 'Bmi1', 'bmp4', 'bmpr1a', 'Brachyury', 'Cdh1',
    'Col5a2', 'ctcf', 'dazl', 'Dnmt1', 'dnmt3a', 'Dnmt3b', 'Dnmt3L', 'Dppa2',
    'Dppa3', 'Dppa4', 'EED', 'ehmt2', 'Esrrb', 'Ezh2', 'Fbxo15',
    'Fgf4', 'FGF5', 'Fgfr2', 'Foxa2', 'FoxO1', 'pax6', 'Pecam1', 'Pou5f1',
    'Prdm14', 'Rest', 'Sall4', 'Sdha',
    'setdb1', 'Smad1', 'Smad4', 'Smad5', 'Socs3', 'Sox2',
    'Sp1', 'Suz12', 'sycp3', 'Tbx3', 'Tcl1', 'Tet1', 'tet2', 'tet3', 'Thy1',
    'Trim28', 'Utf1', 'Xist', 'Zfp42', 'Mecp2',
    'MBD3', 'Foxa2', 'Gata4', 'Ptf1a', 'Cdx2', 'Eomes', 'Gcm1', 'Krt1', 'Afp',
    'Serpina1a', 'Fn1', 'Lama1', 'Lamb1', 'Lamc1', 'Sox17', 'T', 'Wt1', 'Des',
    'Myf5', 'Myod1', 'Hba-x', 'Hbb-y', 'Col1a1', 'Runx2', 'Nes', 'Neurod1',
    'Pax6', 'Cd34', 'Cdh5', 'Flt1', 'Pecam1', 'Ddx4', 'Sycp3', 'Gcg', 'Iapp',
    'Ins2', 'Pax4', 'Pdx1', 'Sst', 'Olig2', 'Tat', 'Foxd3', 'Gata6',
    'Gbx2', 'Nanog', 'Nr5a2', 'Nr6a1', 'Pou5f1', 'Sox2', 'Tcfcp2l1',
    'Utf1', 'Zfp42', 'Commd3', 'Crabp2', 'Ednrb', 'Fgf4', 'Fgf5', 'Gabrb3',
    'Gal', 'Grb7', 'Hck', 'Ifitm1', 'Il6st', 'Kit', 'Lefty1', 'Lefty2',
    'Lifr', 'Nodal', 'Nog', 'Numb', 'Pten', 'Sfrp2', 'Tdgf1', 'Fgf4', 'Fgf5',
    'Gdf3', 'Lefty1', 'Lefty2', 'Nodal', 'Brix1', 'Cd9', 'Diap2', 'Ifitm2',
    'Igf2bp2', 'Lin28a', 'Podxl', 'Rest', 'Sema3a', 'Tert', 'Il6st', 'Lifr',
    'Socs3', 'Stat3', 'Fgfr1', 'Fgfr2', 'Fgfr3', 'Fgfr4', 'Ptch1', 'Ptchd2',
    'Smo', 'Gli1', 'Gli2', 'Gli3', 'Sufu', 'Ncstn', 'Notch1', 'Notch2',
    'Notch3', 'Notch4', 'Psen1', 'Psen2', 'Psenen', 'Rbpjl', 'Hes1', 'Hey1',
    'Acvr1', 'Acvr1b', 'Acvr1c', 'Acvr2a', 'Acvr2b', 'Acvrl1', 'Amhr2',
    'Bmpr1a', 'Bmpr1b', 'Bmpr2', 'Eng', 'Ltbp1', 'Ltbp2', 'Ltbp3', 'Ltbp4',
    'Rgma', 'Tgfbr1', 'Tgfbr2', 'Tgfbr3', 'Tgfbrap1', 'Crebbp', 'E2f5',
    'Ep300', 'Rbl1', 'Rbl2', 'Smad1', 'Smad2', 'Smad3', 'Smad4', 'Smad5',
    'Smad6', 'Smad7', 'Smad9', 'Sp1', 'Zeb2', 'Fzd1', 'Fzd2', 'Fzd3', 'Fzd4',
    'Fzd5', 'Fzd6', 'Fzd7', 'Fzd8', 'Fzd9', 'Lrp5', 'Lrp6', 'Vangl2', 'Bcl9',
    'Bcl9l', 'Ctnnb1', 'Lef1', 'Nfat5', 'Nfatc1', 'Nfatc2', 'Nfatc3', 'Nfatc4',
    'Pygo2', 'Tcf7l1', 'Tcf7', 'Tcf7l2', 'Mecp2', 'MBD3'
]

zak_genes = map(lambda x: x.lower(), zak_genes)

from Bio.SeqRecord import SeqRecord
from pathlib import Path


class IntronRetriever(object):
    def get_single_intron(self, gene, chunk_size=1000):
        row = self.table.find_one(name2=gene)
        chrom_code = row['chrom']
        chrom_path = self.chromo_folder.joinpath(
                        "{}.fa.masked".format(chrom_code))
        chrom_seq = SeqIO.read(str(chrom_path), 'fasta')
        return self.get_intron(row, chrom_seq, chunk_size)
        
    def get_intron(self, row, chrom_seq, chunk_size=1000):
        name2 = row['name2']
        coords_str = row['coords'].strip("[]")
        # Intron doesn't have assigned coords
        if not coords_str:
            yield None
            return
        for coord_pair in coords_str.split(", ("):
            start, end = map(int, coord_pair.strip("()").split(", "))
            seq = chrom_seq[start:end + 1]
            if row['strand'] == "+": #TODO: CHECK THAT THIS IS CORRECT
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
        for chromosome in self.table.distinct('chrom'):
            chrom_code = chromosome['chrom']
            chrom_path = self.chromo_folder.joinpath(
                "{}.fa.masked".format(chrom_code))
            chrom_seq = SeqIO.read(str(chrom_path), 'fasta')
            for row in self.table.find(chrom=chrom_code):
                gene_records = self.get_intron(row, chrom_seq)
                records = list(gene_records)
                if any(records):
                    yield records
                # for intron_chunk in gene_records:
                #     yield intron_chunk

    def __init__(self, organism='mouse'):
        self.organism = organism
        self.db = dataset.connect("sqlite:///db/intron1_coords.db")
        self.table = self.db[self.organism]
        self.chromo_folder = Path("db/chromFaMasked/")
        self.tot_introns = len(list(self.table.distinct("name2")))


def n_probes(chunk_list, probe_size=35):
    # Get # of probes that could be made from set of gene chunks
    return sum(len(chunk) // probe_size for chunk in chunk_list)

class Designer(object):

    def run_oligoarray(self, chunks, max_probes=150, probe_size=None):
        """
        Chunks is a list of gene fragments lists. [Gene[fragment_1...fragment_m],...Gene_n]
        """
        assert(isinstance(chunks[0], list))
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
        used_probes = set([row['Name'] for row in self.probe_db.distinct("Name")])
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
        intron_db = dataset.connect("sqlite:///db/intron_probes2.db")
        self.probe_db = intron_db[organism]
        self.o = OligoarrayDesigner(
            blast_db ='/home/eric/blastdb/old_format/transcribed_mouse2')
        self.probe_size = probe_size

def design_introns():
    # First get used probes
    intron_db = dataset.connect("sqlite:///db/intron_probes2.db")
    probe_db = intron_db['mouse']
    used_probes = set([row['Name'] for row in probe_db.distinct("Name")])
    o = OligoarrayDesigner(
        blast_db='/home/eric/blastdb/old_format/transcribed_mouse2')
    chunksize = 12
    probe_size = 35
    chunks = []
    intron_retriever = IntronRetriever()
    p_bar = ProgressBar(maxval=intron_retriever.tot_introns).start()
    for n, gene in enumerate(intron_retriever):
        if gene[0].id in used_probes:
            continue
        name = gene[0].id.lower()
        chunks.append(gene)
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

# design_introns()
# chunks = IntronRetriever().get_single_intron("Pgk1")
# o = OligoarrayDesigner(
#     blast_db='/home/eric/blastdb/old_format/transcribed_mouse2')
# probe_size = 35
# with NamedTemporaryFile("w") as fasta_input:
#     # Break records into 1000nt chunks
#     for gene_chunks in chunks:
#         if n_probes(gene_chunks) > 150:
#             # Reduce potential probeset size to drop design time
#             while n_probes(gene_chunks, probe_size) > 150:
#                 gene_chunks = random.sample(gene_chunks,
#                                             len(gene_chunks) - 1)
#         SeqIO.write(gene_chunks, fasta_input,
#                     'fasta')  # write all sequences to file
#     fasta_input.flush()
#     res = o.run(fasta_input.name,
#                 max_dist=1000,
#                 min_length=probe_size,
#                 max_length=probe_size,
#                 max_tm=100,
#                 cross_hyb_temp=72,
#                 min_tm=74,
#                 secondary_struct_temp=76,
#                 max_oligos=100,
#                 timeout=15)
# for name, probes in res:
#     probes['Percent GC'] = 100 * probes["Probe (5'-> 3')"].map(
#         gc_count)

# import mRNA_designer
# flat_probes = probes.T.to_dict().values()

class ProbeFilter(object):
    def run_batch(self):
        hit_val2 = sorted(hit_vals, key=lambda x: x[0].split(',')[0])
        finished_probes = {}
        # Iterate through every probeset
        for target, probe_set in groupby(hit_val2,
                                            key=lambda x: x[0].split(',')[0]):
            pass

    def filter_probes(self, probe_set, target, hits, probe_lookup, max_hits=10000, off_target_thresh=7, probe_num=48,):
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
        while len(passed_probes) < probe_num and choice_index != len(remaining_probes):
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

        
    
    def blast2copynum(self, hits):
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
                hit_vals.append((probe_name, self.get_copynum(gencode_id)))
            except:
                print("Failed @ {}".format(gencode_id))
        hit_vals = sorted(hit_vals, key=lambda x: x[1], reverse=True)
        return hit_vals

    def run_blast(self, probe_lookup, match_thresh, strand=None):
        if not strand:
            strand = self.strand
        fasta_str = "\n".join(">{}\n{}".format(*items) for items in probe_lookup.iteritems())
        res = blaster2.local_blast_query(fasta_str, db='gencode_tracks_reversed')
        hits = blaster2.parse_hits(res, strand=strand, match_thresh=18)
        return hits

    def run(self, flat_probes, gene_name, match_thresh=18, n_probes=48, max_off_target=10000):
        # probe_lookup = {
        #     "{Name},{Probe #}".format(**probe): probe["Probe (5'-> 3')"]
        #     for probe in flat_probes}
        probe_lookup = {"{},{}".format(gene_name, n): probe for n, probe in enumerate(flat_probes)}
        res = self.run_blast(probe_lookup, match_thresh)
        hit_vals = self.blast2copynum(res)
        filtered_probes = self.filter_probes(hit_vals, gene_name, res, probe_lookup, probe_num=n_probes, max_hits=max_off_target)
        finished_probes = [probe_lookup[probe]
                                    for probe in filtered_probes]
        return finished_probes

    def get_copynum(self, hits):
        off_target = {name: self.counts[name] for name in hits if name in self.counts.keys()}
        false_hits = sum(off_target.values())
        return false_hits

    def __init__(self, db='gencode_tracks_reversed', strand='+', copy_num='embryonic11.5'):
        self.db = db
        if strand == '+':
            self.strand = 1
        else:
            self.strand = -1
        if copy_num == 'embryonic11.5':
            # Get merged embryonic 11.5 encode data
            with open('db/encode_counts.csv', 'r') as f:
                self.counts = {line[0]: float(line[1]) for line in csv.reader(f)}
        elif copy_num == 'brain':
            with open('db/brain_counts.csv', 'r') as f:
                self.counts = {line[0]: float(line[1]) for line in csv.reader(f)}
        else:
            raise Exception('need a valid copy_num database')


# intron_db = dataset.connect("sqlite:///db/intron_probes2.db")
# intron_db_filtered = dataset.connect("sqlite:///db/intron_probes_filtered.db")
# probe_db = intron_db['mouse']
# filtered_probe_table = intron_db_filtered['mouse']

# pf = ProbeFilter()
# p_num = ProgressBar(maxval=len(list(probe_db.distinct("Name")))).start()
# used = [p['target'] for p in filtered_probe_table.distinct("target")]
# for n, gene in enumerate(probe_db.distinct("Name")):
#     p_num.update(n)
#     if gene in used: continue
#     gene = gene['Name']
#     p_set = [probe["Probe (5'-> 3')"] for probe in probe_db.find(Name=gene)]
#     if len(p_set) < 20: continue
#     passed_probes = pf.run(p_set, gene, max_off_target=5000, n_probes=60)
#     for probe in passed_probes:
#         filtered_probe_table.insert({'target':gene, 'seq':probe})
# p_num.finish()


# import blaster2
# strand = 1
# res = blaster2.local_blast_query(fasta_str, db='gencode_tracks_reversed')
# hits = blaster2.parse_hits(res, strand=strand, match_thresh=18)

# hit_vals = []
# for probe_name, matches in hits.iteritems():
#     gene_name = probe_name.split(",")[0]
#     gencode_id, refseq = map(list, zip(*[match.split(',')
#                                             for match in matches]))
#     if not gene_name in refseq:
#         continue
#         bad_count += 1
#         # raise Exception("Query not found in hits")
#     del gencode_id[refseq.index(gene_name)]
#     try:
#         hit_vals.append((probe_name, blaster2.get_copynum(gencode_id)))
#     except:
#         print("Failed @ {}".format(gencode_id))
# hit_vals = sorted(hit_vals, key=lambda x: x[1], reverse=True)

# def gc_count(probe):
#     return len([1 for c in probe.lower() if c in ['c', 'g']]) / len(probe)

# # Iterate everyg ene
# from itertools import groupby
# from collections import Counter
# hit_val2 = sorted(hit_vals, key=lambda x: x[0].split(',')[0])
# off_target_thresh = 7
# finished_probes = {}
# # Iterate through every probeset
# for target, probe_set in groupby(hit_val2,
#                                     key=lambda x: x[0].split(',')[0]):
#     probe_set = list(probe_set)

#     # First pick probes with no off target hits
#     passed_probes = [probe for probe in probe_set if probe[1] == 0]
#     remaining_probes = set(probe_set)
#     remaining_probes.difference_update(passed_probes)

#     # Sort remaining probes by off target hits
#     choices = sorted(list(remaining_probes), key=lambda x: x[1])
#     choice_index = 0
#     off_target = Counter([])
#     probe_num = 48 

#     # Iterate until at least 24 probes are picked or no more probes are available
#     while len(passed_probes) < probe_num and choice_index != len(remaining_probes):
#         selected = choices[choice_index]
#         # Get off target hits for chosen probe
#         off_target_hit = [k for k in hits[selected[0]]
#                             if target != k.split(',')[1]]
#         test_counter = off_target + Counter(off_target_hit)

#         # Check if adding this probe makes anything go over off target threshold
#         over_thresh = [True for k in test_counter.values()
#                         if k >= off_target_thresh]
#         if not any(over_thresh):
#             passed_probes.append(selected)
#             off_target = test_counter
#         choice_index += 1

#     # If more than 24 probes chosen optimize on GC
#     if len(passed_probes) > probe_num:
#         # Get GC counts of every probe
#         probe_gc = [(probe[0], gc_count(probe_lookup[probe[0]]))
#                     for probe in passed_probes]
#         gc_target = 0.55
#         gc_range = 0.01
#         multiplier = 1
#         chosen_gc = []
#         # Get closest probe set to 0.55 gc
#         while len(chosen_gc) < probe_num:
#             gc_min = gc_target - gc_range * multiplier
#             gc_max = gc_target + gc_range * multiplier
#             chosen_gc = [probe for probe, gc in probe_gc
#                             if gc_max >= gc >= gc_min]
#             multiplier += 1

#         # If too many probes still exists choose a random subset
#         if len(chosen_gc) != probe_num:
#             chosen_gc = random.sample(chosen_gc, probe_num)
#         passed_probes = chosen_gc
#     else:
#         passed_probes = [probe for probe, n in passed_probes]
#     finished_probes[target] = [probe_lookup[probe]
#                                 for probe in passed_probes]

# with open('/home/eric/Downloads/PGK1_intron_48.json', 'w') as fout:
#     import json
#     json.dump(finished_probes, fout)

# {u'Pgk1': ['CCTAGAGATTTGCAGGATAAGTACTACGTCCATTC',
#   'AATGGCCCACTGTCTTGTCTTTGAGCAGCTGCATG',
#   'TGGTGTCTGCTCTCCGTTTTCCTCTTCCCACCAGA',
#   'TTGTAAAATCAACTAGAGTGAGGCGGAACCTGTTG',
#   'GAACAGGTAATAAGTAACGCCCCTTGATTCTAGTT',
#   'CCCACTGACAATAGATAAGTGGGAGGCATATGATC',
#   'ACAAAAGGGTTAACTGCACAGGGTGGATAGTGCTG',
#   'GGCTTCCTCTCTACCCTGAAGACCAATCTACCGCA',
#   'CACTGGTCTTGATTATGAAACTCCATCTTTCCCAT',
#   'CAAGCTGTGATGACATCAGCTTACCGATGAAAACT',
#   'GAAGAAAACGTACAAGAAAGAGAAGTGCCAGCAAG',
#   'CCATCCTTGGGTCCTATAGTTGTTACATGAACAAC',
#   'ATTAGTTCTCACTGACCTCTTGATTGCCTTTGACC',
#   'TGATGCAAACCCCAGGGCCAGGAGTCCAACCGCAA',
#   'TCGTCACGGATGGTCCAAGAATCTTGGGGGAAGGG',
#   'TGTTTAGTGTGATGGTCTCTGTAGCAGTCTGAGTG',
#   'CCATCCTGATTATTCTATCCCATTTGGGGTAACAT',
#   'AGAAAGGAGGTGAGTGAACGGGCCCAGGATGCAGG',
#   'CCTAAAGGACACATTTGGTCGCTTTTTAAGCTAGA',
#   'TGCCTTGTGGACTGATTGCTAAAGTCTGCTAGTAC',
#   'TGAAGCACTGCCAGCAACTACCACAAAAGGCACTG',
#   'CTTTTCCAGGCTTAAAGGTTAACCTGGGTTAAATG',
#   'CCCTCACAGGATGAATCTTTTGTCCTTAGAACAGA',
#   'TTTGGGCTAAAATATCCAACCATCTGTGGCTTCTA']}