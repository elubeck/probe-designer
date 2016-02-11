import copy
import csv
import sys
from collections import defaultdict
from itertools import groupby
from pathlib import Path

import dataset
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from progressbar import ProgressBar
from probe_designer.sequence_getter import SequenceGetter
from probe_designer.utils.timeout import timeout
from probe_designer.utils.misc import gc_count

csv.field_size_limit(sys.maxsize)  # Prevent field size overflow.


def sub_seq_splitter(seq,
                     size,
                     gc_target=0.55,
                     stride=1,
                     spacing=2,
                     gc_min=None,
                     gc_max=None,
                     debug=True):
    """
    Takes sequence and designs best probes of given size closests to gc_target.
    """
    assert (isinstance(seq, str))
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
    seq_lookup = {probe['seq']: probe for probe in probes}
    overlapping = defaultdict(list)
    for probe_n, probe_ind in enumerate(probe_indices):
        for i in range(probe_n, len(probe_indices)):
            if i != probe_n:
                if any(probe_ind & probe_indices[i]):
                    overlapping[probes[probe_n]['seq']].append(probes[i])
    overlapping_b = copy.deepcopy(overlapping)
    for k, v in overlapping.items():
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
                if not any(om) or not any(om[0][1]):
                    break
                # if len(om[0][1]) == 0:
                #     break
                # Otherwise throw out worst probe
                else:
                    probe_group = [probe_seq for probe_seq, overind in om[1:]]
                if not any(probe_group):  # Is this correct?
                    break
                n_iters += 1
            probe_group = [probe_seq for probe_seq, overind in om[1:]]
        for probe_seq in probe_group:
            passed_probes.append(probe_seq)
            for overlapped in overlapping_2[probe_seq]:
                disallowed_probes.append(overlapped['seq'])
        probe_index += 1

    return passed_probes


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
                    for probe in sub_seq_splitter(str(chunk),
                                                  35,
                                                  gc_min=0.35, ):
                        gene_probes[gene].append(probe)
                    flat_p = [{'target': gene,
                               'seq': probe} for probe in gene_probes[gene]]
                    p_table.insert_many(flat_p)
        except TimeoutError:
            print("Timedout at {}".format(gene))
    pbar.finish()
    return gene_probes


def design_step(gene, max_time=180, cds_only=False, length=35, debug=True):
    rr = RNARetriever2()
    try:
        with timeout(seconds=max_time):
            try:
                gene_records = [
                    str(record.seq)
                    for record in rr.get_single_rna(gene,
                                                    chunk_size=100000,
                                                    cds_only=cds_only)
                ]
            except:
                print("Record lookup fail at {}".format(gene))
                return ("FAILED", [], [])
    except TimeoutError:
        print("Timedout at {}".format(gene))
        return ("FAILED", [], [])
    probes = []
    for chunk in gene_records:
        for probe in sub_seq_splitter(str(chunk),
                                      length,
                                      gc_min=0.35,
                                      debug=debug):
            probes.append(probe)
    return gene, probes, gene_records


def batch_design2(genes,
                  max_time=180,
                  cds_only=False,
                  db_name='mrna_probes_full_tx',
                  length=35,
                  debug=True):
    db = dataset.connect("sqlite:///db/{}.db".format(db_name))
    p_table = db['unfiltered_probes']
    cds_table = db['cds_table']
    gene_probes = defaultdict(list)
    pbar = ProgressBar(maxval=len(genes)).start()
    used = {g['target'] for g in p_table.distinct('target')}
    rr = RNARetriever2()
    gi = (gene for gene in genes if gene not in used)
    gene_probes = {}
    from multiprocessing import Pool, cpu_count
    from functools import partial
    d_wrap = partial(design_step,
                     max_time=max_time,
                     cds_only=cds_only,
                     length=length,
                     debug=debug)
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


class RNAGetter(SequenceGetter):
    def get_rna(self, gene_isoforms, cds_only=False):
        """
        Get all rnas from a gene and return in transcribed order.
        Only chooses exons in every isoform.  If cds_only is true,
        only choose regions from CDS.
        """

        def range_unpacker(vals):
            v0, v1 = sorted(vals)
            return set(range(v0 + 1, v1 + 1))

        rna_set = []
        exon_unions = []
        for isoform in gene_isoforms:
            tx_start, tx_end, e_start, e_end, cds_start, cds_end = self.get_positions(
                isoform)
            if cds_only is True:
                span = range_unpacker([cds_end, cds_start])
            else:
                span = range_unpacker([tx_start, tx_end])
            exon_pos = [
                p for start, end in zip(e_start, e_end)
                for p in range_unpacker([start, end]) if p in span
            ]
            rna_set.append(exon_pos)
            exon_unions.append(
                {tuple(sorted(pos))
                 for pos in zip(e_end[:-1], e_start[1:])})

        # Get Regions in Every Isoform
        conserved_regions = sorted(set.intersection(*(set(v) for v in rna_set
                                                      )))
        all_unions = set.intersection(*exon_unions)
        contigs = sorted(self.get_range(conserved_regions))
        union_lookup = dict((k, v + 1) for k, v in all_unions)

        # Combine regions that are seperated by an intron but contiguous for RNA
        i_frags = iter(contigs)

        frags = [next(i_frags, None)]
        while True:
            reg = next(i_frags, None)
            if reg is None:
                break
            if frags[-1][-1] in union_lookup:
                first_ind = union_lookup[frags[-1][-1]]
                if first_ind == reg[0]:
                    frags[-1] = frags[-1] + reg
            else:
                frags.append(reg)
        rnas = sorted(sorted(zip(r[0::2], r[1::2])) for r in frags)

        # Remove duplicate rnas from spliceoforms and return sorted result
        if isoform['strand'] == "+":
            return rnas
        elif isoform['strand'] == '-':
            return [r[::-1] for r in rnas][::-1]

    def run_gene(self, gene, cds_only=False):
        """
        Wrapper for finding results for a gene
        """
        hits = list(self.table.find(name2=gene))
        return {
            'chrom': hits[0]['chrom'],
            'strand': hits[0]['strand'],
            'name2': gene,
            'exons': self.get_rna(hits, cds_only)
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


class RNARetriever2(object):
    def get_single_rna(self, gene, chunk_size=1000, cds_only=False):
        row = self.rna_getter.run_gene(gene, cds_only)
        chrom_code = row['chrom']
        chrom_path = self.chromo_folder.joinpath(
            "{}.fa.masked".format(chrom_code))
        chrom_seq = SeqIO.read(str(chrom_path), 'fasta')
        return self.get_rna(row, chrom_seq, chunk_size)

    def get_rna(self, row, chrom_seq, chunk_size=10000000):
        name2 = row['name2']
        if not row['exons']:
            yield None
            return
        for exon in row['exons']:
            exon_seq = ""
            for start, end in exon:
                seq = chrom_seq[start:end]
                if row['strand'] == "+":  #TODO: CHECK THAT THIS IS CORRECT
                    seq = seq.reverse_complement()
                exon_seq += seq.seq
                # Drop masked sequences
            for sub_seq in exon_seq.split("N"):
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
