from __future__ import print_function, division
import dataset
from operator import itemgetter
from itertools import groupby
from collections import defaultdict
import json


def get_exon_pos():
    """
    Gets position of all exons from table
    """
    db = dataset.connect("sqlite:///db/refGene.db")
    table = db['mouse']
    exon_space = defaultdict(set)
    for line in table:
        chrom = line['chrom']
        e_start = map(int, line['exonStarts'].strip(',').split(","))
        e_end = map(int, line['exonEnds'].strip(',').split(","))
        for start, end in zip(e_start, e_end):
            if end <= start: raise Exception()
            for i in range(start, end + 1):
                exon_space[chrom].add(i)
    return exon_space


class IntronGetter(object):
    """
    Gets introns from refgenome table
    """

    def get_range(self, sequence, min_size=20):
        #Find longest contiguous sequences
        ranges = []
        for k, g in groupby(enumerate(sequence), lambda (i, x): i - x):
            group = map(itemgetter(1), g)
            start = group[0]
            end = group[-1]
            if end - start > min_size:
                ranges.append((start, end))
        return ranges

    def get_positions(self, gene):
        """
        Finds position of exons for gene
        """
        if gene['strand'] == '+':
            start = int(gene['txStart'])
            end = int(gene['txEnd'])
            e_start = map(int, gene['exonStarts'].strip(',').split(","))
            e_end = map(int, gene['exonEnds'].strip(',').split(","))
        elif gene['strand'] == '-':
            # Flip things around
            start = int(gene['txEnd'])
            end = int(gene['txStart'])
            e_end = map(int, gene['exonStarts'].strip(',').split(","))[::-1]
            e_start = map(int, gene['exonEnds'].strip(',').split(","))[::-1]
        return (start, end, e_start, e_end)

    def get_intron(self, gene_isoforms, in_all=False):
        rna_set = []
        for gene in gene_isoforms:
            to_flatten = []
            tx_start, tx_end, e_start, e_end = self.get_positions(gene)
            to_flatten.append((tx_start, e_start[0]))  # Add first intron
            to_flatten.append((e_end[-1], tx_end))  # Add last intron
            to_flatten += list(zip(e_end[:-1],
                                   e_start[1:]))  # Add all remaining introns
            intron_pos = {
                p
                for end, start in to_flatten for p in range(*sorted([end,
                                                                     start]))
            }
            rna_set.append(intron_pos)

        # Get Regions in Every Isoform
        conserved_regions = set.intersection(*rna_set)

        # Make set of points that are only intron1, not any other exons
        exon_pos = self.exons[gene['chrom']
                              ]  # Get all the exons for chromosome
        non_overlapping = sorted(conserved_regions.difference(exon_pos))

        rna_regions = [contiguous_frag
                       for contiguous_frag in self.get_range(non_overlapping)]
        # Remove duplicate rnas from spliceoforms and return sorted result
        rnas = sorted(set(rna_regions), key=lambda x: x[0])
        if gene['strand'] == "+":
            return rnas
        elif gene['strand'] == '-':
            return rnas[::-1]

    def get_intron2(self, gene_isoforms, in_all=False):
        """
        Get all introns from a gene and return in transcribed order.
        Removes regions that overlap exons
        If in_all=True makes sure selected introns are in all isoforms.
        """
        intron_set = []
        tx_starts = []
        tx_ends = []
        for gene in gene_isoforms:
            gene_introns = []
            tx_start, tx_end, e_start, e_end = self.get_positions(gene)
            tx_starts.append(tx_start)
            tx_ends.append(tx_end)
            # if gap between start and 1st exon
            gene_introns.append(sorted([tx_start, e_start[0] - 1]))
            # Iterate all middle introns
            for end, start in zip(e_end[:-1], e_start[1:]):
                gene_introns.append(sorted([end + 1, start - 1]))
            # Add last exon
            gene_introns.append(sorted([e_end[0], tx_end]))
            intron_set.append(
                [{'start': start,
                  'end': end,
                  'length': end - start} for start, end in gene_introns
                 if end - start != 0])
        if in_all:
            if gene['strand'] == "+":
                intron_set = [[intron for intron in isoform
                               if intron['start'] >= max(tx_starts)
                               if intron['end'] <= min(tx_ends)]
                              for isoform in intron_set]
            if gene['strand'] == '-':
                intron_set = [[intron for intron in isoform
                               if intron['start'] <= min(tx_starts)
                               if intron['end'] >= max(tx_ends)
                               if intron['end'] <= min(tx_starts)]
                              for isoform in intron_set]
        # Make set of points that are only intron1, not any other exons
        exon_pos = self.exons[gene['chrom']
                              ]  # Get all the exons for chromosome
        isoform_regions = []
        for isoform in intron_set:
            reg = []
            for intron in isoform:
                intron_range = set(range(intron['start'], intron['end']))
                non_overlapping = sorted(intron_range.difference(exon_pos))
                for r_vals in self.get_range(non_overlapping):
                    for p in range(*r_vals):
                        reg.append(p)
            isoform_regions.append(set(reg))
        overlapping_isoforms = set.intersection(*isoform_regions)
        intron_regions = [reg for reg in self.get_range(overlapping_isoforms)]
        import ipdb
        ipdb.set_trace()
        # Remove duplicate introns from spliceoforms and return sorted result
        introns = sorted(set(intron_regions), key=lambda x: x[0])
        if gene['strand'] == "+":
            return introns
        elif gene['strand'] == '-':
            return introns[::-1]

    def run_gene(self, gene):
        """
        Wrapper for finding results for a gene
        """
        hits = list(self.table.find(name2=gene))
        return {
            'chrom': hits[0]['chrom'],
            'strand': hits[0]['strand'],
            'name2': gene,
            'introns': self.get_intron(hits, in_all = self.in_all)
        }

    def __init__(self, table=None, in_all=False):
        if table is None:
            db = dataset.connect("sqlite:///db/refGene.db")
            self.table = db['mouse']
        else:
            self.table = table
        self.tot_records = len(list(self.table.distinct('name2')))
        self.in_all = in_all
        with open('db/exon_pos.json', 'r') as f_in:
            self.exons = json.load(f_in)


class IntronIterator(IntronGetter):
    def __iter__(self):
        for chromosome in self.table.distinct('chrom'):
            chrom_code = chromosome['chrom']
            for row in self.table.distinct('name2', chrom=chrom_code):
                name2 = row['name2']
                yield self.run_gene(name2)

    def __init__(self, ):
        super(self.__class__, self).__init__()


class ChromIntronIterator(IntronGetter):
    def __iter__(self):
        for row in self.table.distinct('name2', chrom=self.chrom):
            name2 = row['name2']
            if name2 in self.skip_introns: continue
            yield self.run_gene(name2)

    def __init__(self, chrom, skip_introns=None):
        super(self.__class__, self).__init__()
        self.chrom = chrom
        if skip_introns is None:
            self.skip_introns = []
        else:
            self.skip_introns = skip_introns
