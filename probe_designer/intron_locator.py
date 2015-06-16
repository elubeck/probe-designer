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

    def get_intron(self, gene_isoforms):
        """
        Get all introns from a gene and return in transcribed order.
        Removes regions that overlap exons
        """
        intron_set = []
        for gene in gene_isoforms:
            gene_introns = []
            tx_start, tx_end, e_start, e_end = self.get_positions(gene)
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

        # Make set of points that are only intron1, not any other exons
        exon_pos = self.exons[gene['chrom']]  # Get all the exons for chromosome
        intron_regions = []
        for isoform in intron_set:
            for intron in isoform:
                intron_range = set(range(intron['start'], intron['end']))
                non_overlapping = sorted(intron_range.difference(exon_pos))
                for contiguous_frag in self.get_range(non_overlapping):
                    intron_regions.append(contiguous_frag)
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
            'introns': self.get_intron(hits)
        }

    def __init__(self, table=None):
        if table is None:
            db = dataset.connect("sqlite:///db/refGene.db")
            self.table = db['mouse']
        else:
            self.table = table
        self.tot_records = len(list(self.table.distinct('name2')))
        with open('db/exon_pos.json', 'r') as f_in:
            self.exons = json.load(f_in)


class IntronIterator(IntronGetter):
    def __iter__(self):
        for chromosome in self.table.distinct('chrom'):
            chrom_code = chromosome['chrom']
            for row in self.table.distinct('name2', chrom=chrom_code):
                name2 = row['name2']
                yield self.run_gene(name2)

    def __init__(self, chrom):
        super(self.__class__, self).__init__()
        self.chrom = chrom
