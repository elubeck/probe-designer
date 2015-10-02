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


from sequence_getter import SequenceGetter
class IntronGetter(SequenceGetter):
    """
    Gets introns from refgenome table
    """

    def get_intron(self, gene_isoforms, in_all=False):
        rna_set = []
        for gene in gene_isoforms:
            to_flatten = []
            tx_start, tx_end, e_start, e_end, cds_start, cds_end = self.get_positions(gene)
            to_flatten.append((tx_start, e_start[0]))  # Add first intron
            to_flatten.append((e_end[-1], tx_end))  # Add last intron
            to_flatten += list(zip(e_end[:-1],
                                   e_start[1:]))  # Add all remaining introns
            intron_pos = {
                p
                for end, start in to_flatten
                for p in range(*sorted([end + 1, start + 1]))
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
        super(IntronGetter, self).__init__(table=table)
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
