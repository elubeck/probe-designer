from builtins import object, sorted, enumerate, int, len, list
from itertools import groupby
from operator import itemgetter
import dataset

__author__ = 'Eric'


class SequenceGetter(object):
    """
    Gets sequences from refgenome table
    """

    def get_range(self, sequence, min_size=20):
        #Find longest contiguous sequences
        sequence = sorted(sequence)
        ranges = []
        for k, g in groupby(enumerate(sequence), lambda x: x[0] - x[1]):
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
            cds_start = int(gene['cdsStart'])
            cds_end = int(gene['cdsEnd'])
        elif gene['strand'] == '-':
            # Flip things around
            start = int(gene['txEnd'])
            end = int(gene['txStart'])
            e_start = splitter(gene['exonStarts'])[::-1]
            e_end = splitter(gene['exonEnds'])[::-1]
            cds_end = int(gene['cdsStart'])
            cds_start = int(gene['cdsEnd'])
        return (start, end, e_start, e_end, cds_start, cds_end)

    def __init__(self, table=None, ):
        if table is None:
            db = dataset.connect("sqlite:///db/refGene.db")
            self.table = db['mouse']
        else:
            self.table = table
        self.tot_records = len(list(self.table.distinct('name2')))
