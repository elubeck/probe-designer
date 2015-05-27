from __future__ import division, print_function, with_statement

import os
import re
import time
from collections import defaultdict
from itertools import groupby
from operator import itemgetter
import tempfile

import arrow
import dataset
import numpy as np
from Bio import AlignIO, Entrez, SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline, MuscleCommandline
from Bio.Application import ApplicationError
from future.builtins import map, object, range, str

Entrez.email = 'elubeck@caltech.edu'

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement[c] for c in seq.upper())[::-1]


class CDSError(Exception):
    def __init__(self, exception):
        Exception.__init__(self, exception)

class CDS(object):
    """
    CDS object for designing a gene.
    If variants = True will fine only the probes that hit all variants of the gene.  This prevents probes from favoring
    any one spliceoform.
    """

    def mRNA_cds(self):
        """
        Searches Entrez for mRNAs for a given reference.
        @param gene: Name of gene
        @param organism: A string to restrict entrez queries to.
        @return: Handle to results @raise Exception:
        """
        handle = None
        iterations = 0
        max_it = 10
        while True:
            handle1 = Entrez.esearch(db="nucleotide", term='(%s) AND %s AND mRNA' % (self.gene, self.organism),
                                    retmax=500)
            record = Entrez.read(handle1)
            handle = Entrez.efetch(db='nucleotide', id=record['IdList'], rettype='gb')
            if handle is not None:
                break
            if iterations > max_it:
                raise Exception("CDS not retrieved")
            iterations += 1
            time.sleep(4)
        res = list(SeqIO.parse(handle, format='gb'))
        if len(res) == 0:
            print("No Records found for %s" % self.gene)
        return res

    def align_seqs(self, seqs, aligner='muscle'):
        """
        Given a number of CDS variants for a given gene, aligns the variants to find conserved sequences.

        @param seqs: list of CDS variants as gi records
        @param aligner: program to use for alignment
        @return: list of CDS
        """
        if len(seqs) == 1:
            return [str(seqs[0].seq)]
        with tempfile.NamedTemporaryFile() as out_file:
            with tempfile.NamedTemporaryFile('w') as f_in:
                SeqIO.write(seqs, f_in, 'fasta')
                f_in.flush()
                if aligner == 'muscle':
                    cline = MuscleCommandline(input=f_in.name,
                                            out=out_file.name,)
                elif aligner == 'clustalo':
                    cline = ClustalOmegaCommandline(infile=f_in.name,
                                                    outfile=out_file.name,)
                try:
                    stdout, stderr = cline()
                except:
                    raise CDSError("Muscle failed for {}".format(self.gene))
            ali = AlignIO.read(out_file.name, format='fasta')
        ali_arr = np.array([list(rec) for rec in ali], np.character)
        letter_count = np.sum(ali_arr != '-',0)
        no_blanks = np.nonzero(letter_count==len(ali_arr))[0]
        passed = [i for i in no_blanks
                  if len(np.unique(ali_arr[:, i])) == 1]
        contiguous = []
        for k, g in groupby(enumerate(passed), lambda i_x:i_x[0]-i_x[1]):
            item = list(map(itemgetter(1), g))
            if len(item) >= 20:
                assert((ali_arr[0, item] == ali_arr[1, item]).all())
                contiguous.append(''.join(ali_arr[0, item]))
        return contiguous

    def filter_seqs(self, handle):
        """
        Filters the entrez query to return only matching gene records
        @param handle: list of gene records of several genes
        @return: list of gene records for same gene
        """
        reg_exp = re.compile('\(([^\)]+)\),')
        match = ['NM']
        if self.variants == True:
            match = ['NM', 'XM']
        passed = {}
        for record in handle:
            if any(record.id.startswith(cond) for cond in match):
                if self.variants is False and 'variant' in record.description:
                    continue
                passed[record.description] = record
        #Attempts to only find records that have the given gene name
        #Flattens nested list of grouped genes
        #gene_query =  {k.lower(): [v1 for k1, v1 in v]
        #               for k,v in groupby(passed.iteritems(),
        #                    lambda x: reg_exp.findall(x[0])[-1])}
        gene_query = defaultdict(list)
        for desc, record in list(passed.items()):
            gene_name = reg_exp.findall(desc)[-1]
            gene_query[gene_name.lower()].append(record)
        if len(gene_query) == 1:
            return list(list(gene_query.values())[0])
        elif len(gene_query) > 1:
            replaced_name = self.gene.lower().replace(".", "-")
            if self.gene.lower() in list(gene_query.keys()):
                return list(gene_query[self.gene.lower()])
            # This is a hacky fix.  Sometimes . is turned into -
            elif replaced_name in list(gene_query.keys()):
                return list(gene_query[replaced_name])
            else:
                raise CDSError("Too many genes returned from Entrez query for {}".format(self.gene))
        elif len(gene_query) == 0:
            raise CDSError("No genes found for {}".format(self.gene))
        raise Exception("This shouldn't have happened for {}".format(self.gene))

    def extract_feature(self, cds_list, feature_query='CDS'):
        """
        Extracts a feature from a list of gene records.  Defaults to extracting only CDS
        @param cds_list: list of gene records
        @param feature_query: the feature to look for.
        @return: an iterator that generates gi records for the given feature query.
        """
        for record in cds_list:
            for feature in record.features:
                if feature.type == feature_query:
                    yield feature.extract(record)


    def merge_cds(self, handle,):
        """
        @type self: object
        @param handle:
        @return: @raise Exception:
        """
        cds = self.filter_seqs(handle)
        features = list(self.extract_feature(cds, 'CDS'))
        if len(cds) > 1 and self.variants is False:
            raise Exception("Variants should not have been found.")
        if len(cds) > 1:
            contiguous = self.align_seqs(features)
        elif any(cds):
            contiguous = [str(cds[0].seq)]
        else:
            raise Exception("Failed finding CDS")
        return {"CDS List": contiguous,
                '# Isoforms': len(cds)}

    def check_db(self):
        res = self.table.find(gene=self.gene, variants=self.variants)
        for item in res:
            if arrow.get(item['date']) >= self.date.replace(years = -1):
                item['CDS List'] = item['CDS List'].split(",")
                return item
        else:
            return None

    def save_db(self, cds):
        cds['date'] = arrow.now().datetime
        cds['variants'] = self.variants
        cds['gene'] = self.gene
        cds['CDS List'] = ",".join(cds['CDS List'])
        self.table.insert(cds)


    def run(self):
        db_ref = self.check_db()
        if db_ref is not None:
            print("Getting Entry from DB for {}".format(self.gene))
            return db_ref
        handle = self.mRNA_cds()
        cds = self.merge_cds(handle,)
        self.save_db(cds)
        print("Storing Entry for {} in DB table for {}".format(self.gene, self.organism))
        return cds

    def __init__(self, gene, organism='"Mus musculus"[porgn:__txid10090]', variants=True):
        self.gene = gene
        self.organism = organism
        self.variants = variants
        self.db = dataset.connect("sqlite:///db/cds.db")
        self.table = self.db[self.organism]
        self.date = arrow.utcnow()

def dump_db_to_fasta(file_path, organism='"Mus musculus"[porgn:__txid10090]'):
    db = dataset.connect("sqlite:///db/cds.db")
    table = db[organism]
    with open(file_path, 'w') as f:
        for line in table.find():
            for n,sub_seq in enumerate(line['CDS List'].split(",")):
                fasta_str = ">{}$Block#{}\n{}\n".format(line['gene'],n, sub_seq)
                f.write(fasta_str)
    print("DONE")
        
# dump_db_to_fasta("cds.fasta")            
if __name__ == '__main__':
    pass
    # cds = CDS("dazl").run()
    # print(len(''.join(cds['CDS List'])))
