from __future__ import print_function, with_statement, division
from Bio.Align.Applications import ClustalOmegaCommandline, MuscleCommandline
from Bio import Entrez
from Bio import SeqIO
from Bio import AlignIO
import numpy as np
from operator import itemgetter
from itertools import groupby
import os
import time

Entrez.email = 'elubeck@caltech.edu'


def get_mRNA_cdss(gene, organism='"Mus musculus"[porgn:__txid10090]', conserved_region=True):
    handle = None
    iterations = 0
    max_it = 10
    while True:
        handle = Entrez.esearch(db="nucleotide", term='(%s) AND %s AND mRNA' % (gene, organism))
        record = Entrez.read(handle)
        handle = Entrez.efetch(db='nucleotide', id=record['IdList'], rettype='gb')
        if handle is not None:
            break
        if iterations > max_it:
            raise Exception("CDS not retrieved")
        iterations += 1
        time.sleep(4)
    return handle


def align_seqs(seqs, aligner='muscle'):
    import random
    r_num = ''.join(map(str, random.sample(range(100), 10)))
    in_file = "%s_temp_variant_cds.fasta" % r_num
    out_file = "%s_temp_variant_cds_aligned.fasta" % r_num
    with open(in_file, 'w') as f:
        SeqIO.write(seqs, f, 'fasta')
    if aligner == 'muscle':
        cline = MuscleCommandline(input=in_file,
                                  out=out_file,)
    elif aligner == 'clustalo':
        cline = ClustalOmegaCommandline(infile=in_file,
                                        outfile=out_file,)
    cline()
    ali = AlignIO.read(out_file, format='fasta')
    ali_arr = np.array([list(rec) for rec in ali], np.character)
    letter_count = np.sum(ali_arr != '-',0)
    no_blanks = np.nonzero(letter_count==len(ali_arr))[0]
    passed = [i for i in no_blanks
              if len(np.unique(ali_arr[:, i])) == 1]
    contiguous = []
    for k, g in groupby(enumerate(passed), lambda (i,x):i-x):
        item = map(itemgetter(1), g)
        if len(item) >= 20:
            assert((ali_arr[0, item] == ali_arr[1, item]).all())
            contiguous.append(''.join(ali_arr[0, item]))
    os.remove(in_file)
    os.remove(out_file)
    return contiguous


def merge_cds(handle, variants=True, query=None):
    cds = []
    match = ['NM']
    if variants == True:
        match = ['NM', 'XM']
    # for record in SeqIO.parse(handle, format='gb'):
    for record in handle:
        if any(record.id.startswith(cond) for cond in match):
            if variants is False and 'variant' in record.description:
                continue
            if query is not None and "(%s)" % query.lower() not in record.description.lower():
                continue
            for feature in record.features:
                if feature.type=='CDS':
                    cds.append(feature.extract(record))
                    break
    if len(cds) > 1 and variants is False:
        raise Exception("Variants should not have been found.")
    if len(cds) > 1:
        contiguous = align_seqs(cds)
    elif any(cds):
        contiguous = [str(cds[0].seq)]
    else:
        raise Exception("Failed finding CDS")
    return {"CDS List": contiguous,
            '# Isoforms': len(cds)}


if __name__ == '__main__':
    gene = "RUFY4"
    handle = get_mRNA_cdss(gene)
    cdss = merge_cds(handle)
