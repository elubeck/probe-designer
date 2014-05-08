from __future__ import print_function, with_statement, division
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import Entrez
from Bio import SeqIO
from Bio import AlignIO
import numpy as np
from operator import itemgetter
from itertools import groupby
import os

Entrez.email = 'elubeck@caltech.edu'


def get_mRNA_cdss(gene, organism='"Mus musculus"[porgn:__txid10090]', conserved_region=True):
    handle = Entrez.esearch(db="nucleotide", term='(%s) AND %s AND mRNA' % (gene, organism))
    record = Entrez.read(handle)
    handle = Entrez.efetch(db='nucleotide', id=record['IdList'], rettype='gb')
    return handle


def align_seqs(seqs):
    in_file = "temp_variant_cds.fasta"
    out_file = "temp_variant_cds_aligned.fasta"
    with open(in_file, 'w') as f:
        SeqIO.write(seqs, f, 'fasta')
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


def merge_cds(handle, variants=True):
    cds = []
    match = ['NM']
    if variants == True:
        match = ['NM', 'XM']
    for record in SeqIO.parse(handle, format='gb'):
        if any(record.id.startswith(cond) for cond in match):
            if variants is False and 'variant' in record.id:
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
