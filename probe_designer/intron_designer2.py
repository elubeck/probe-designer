from __future__ import division, print_function

import csv
import gzip
import random
from collections import Counter
from math import ceil

import dataset
from Bio import SeqIO
from Bio.Seq import Seq

from oligoarray_designer import Oligoarray, OligoarrayDesigner


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
    def get_intron(self, row, chrom_seq):
        name2 = row['name2']
        coords_str = row['coords'].strip("[]")
        # Intron doesn't have assigned coords
        if not coords_str:
            yield None
            return
        for coord_pair in coords_str.split(", ("):
            start, end = map(int, coord_pair.strip("()").split(", "))
            seq = chrom_seq[start:end + 1]
            if row['strand'] == "-":
                seq = seq.reverse_complement()
            # Drop masked sequences
            for sub_seq in seq.seq.split("N"):
                if len(sub_seq) >= 20:
                    record = SeqRecord(sub_seq,
                                       name=name2,
                                       id=name2,
                                       description='')
                    yield record

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


if __name__ == "__main__":
    # First get used probes
    intron_db = dataset.connect("sqlite:///db/intron_probes.db")
    probe_db = intron_db['mouse']
    used_probes = set([row['Name'] for row in probe_db.distinct("Name")])
    from tempfile import NamedTemporaryFile
    o = OligoarrayDesigner(
        blast_db='/home/eric/blastdb/old_format/transcribed_mouse')
    # TODO: Implement Chunking of oligoarray runs
    chunksize = 5
    chunks = []
    for n, gene in enumerate(IntronRetriever()):
        name = gene[0].id.lower()
        # if name not in zak_genes:
        #     continue
        chunks.append(gene)
        if len(chunks) == chunksize:
            with NamedTemporaryFile("w") as fasta_input:
                for gene in chunks:
                    SeqIO.write(gene, fasta_input,
                                'fasta')  # write all sequences to file
                fasta_input.flush()
                res = o.run(fasta_input.name,
                            max_dist=1000,
                            min_length=35,
                            max_length=35,
                            max_tm=100,
                            cross_hyb_temp=72,
                            min_tm=74,
                            secondary_struct_temp=76,
                            max_oligos=1000,
                            timeout=10)
            # TODO: CHECK THAT GROUPBY is working
            for name, probes in res:
                probes['Percent GC'] = 100 * probes["Probe (5'-> 3')"].map(
                    gc_count)
                print(name, len(probes))
                probe_db.insert_many(probes.T.to_dict().values())
            chunks = []
