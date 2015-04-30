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


class Introns(object):
    """
    Iterator that iterates over introns retrieved from UCSC.
    """

    def get_ids(self):
        # Table of every RefSeqID -> GeneName
        with open("db/Refseq2Gene.tsv", "r") as tsvin:
            return dict(csv.reader(tsvin, delimiter='\t'))

    def get_introns(self, repeat_masking=True, count_introns=True):
        # Refseq intron fasta indexed by ID downloaded from UCSC genome browser
        # Repetive elements are masked in lowercase
        # By default repeats are masked
        used_names = []
        with gzip.open("db/refseq_intron_maskedLower.gzip.gz", 'rb') as f:
            for fasta in SeqIO.parse(f, 'fasta'):
                id = fasta.id.split("refGene_")[1]
                ref_name = self.refseq_id[id]
                if count_introns:
                    used_names.append(ref_name)
                    intron_num = used_names.count(ref_name)
                    fasta.name = "{}_Intron#{}".format(ref_name, intron_num)
                else:
                    fasta.name = ref_name
                if repeat_masking:
                    fasta.seq = Seq("".join([c for c in fasta.seq
                                             if c.isupper()]))
                yield fasta

    def __iter__(self):
        for intron in self.get_introns(repeat_masking=self.repeat_masking):
            yield intron

    def __init__(self, repeat_masking=True, count_introns=True):
        self.refseq_id = self.get_ids()
        self.repeat_masking = repeat_masking
        self.count_introns = count_introns
        self.organism = "mouse"


class Intron(object):
    def to_fasta(self, chunk_size=1000):
        rev_comp = self.record.seq.reverse_complement()
        for chunk_n, chunk_start in enumerate(range(0, len(rev_comp), 1000)):
            seq_chunk = rev_comp[chunk_start:chunk_start + 1000]
            fasta_str = ">{}=Chunk:{}\n{}\n".format(self.name, chunk_n,
                                                    seq_chunk)
            yield fasta_str

    def __init__(self, record):
        self.record = record
        # TODO: Embed Intron number in name
        self.name = "{}_{}".format(self.record.name, "intron")


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

if __name__ == "__main__":
    # Shortcut code to count total # of introns
    # tot_introns = len([i for i in Introns(repeat_masking=False, count_introns=False)])
    tot_introns = 30759
    introns = Introns()
    intron_db = dataset.connect("sqlite:///db/intron_probes.db")
    probe_db = intron_db['mouse']
    used_probes = set([row['Name'] for row in probe_db.distinct("Name")])
    queue = []
    # Temporarily changed for zak genes
    undesigned_introns = ((n, i) for n, i in enumerate(introns)
                          if i.name not in used_probes)
    # if i.name.lower().split("_intron")[0] in zak_genes)
    intron_1 = ((n, i) for n, i in undesigned_introns
                if i.name.endswith("Intron#1"))
    o = OligoarrayDesigner(
        blast_db='/home/eric/blastdb/old_format/mouse_introns_revcomp')
    for design_num, (n, i) in enumerate(intron_1):
        print(design_num, n, i.name)
        f_loc = "temp/temp.fasta"
        n_chunks = len(i.seq) // 1000
        # If too many chunks exist pick random ones
        tot_chunks = 150
        if n_chunks > tot_chunks:
            sub_set = random.sample(range(0, len(i.seq), 1000), tot_chunks)
        used = []
        with open(f_loc, "w") as f:
            rev_comp = i.seq.reverse_complement()
            for chunk_n, chunk_start in enumerate(range(0, len(rev_comp),
                                                        1000)):
                # Skip chunks not in list
                if n_chunks > tot_chunks and chunk_start not in sub_set:
                    continue
                seq_chunk = rev_comp[chunk_start:chunk_start + 1000]
                fasta_str = ">{}=Chunk:{}\n{}\n".format(i.name, chunk_n,
                                                        seq_chunk)
                f.write(fasta_str)
                used.append(chunk_n)
        n_chunks = len(used)
        n_probes = 0
        n_iterations = 0
        # Set oligo number threshold to minimum to hit n_chunks
        max_oligos = int(ceil(tot_chunks / n_chunks))
        # Reduce Computation Time: Only design small set of oligos for every chunk
        # Iterate this until at least 100 probes designed
        p_num = 9553553535  # Just a random number for initiatlization
        while n_probes < 100:
            print(max_oligos, n_iterations, n_chunks * max_oligos)
            res = o.run(f_loc,
                        max_dist=1000,
                        min_length=35,
                        max_length=35,
                        max_tm=100,
                        cross_hyb_temp=72,
                        min_tm=74,
                        secondary_struct_temp=76,
                        max_oligos=max_oligos,
                        timeout=10)
            if len(res) > 1:
                raise Exception("This shouldn't happen")
            try:
                name, probes = res[0]
            except IndexError as e:
                print("No Probes returned")
                if n_iterations != 0:
                    print("Giving up making probes for {}".format(i.name))
                    break
                n_iterations += 1
                max_oligos *= 2
                continue
            n_probes = len(probes)
            if n_probes == p_num and n_iterations != 0:
                print("Can't make more probes for {}".format(i.name))
            elif n_probes < 100:
                if n_iterations > 10:
                    raise Exception("Too many iterations")
                n_iterations += 1
                max_oligos *= 2
                p_num = n_probes
                print(
                    "Trying to design more probes.  Last cycle designed {}".format(
                        len(probes)))
                continue
            progress = 100 * (n / tot_introns)
            print("{:.02f}%: {} with {} probes after {} iterations".format(
                progress, name, len(probes), n_iterations))
            probes['Percent GC'] = 100 * probes["Probe (5'-> 3')"].map(gc_count)
            probe_db.insert_many(probes.T.to_dict().values())
            break
    print(n)
