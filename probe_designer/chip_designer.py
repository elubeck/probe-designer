__author__ = 'Eric'
from utils.misc import reverse_complement
from pathlib import Path
import zipfile
from tempfile import TemporaryDirectory, NamedTemporaryFile
import csv
from itertools import groupby
import datetime

pgk1_bridges = ['ATTGAGCCAGCAGAAAATGG', 'ATTGAGGGTGTGCTTCGCAC']
split_pgk1 = [
    ['gcttttgaagcgtgcagaat',
     'aacaccgtgaggtcgaaagg',
     'tcagcttgttggaaagcgac',
     'taggaacgttgaagtccacc',
     'gtttgttatctggttgttct',
     'ttggaacagcagccttgatc',
     'ccattgtccaagcagaattt',
     'ggctcataaggacaacggac',
     'ctggctctaaggagtacttg',
     'gcagagatttgagttcagca',
     'cccacacaatccttcaagaa',
     'acaggcattctcgacttctg',
     'catgaaagcggaggttttcc',
     'ctcagctttaaccttgtttc',
     'tgaggctcggaaagcatcaa',
     'agacatctcctagtttggac',
     'gtcccaaaagcatcattgac',
     'ccaaagccttggcaaagtag',
     'caagatagccaggaagggtc',
     'attgatcagctggatcttgt',
     'ccttaaggaaggtaaaggcc',
     'ccaatctccatgttgttgag',
     'ctccttcttcatcatacaga',
     'tttctcagctttggacatga', ],
    ['gcaaggtaatcttcacacca', 'tgtcagcagtgacaaagtca', 'tataccagaggccacagtag',
     'atatttcttgctgctctcag', 'ccattccaaacaatctgctt', 'ccattcaaataccccaacag',
     'atccatgagtgacttggttc', 'ctagaagtggctttcaccac', 'caccacctatgatagtgatg',
     'ctcacatggctgactttatc', 'aaggactttaccttccagga', 'ggaaccaaaggcaggaaaga',
     'actaggttgacttaggagca', 'gcgctaacaccaaatggaga', 'atctcagccactagctgaat',
     'aactgtttaagggttcctgg', 'aagacgagctgagatgctgt', 'cacagacaaatcctgatgca',
     'gcacaatggttttagtcact', 'ataaatagacgccctctaca', 'tttacagctcacttcctttc',
     'aagctaaccagaggctacat', 'atttcatcagattgccatgc', 'atttctgcagacttacagct']
]

class OligoChip(object):
    f_cutting = {"ScaI": "AGTACT",}
    r_cutting = {"EcoRI": "GAATTC"}

    f_cut_spacer = 'TAG'
    r_cut_spacer = 'GAT'

    t7_sequence = "TAATACGACTCACTATAGGG"

    control_seq = 'GACGCACATATGCGGGCAAGATGCATGCATGCATGCATGCGCTTGCAAGCTTGCAAGCTTGCAAGCTTGCAAGCGGCATCTTCGTGACTGCGGA'


    primers = [
        ('AAGCGCCACGAGTTGTCACG', 'GCATCACTTTGGGCTCGGCT'),
        ('GCCCCATCATGTGCCTTTCC', 'GCCTGCGAGAGGGAAATCCA'),
        ('TGTGCGCTCCGATTGTCCTC', 'GGCCAACAGACCCCATTTGC'),
        ('ATTGAGGGTCTTCGCGTGCC', 'GGTTGCAAAGCGCCGGTTAC'),
        ('AATTGAGCAGCTCGGGCCAC', 'AGTTGCAGGCTTCCATCGCC'),
        ('TGAAGGTGCTGGATGAGCCG', 'GAAAGGACAGCTGCCCGGAA'),
        ('GATGGACTGGCAATTCCGGG', 'GAATATGCCCGGCGAGAACG'),
        ('ACAGTGGAACGCGCCATGAG', 'CGCTTGTTGGCTGGTATGCG'),
        ('CCGCACGCCGTCCTTAAATC', 'AGATCCGGCAGCACGGAAAG'),
        ('GACGCACATATGCGGGCAAG', 'TCCGCAGTCACGAAGATGCC'),
        ('TGCAGCTCCGCGAAATGAAG', 'AATGGCACAGACAGGCAGCG'),
        ('ATCGCTTCGGAGCCTTGGTG', 'CCCAAGGACAAGCAAACCGG'),
        ('ATGCGCTGCAACTGAGACCG', 'CTCGACCAAGGCTGGCACAA'),
    ]



    used_primers = {}
    cutting_oligos = {}
    seqs = {}
    primer_ind = 0

    def add_oligos(self, oligo_list, name, primer_n=None, f_cut_enzyme='',
                   r_cut_enzyme=''):
        if primer_n is None:
            primer_pair = self.primers[self.primer_ind]
            self.primer_ind += 1
        else:
            primer_pair = self.primers[primer_n]
        primer_pair = (primer_pair[0], reverse_complement(primer_pair[1]))
        self.used_primers[name] = primer_pair

        f_cut = ""
        r_cut = ""
        f_spacer = ""
        r_spacer = ""
        f_cut_primer = ""
        r_cut_primer = ""
        if f_cut_enzyme:
            f_cut = self.f_cutting[f_cut_enzyme]
            f_spacer = self.f_cut_spacer
            f_cut_primer = reverse_complement(primer_pair[0] + f_cut + f_spacer)
        if r_cut_enzyme:
            r_cut = self.r_cutting[r_cut_enzyme]
            r_spacer = self.r_cut_spacer
            r_cut_primer = reverse_complement(r_spacer + r_cut + primer_pair[1])
        if any((f_cut_primer, r_cut_primer)):
            self.cutting_oligos[name] = (f_cut_primer, r_cut_primer)
        seqs = []
        for name, oligo in oligo_list:
            vars = {'f_primer': primer_pair[0],
                    'f_cutting': f_cut,
                    'f_spacer': f_spacer,
                    'seq': oligo,
                    'r_spacer': r_spacer,
                    'r_cutting': r_cut,
                    'r_primer': primer_pair[1]}
            seq = "{f_primer} {f_cutting} {f_spacer} {seq} {r_spacer} {r_cutting} {r_primer}".format(**vars)
            seqs.append((name, seq))
        self.seqs[name] = seqs

    def add_pgk1(self, oligos):
        seq_template = "{f_primer} {bridge} TATA {probe} {r_primer}"

        split_pgk1_template = []
        for cn, color_block in enumerate(split_pgk1):
            pgk1_bridge = pgk1_bridges[cn]
            for n, probe in enumerate(color_block):
                split_pgk1_template.append(({
                                        'bridge': reverse_complement(pgk1_bridge),
                                        'probe': probe,
                                    }, 'Split PGK1-Color{}-#{}'.format(cn, n)))

        for name, primer_pair in self.used_primers.items():
            f_primer, r_seq = primer_pair[0], reverse_complement(primer_pair[1])
            primer_dict = {'f_primer': f_primer, 'r_primer': r_seq}
            for template, name in split_pgk1_template:
                template.update(primer_dict)
                oligos.append((name, seq_template.format(**template)))
        return oligos


    def output_chip(self, folder_name="", zip=False):
        # Make an archive storing all files
        if not folder_name:
            folder_name = str(datetime.datetime.now().date())
        out_dir = Path("output").joinpath(folder_name)
        try:
            out_dir.mkdir()
        except FileExistsError:
            pass
        primer_file = out_dir.joinpath("primers.csv")

        def sort_key(i):
            return i[0]

        # Write out primers
        with primer_file.open("w", encoding='utf-8') as f_out:
            p_file = csv.writer(f_out)
            primers = sorted(list(self.used_primers.items()) + list(self.cutting_oligos.items()), key=sort_key)
            p_file.writerow(['F', 'R', 'T7', 'F Cutting', 'R Cutting', ])
            for name, oligos in groupby(primers, key=sort_key):
                oligos = [oligo for nv, oligo_pair in oligos for oligo in oligo_pair]
                t7_primer = (self.t7_sequence + oligos[1]) # T7 + Reverse Primer
                oligo_seq = oligos[:2] + [t7_primer] + oligos[2:]
                p_file.writerow(oligo_seq)

        # Write out unformatted sequences
        for name, s_seqs in self.seqs.items():
            f_name = out_dir.joinpath(("{}.csv".format(name)))
            with f_name.open("w", encoding='utf-8') as f_out:
                cv = csv.writer(f_out)
                cv.writerow(["name", "sequence"])
                cv.writerows(s_seqs)

        # Make order
        order_file = out_dir.joinpath("order.csv")
        with order_file.open("w", encoding='utf-8') as f_out:
            cv = csv.writer(f_out)
            for name, s_seqs in self.seqs.items():
                for n, o in self.add_pgk1(s_seqs):
                    cv.writerow([n, "".join(o.split())])
            cv.writerow(['Control', self.control_seq])



    def __init__(self, primers=None):
        if primers:
            self.primers = primers

if __name__ == "__main__":
    chip = OligoChip()
    chip.add_oligos([("DOG", "GEFA"), ("CAT", "DEAAA")], name='BOBDOLE', f_cut_enzyme="ScaI", r_cut_enzyme="EcoRI")
    chip.output_chip()