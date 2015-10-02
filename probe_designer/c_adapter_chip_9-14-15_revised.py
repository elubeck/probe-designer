from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
# from builtins import *
import dataset
import intron_designer2
from collections import defaultdict
from itertools import groupby
import blaster2
from progressbar import ProgressBar
import mRNA_designer
import blaster2
import csv
import json
from get_seq import reverse_complement
from collections import Counter
import zipfile
from tempfile import NamedTemporaryFile


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement[c] for c in seq.upper())[::-1]


files = []

seq_template = "{f_primer} {f_cutting} TAG {bridge} TATA {probe} GAT {r_cutting} {r_primer}"
f_cutting = 'AGTACT'  # ScaI
r_cutting = 'GAATTC'  # EcoRI

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
    # ('GACGCACATATGCGGGCAAG', 'TCCGCAGTCACGAAGATGCC'), # Commented out due to presence in control
    ('TGCAGCTCCGCGAAATGAAG', 'AATGGCACAGACAGGCAGCG'),
    ('ATCGCTTCGGAGCCTTGGTG', 'CCCAAGGACAAGCAAACCGG'),
    ('GACGCACATATGCGGGCAAG', 'GGCATCTTCGTGACTGCGGA'),
    ('ATGCGCTGCAACTGAGACCG', 'CTCGACCAAGGCTGGCACAA')
]  # yapf: disable

pgk1_bridges = ['ATTGAGCCAGCAGAAAATGG', 'ATTGAGGGTGTGCTTCGCAC']

primers = [(f_primer, reverse_complement(r_primer))
           for f_primer, r_primer in primers]

prime_iter = iter(primers)

initiators = ["gCATTCTTTCTTgAggAgggCAgCAAACgggAAgAg",
              "AgCTCAgTCCATCCTCgTAAATCCTCATCAATCATC",
              "AAAgTCTAATCCgTCCCTgCCTCTATATCTCCACTC",
              "CACATTTACAgACCTCAACCTACCTCCAACTCTCAC",
              "CACTTCATATCACTCACTCCCAATCTCTATCTACCC"]

all_probes = []

used_primers = []

with open('temp/validated_bridges-1k_old.txt', 'r') as f_in:
    bridges = [line.strip('\n') for line in f_in]
iter_bridges_old = iter(bridges)

from collections import defaultdict
brain_genes = defaultdict(list)
with open("temp/brain_mrna_9-14-15.csv", "rb") as f_in:
    c = csv.reader(f_in)
    for line in c:
        gene, n, probe = line
        brain_genes[gene].append(probe)

barcodes = []
with open("db/barcode_5hyb.csv", "r") as f_in:
    c = csv.reader(f_in)
    barcodes = [[int(element) - 1 for element in line] for line in c]

from random import shuffle
shuffle(barcodes)
brain_barcodes = dict(zip(brain_genes.keys(), barcodes))

brain_primers = [prime_iter.next() for n in range(5)]
brain_probes = []
for gene, code in brain_barcodes.items():
    for n, digit in enumerate(code):
        init = initiators[digit]
        f_primer, r_primer = brain_primers[n]
        for pn, probe in enumerate(brain_genes[gene]):
            vals = {
                'f_primer': f_primer,
                'f_cutting': f_cutting,
                'bridge': init,
                'probe': probe,
                'r_cutting': r_cutting,
                'r_primer': r_primer
            }
            brain_probes.append((
                "{}-{}-{}".format(gene, n + 1, pn), seq_template.format(**vals)
            ))
        used_primers.append(
            ("Brain hyb {}".format(n + 1), (f_primer, r_primer)))

all_probes += brain_probes

t_dir = "9-14-15"
import os
# os.mkdir('temp/{}'.format(t_dir))

# Write barcodes to file
with open("temp/{}/barcodes_8-27-15_brain.csv".format(t_dir), "wb") as f_out:
    cv = csv.writer(f_out)
    cv.writerow(['name', 'Hyb1', 'Hyb2', 'Hyb3', 'Hyb4', 'Hyb5', ])
    for name, barcode in brain_barcodes.items():
        cv.writerow([name] + list(barcode))

#### Get existing 10K Oligos
introns_10k = []
with open('temp/barcodes_10k_introns_6-26-15.csv', 'r') as f_in:
    for n, line in enumerate(csv.reader(f_in)):
        if n == 0: continue
        name = line[0]
        # These are actual bridge sequences, not reverse complement
        bridge = line[1]
        barcode = [int(i) for i in line[2:]]
        introns_10k.append((name, bridge, barcode))

### Get passed bridges
with open("db/passed_bridges_mouse.txt", 'r') as f_in:
    bridges = {line.strip('\n') for line in f_in}

used_bridges = {
    b.lower()
    for i in introns_10k for b in (i[1], reverse_complement(i[1]))
}

available_bridges = bridges.difference(used_bridges)

# Check that PGK1 bridges not used
assert (len([True for bridge in pgk1_bridges
             if bridge.upper() in available_bridges]) == 0)

# Check that nothing binds itself
bridge_counter = Counter(
    list(available_bridges) + [reverse_complement(bridge).lower()
                               for bridge in available_bridges])
too_many = [(k, v) for k, v in bridge_counter.iteritems() if v != 1]
assert (len(too_many) == 0)

# Randomly scramble bridges
from random import shuffle
bridges = list(available_bridges)
shuffle(bridges)

# Make bridge iterator so no bridge is used twice
iter_bridges = iter(available_bridges)

##### Get Tertiary Oligos.  These will be dye labelled
n_hybridizations = 8
n_colors = 4
readout_oligos = [[iter_bridges.next() for col in range(n_colors)]
                  for hyb in range(n_hybridizations)]

# Make bridge file
path = 'temp/{}/readout_oligos.csv'.format(t_dir)
with open(path, 'w') as f_out:
    cv = csv.writer(f_out)
    cv.writerow(["Hybridization"] + ["Color {}".format(n + 1)
                                     for n in range(n_colors)])
    for n, line in enumerate(readout_oligos):
        cv.writerow([n + 1] + line)

files.append(path)

secondary_template = "{f_primer_1st} {b1} TATA {bridge} ATAT {b2} {r_primer_1st}"
tertiary_template_1 = "{f_primer_2nd} {b1_rc} {h1} {h2} {h3} {h4} {r_primer_2nd}"
tertiary_template_2 = "{f_primer_2nd} {h5} {h6} {h7} {h8} {b2_rc} {r_primer_2nd}"

f_outn = "temp/{}/barcode_c_primers.csv".format(t_dir)
with open(f_outn, 'w') as f_out:
    cv = csv.writer(f_out)
    cv.writerow(['Name', 'Bridge', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7',
                 'B8'])
    oligos = []
    for gene in introns_10k:
        name = gene[0]
        bridge = gene[1]
        barcode = gene[2]
        code_seqs = [reverse_complement(readout_oligos[n_hyb][digit - 1])
                     for n_hyb, digit in enumerate(barcode)]
        # Get bridge sequences
        b = [iter_bridges.next() for i in range(2)]
        b_rev = [reverse_complement(i) for i in b]

        # Get Primers
        f_primer_1st, r_primer_1st = primers[0]
        used_primers.append(("Secondary", (f_primer_1st, r_primer_1st)))

        # Get Primers
        f_primer_2nd, r_primer_2nd = primers[1]
        used_primers.append(("Tertiary", (f_primer_2nd, r_primer_2nd)))

        vars = {
            'f_primer_1st': f_primer_1st,
            'r_primer_1st': r_primer_1st,
            'f_primer_2nd': f_primer_2nd,
            'r_primer_2nd': r_primer_2nd,
            'bridge': bridge,
            'b1': b[0],
            'b2': b[1],
            'b1_rc': reverse_complement(b[0]),
            'b2_rc': reverse_complement(b[1]),
            'h1': code_seqs[0],
            'h2': code_seqs[1],
            'h3': code_seqs[2],
            'h4': code_seqs[3],
            'h5': code_seqs[4],
            'h6': code_seqs[5],
            'h7': code_seqs[6],
            'h8': code_seqs[7],
        }
        oligos.append(
            ("{}-Secondary".format(name), secondary_template.format(**vars)))
        oligos.append(
            ("{}-Tertiary-1".format(name), tertiary_template_1.format(**vars)))
        oligos.append(
            ("{}-Tertiary-2".format(name), tertiary_template_2.format(**vars)))
        cv.writerow([name, vars['bridge'], vars['h1'], vars['h2'], vars['h3'],
                     vars['h4'], vars['h5'], vars['h6'], vars['h7'],
                     vars['h8']])
files.append(f_outn)

##### Setup split PGK1
split_pgk1 = [
    ['gcttttgaagcgtgcagaat', 'aacaccgtgaggtcgaaagg', 'tcagcttgttggaaagcgac',
     'taggaacgttgaagtccacc', 'gtttgttatctggttgttct', 'ttggaacagcagccttgatc',
     'ccattgtccaagcagaattt', 'ggctcataaggacaacggac', 'ctggctctaaggagtacttg',
     'gcagagatttgagttcagca', 'cccacacaatccttcaagaa', 'acaggcattctcgacttctg',
     'catgaaagcggaggttttcc', 'ctcagctttaaccttgtttc', 'tgaggctcggaaagcatcaa',
     'agacatctcctagtttggac', 'gtcccaaaagcatcattgac', 'ccaaagccttggcaaagtag',
     'caagatagccaggaagggtc', 'attgatcagctggatcttgt', 'ccttaaggaaggtaaaggcc',
     'ccaatctccatgttgttgag', 'ctccttcttcatcatacaga', 'tttctcagctttggacatga', ],
    ['gcaaggtaatcttcacacca', 'tgtcagcagtgacaaagtca', 'tataccagaggccacagtag',
     'atatttcttgctgctctcag', 'ccattccaaacaatctgctt', 'ccattcaaataccccaacag',
     'atccatgagtgacttggttc', 'ctagaagtggctttcaccac', 'caccacctatgatagtgatg',
     'ctcacatggctgactttatc', 'aaggactttaccttccagga', 'ggaaccaaaggcaggaaaga',
     'actaggttgacttaggagca', 'gcgctaacaccaaatggaga', 'atctcagccactagctgaat',
     'aactgtttaagggttcctgg', 'aagacgagctgagatgctgt', 'cacagacaaatcctgatgca',
     'gcacaatggttttagtcact', 'ataaatagacgccctctaca', 'tttacagctcacttcctttc',
     'aagctaaccagaggctacat', 'atttcatcagattgccatgc', 'atttctgcagacttacagct']
]  # yapf: disable

control_oligos = []
# Make split PGK1 probes except for template
split_pgk1_template = []
for cn, color_block in enumerate(split_pgk1):
    pgk1_bridge = pgk1_bridges[cn]
    for n, probe in enumerate(color_block):
        split_pgk1_template.append(({
            'bridge': reverse_complement(pgk1_bridge),
            'probe': probe,
        }, 'Split PGK1-Color{}-#{}'.format(cn, n)))

seq_template = "{f_primer} {bridge} TATA {probe} {r_primer}"

for name, (f_primer, r_primer) in set(used_primers):
    primer_dict = {'f_primer': f_primer, 'r_primer': r_primer}
    for template, name2 in split_pgk1_template:
        template.update(primer_dict)
        control_oligos.append(
            (name2 + '-{}'.format(name), seq_template.format(**template)))

control_seq = seq_template.format(**{
    'f_primer': 'GACGCACATATGCGGGCAAG',
    'bridge': 'ATGCATGCATGCATGCATGC',
    'probe': 'GCTTGCAAGCTTGCAAGCTTGCAAGCTTGCAAGCA',
    'r_primer': reverse_complement('TCCGCAGTCACGAAGATGCC')
})

used_primers.append(
    (('GACGCACATATGCGGGCAAG', 'GGCATCTTCGTGACTGCGGA'), "Control"))
control_oligos += [(control_seq, 'Control-Exclude')]

f_name = "temp/{}/used_primers.csv".format(t_dir)
with open(f_name, 'wb') as f_out:
    bob = csv.writer(f_out)
    bob.writerow(['Name', 'Forward', 'Reverse', ])
    pr = list(set(used_primers))
    up = [p for p in pr if 'Control' not in p[1]]
    for name, (f, r) in up:
        bob.writerow([name, f, r])
    control = used_primers[-1]
    bob.writerow([control[1], control[0][0], control[0][1]])
files.append(f_name)

# Make Oligo Order
all_probes = all_probes + oligos
chip_size = 92918
remaining_oligos = chip_size - len(all_probes) - len(control_oligos)

#Randomly sample remaining oligos to fill chip
from random import sample
all_probes += sample(oligos, remaining_oligos)

all_probes += control_oligos

assert (len(all_probes) <= 92918)

f_name = 'temp/{}/order_9-5-15.csv'.format(t_dir)
with open(f_name, 'w') as f_out:
    cv = csv.writer(f_out)
    cv.writerow(['Name', 'Oligo'])
    for line in all_probes:
        cv.writerow(line)
files.append(f_name)

# Make an archive storing all files
for file in ['temp/{}/9-6-15.zip'.format(t_dir),
             'temp/{}/9-6-15_ordering.zip'.format(t_dir)]:
    with zipfile.ZipFile(file, 'w',
                         compression=zipfile.ZIP_DEFLATED) as archive:
        if 'ordering' in file:
            with NamedTemporaryFile() as f_temp:
                cv = csv.writer(f_temp)
                t_chip = [(''.join(seq.split(' ')), ''.join(seq.split(' ')))
                          for seq, namel in all_probes]
                cv.writerows(t_chip)
                f_temp.flush()
                archive.write(f_temp.name, arcname='oligos' + '.csv')
        else:
            archive.write('temp/{}/order_9-5-15.csv'.format(t_dir),
                          arcname='oligos.csv')
        for file in files:
            archive.write(file, arcname=file.split('/')[2])
