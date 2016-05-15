from __future__ import print_function, division
import csv
from collections import Counter
from itertools import groupby
import zipfile
from tempfile import NamedTemporaryFile


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement[c] for c in seq.upper())[::-1]


files = []

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


primers = [(f_primer, reverse_complement(r_primer))
           for f_primer, r_primer in primers]
used_primers = []

pgk1_bridges = ['ATTGAGCCAGCAGAAAATGG', 'ATTGAGGGTGTGCTTCGCAC']

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
path = 'temp/readout_oligos.csv'
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

f_outn = "temp/barcode_c_primers.csv"
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
control_oligos += [('Control-Exclude', control_seq)]

# Make Oligo Order
all_probes = oligos + oligos
chip_size = 92918
remaining_oligos = chip_size - len(all_probes) - len(control_oligos)

#Randomly sample remaining oligos to fill chip
from random import sample
all_probes += sample(oligos, remaining_oligos)

all_probes += control_oligos

assert (len(all_probes) <= 92918)

f_name = 'temp/order_9-5-15.csv'
with open(f_name, 'w') as f_out:
    cv = csv.writer(f_out)
    cv.writerow(['Name', 'Oligo'])
    for line in all_probes:
        cv.writerow(line)

# Make an archive storing all files
for file in ['temp/9-6-15.zip', 'temp/9-6-15_ordering.zip']:
    with zipfile.ZipFile(file, 'w',
                         compression=zipfile.ZIP_DEFLATED) as archive:
        if 'ordering' in file:
            with NamedTemporaryFile() as f_temp:
                cv = csv.writer(f_temp)
                t_chip = [(''.join(seq.split(' ')), namel)
                          for seq, namel in all_probes]
                cv.writerows(t_chip)
                f_temp.flush()
                archive.write(f_temp.name, arcname='oligos' + '.csv')
        else:
            archive.write('temp/order_9-5-15.csv', arcname='oligos.csv')
        for file in files:
            archive.write(file, arcname=file.split('/')[1])
