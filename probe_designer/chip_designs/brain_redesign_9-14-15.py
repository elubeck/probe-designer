from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from builtins import *
import dataset
import intron_designer2
from collections import defaultdict
from itertools import groupby
import blaster2
from progressbar import ProgressBar
import probe_designer.mRNA_designer
import blaster2
import csv
import json
from probe_designer.utils.misc import reverse_complement
from collections import Counter


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement[c] for c in seq.upper())[::-1]


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
with open("temp/brain_8-24-15.csv", "r") as f_out:
    c = csv.reader(f_out)
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
            brain_probes.append((seq_template.format(**vals),
                                 "{}-{}-{}".format(gene, n + 1, pn)))
        used_primers.append(
            ((f_primer, r_primer), "Brain hyb {}".format(n + 1)))

all_probes += brain_probes

# Write barcodes to file
with open("temp/barcodes_8-27-15_brain.csv", "wb") as f_out:
    cv = csv.writer(f_out)
    cv.writerow(['name', 'Hyb1', 'Hyb2', 'Hyb3', 'Hyb4', 'Hyb5', ])
    for name, barcode in brain_barcodes.items():
        cv.writerow([name] + list(barcode))

zak_genes = defaultdict(list)
with open("temp/zak_8-24-15.csv", "r") as f_out:
    c = csv.reader(f_out)
    for line in c:
        gene, n, probe = line
        zak_genes[gene].append(probe)

barcodes = []
with open("db/barcode_4hyb.csv", "r") as f_in:
    c = csv.reader(f_in)
    barcodes = [[int(element) - 1 for element in line] for line in c]

from random import sample, shuffle
zak_picks = {
    'Dazl', 'Dppa2', 'Fgf4', 'Fgf5', 'Kdm4d', 'Klf2', 'Klf4', 'Lin28a',
    'Lin28b', 'Mecp2', 'Mettl14', 'Sall4', 'Socs3', 'Sox2', 'Xist', 'Zfp42',
    'Foxa2', 'Hdac2'
}
necessary_genes = zak_picks.intersection(zak_genes.keys())
remaining_genes = set(zak_genes.keys()).difference(zak_picks)
n_choices = len(barcodes) - len(necessary_genes)
remaining_picks = sample(remaining_genes, n_choices)

# Barcodes of chosen genes
shuffle(barcodes)
zak_barcodes = dict(zip(necessary_genes | set(remaining_picks), barcodes))

zak_primers = [prime_iter.next() for n in range(4)]
zak_probes = []
for gene, code in zak_barcodes.items():
    for n, digit in enumerate(code):
        init = initiators[digit]
        f_primer, r_primer = zak_primers[n]
        for pn, probe in enumerate(zak_genes[gene]):
            vals = {
                'f_primer': f_primer,
                'f_cutting': f_cutting,
                'bridge': init,
                'probe': probe,
                'r_cutting': r_cutting,
                'r_primer': r_primer
            }
            zak_probes.append((seq_template.format(**vals),
                               "{}-{}-{}".format(gene, n + 1, pn)))
            used_primers.append(((f_primer, r_primer), "Zak hyb {}".format(n)))

all_probes += zak_probes

# Write barcodes to file
with open("temp/barcodes_8-27-15_zak.csv", "wb") as f_out:
    cv = csv.writer(f_out)
    cv.writerow(['name', 'Hyb1', 'Hyb2', 'Hyb3', 'Hyb4', ])
    for name, barcode in zak_barcodes.items():
        cv.writerow([name] + list(barcode))

with open('temp/chip_1.csv', 'r') as f_in:
    oligos = [(oligo, name) for oligo, name in csv.reader(f_in)]

from itertools import groupby
primers_to_use = {'AAGCGCCACGAGTTGTCACG', 'ATTGAGGGTCTTCGCGTGCC'}

two_thousand_genes = dict()
two_thousand_bridges = dict()
for f_primer, group in groupby(oligos, lambda x: x[0].split()[0]):
    if f_primer in primers_to_use:
        for gene, oligos in groupby(group, lambda x: x[1].split('-')[0]):
            if 'Split PGK1' not in gene:
                ol = list(oligos)
                bridges = list({b[0].split()[3] for b in ol})
                assert (len(bridges) == 1)
                two_thousand_genes[gene] = ol
                two_thousand_bridges[gene] = reverse_complement(bridges[0])

pgk1_bridges = ['ATTGAGCCAGCAGAAAATGG', 'ATTGAGGGTGTGCTTCGCAC']
u_bridges = set(list(two_thousand_bridges.values()) + pgk1_bridges)

iter_bridges = (b for b in iter_bridges_old if b not in u_bridges)

n_hybridizations = 8
n_colors = 5
tert_seqs = [[iter_bridges.next() for col in range(n_colors)]
             for hyb in range(n_hybridizations)]

with open('temp/quats_8hybs.csv', 'wb') as f_out:
    f_out2 = csv.writer(f_out)
    f_out2.writerow(['Hyb', 'Color 1', 'Color 2', 'Color 3', 'Color 4',
                     'Color 5'])
    for n, hyb in enumerate(tert_seqs):
        f_out2.writerow([n + 1] + list(map(reverse_complement, hyb)))

n_genes = 2000
n_secondary = 2
sec_seqs = [[iter_bridges.next() for c in range(n_secondary)]
            for gene in range(n_genes)]
iter_secs = iter(sec_seqs)

with open('temp/tertiaries_2k_genes.csv', 'wb') as f_out:
    f_out2 = csv.writer(f_out)
    f_out2.writerow(['2ndary', '3rd_1', '3rd_2'])
    for gene, bridges in sec_seqs:
        f_out2.writerow([gene] + list(map(reverse_complement, bridges)))

with open('db/barcode10k.csv', 'r') as f_in:
    replace = {1: 4, 2: 2, 3: 3, 4: 1}  # Swap 1 and 4
    barcodes_10k = list({
        tuple(map(lambda x: replace[int(x)], line))
        for line in csv.reader(f_in)
    })

sorted(barcodes_10k)

secondary = dict()
tertiary = dict()

secondary_template = "{f_primer_2nd} {f_cutting} TAG {b1} {bridge} {b2} GAT {r_cutting} {r_primer_2nd}"
tertiary_template_1 = "{f_primer_3rd} {b1_rc} {h1} {h2} {h3} {h4} {r_primer_3rd}"
tertiary_template_2 = "{f_primer_3rd} {h5} {h6} {h7} {h8} {b2_rc} {r_primer_3rd}"

f_primer_2nd, r_primer_2nd = prime_iter.next()
f_primer_3rd, r_primer_3rd = prime_iter.next()

used_primers.append(((f_primer_2nd, r_primer_2nd), "C Bridge Secondary"))
used_primers.append(((f_primer_3rd, r_primer_3rd), "C Bridge Tertiary"))

branch_secondary = []
branch_tertiaries = []

bc_2k = {}
for n2, (gene, bridge) in enumerate(two_thousand_bridges.iteritems()):
    b1, b2 = iter_secs.next()

    # Get the Barcodes
    gene_code = barcodes_10k[n2]
    code_seqs = [tert_seqs[n_hyb][digit - 1]
                 for n_hyb, digit in enumerate(gene_code)]
    bc_2k[gene] = (reverse_complement(bridge), code_seqs)
    split = len(code_seqs) // 2

    vars = {
        'f_primer_2nd': f_primer_2nd,
        'r_primer_2nd': r_primer_2nd,
        'f_primer_3rd': f_primer_3rd,
        'r_primer_3rd': r_primer_3rd,
        'f_cutting': f_cutting,
        'r_cutting': r_cutting,
        'bridge': bridge,
        'b1': b1,
        'b1_rc': reverse_complement(b1),
        'b2': reverse_complement(b2),
        'b2_rc': b2,
        'h1': code_seqs[0],
        'h2': code_seqs[1],
        'h3': code_seqs[2],
        'h4': code_seqs[3],
        'h5': code_seqs[4],
        'h6': code_seqs[5],
        'h7': code_seqs[6],
        'h8': code_seqs[7],
    }
    branch_secondary.append((secondary_template.format(**vars),
                             "{}-BranchedSecondary".format(gene)))
    branch_tertiaries.append((tertiary_template_1.format(**vars),
                              "{}-1-BranchedTertiary".format(gene)))
    branch_tertiaries.append((tertiary_template_2.format(**vars),
                              "{}-2-BranchedTertiary".format(gene)))

all_probes += branch_secondary
all_probes += branch_tertiaries

# Write barcodes to file
with open("temp/barcodes_8-27-15_2kgenes.csv", "wb") as f_out:
    cv = csv.writer(f_out)
    cv.writerow(['name', '2ndary_bridge', 'Hyb1', 'Hyb2', 'Hyb3', 'Hyb4',
                 'Hyb5', 'Hyb6', 'Hyb7', 'Hyb8'])
    for name, (bridge, barcode) in bc_2k.items():
        cv.writerow([name, bridge] + list(barcode))

#TODO: REMOVE THIS

two_color_coloc = defaultdict(list)
with open('temp/5-27-15_intron_probes2.csv', 'r') as f_in:
    for probe, gene in csv.reader(f_in, delimiter=b'\t'):
        probe_seq = probe[53:88]
        two_color_coloc[gene].append(probe_seq)

select_probes = sample(two_color_coloc.keys(), 98)
#TODO: REMOVE THIS

# select_probes = []
two_color_coloc = defaultdict(list)
with open('temp/5-27-15_intron_probes2.csv', 'r') as f_in:
    for probe, gene in csv.reader(f_in, delimiter=b'\t'):
        if gene in select_probes:
            probe_seq = probe[53:88]
            two_color_coloc[gene].append(probe_seq)

init_2c = [initiators[0], initiators[2]]

f_primer, r_primer = prime_iter.next()
used_primers.append(((f_primer, r_primer), "Colocalizing Introns"))

from random import shuffle
for gene, probes in two_color_coloc.items():
    shuffle(probes)
    mid = len(probes) // 2
    for n, section in enumerate((probes[:mid], probes[mid:])):
        for pn, probe in enumerate(probes):
            vals = {
                'f_primer': f_primer,
                'f_cutting': f_cutting,
                'bridge': init_2c[n],
                'probe': probe,
                'r_cutting': r_cutting,
                'r_primer': r_primer
            }
            all_probes.append(
                (seq_template.format(**vals),
                 "{}-Color{}-{}-Colocalizing".format(gene, n, pn + 1)))

pcdha = [u'AATTTGTCTGGCAACTCACCGGGACCGGAC', u'AAAGGTCCAGCTGTTGCTGTTGACACCGGC',
         u'GATGGAGATGATTGCAGGAGATCCTGGGAT', u'TTGTCAATTTGGTTGTTAGCAGGCTCCTGC',
         u'GCCCACCGGAGGGGACACCTCTCCTGCCTC', u'TTTCCCCCTCCAGATGAAAGACTCCCGGTG',
         u'GAACTGCTGGCTCGTGTGCATACTGCAGGT', u'GCTGTGATTGGCAGGTGGCCGGGATAATAC',
         u'ATGCAAGAGGCAGCCCCAGAGTGACCTGAA', u'CATCCCAGCCAATACCAGGCTCTCCAACAC',
         u'AACAAGGTGGAGGTGTAGGCTAGCACACCG', u'GGCTCATGTCATAGGAGAAAGGAGCGGGGA',
         u'CCCTTGCCAGCAGGATTCTTGGCAAGAGAC', u'ACTGGTCACTGTTGTCCGTCGTGCTGTTCC',
         u'GAAGCACAGGGCAGAGAGTCCAACTCTCAC', u'ATCAGATCTGCTTCTGTTGCCCTGGCCGAC',
         u'AGCCAGCACACTTGAGGGGCTCCTAAAGCT', u'TTCTCCTTTGCGGCGGGGGAGAAAGGAAAC',
         u'GCCGGTGGGGATGGTTTTGTCTGTTAGAGG', u'GCCACACATACCCAGTGACCACAAACACCC',
         u'CTCAACCTGAATACCTTTCAGCCTGCCCTG', u'GGGCTTCAGTGGGGCAAAATGATCACAACC',
         u'CTTTGAGCCCCAGCCCTCCTTGTGTTAACT', u'GGATTCTCTTTTGAGGGAAGGTTGGGGGGA',
         u'AGCGTGAGCACAAGGACGTGTACAGGAAAC', u'GCCAGAGTGGAGAGAGAGAGAGTTAACACC',
         u'GGCTAAGTGGCTTGTTTCCATTTGGTGGCC', u'TGTGCTTGCATACCCTCAGGGACTTTCTTC',
         u'TGTAGCAATGGTTCACATGCGCTGTCAACG', u'CCTGCAGAGCTGAAGTTACACACCCTTCAA',
         u'GCCCGATTCTCCTAGTTCAAAAGGCACGTA', u'CTGGGTCTTGATATCGTGGGTTAGAGAAGG',
         u'CACAGAACCAGAGGAATGAAAAGTCTGCTC', u'CAAGCACAGAGTAAGTATATCTCAGGAGTC',
         u'CAGCAACGATCATTTCAAATGGCTGTAGAC', u'TTAAATGCGATCAAGACAGACACAGGGAAC',
         u'CAAAGCAGTCTGCAAATTCAAAGAGGTCAG', u'GTTGAGAGATTGAGAAAGTCTGGTCCTTGT',
         u'GTTCCAGGCAATTCTCTTAAAGCACTTCTG', u'GCATATGTGTTGTTGGAAAAGGATGAGACC',
         u'GGGAGATTGATTCAGGTAACATCATAGGTG', u'GGTTAGCAATTCAGAGTGTTCTAGAAAAGG',
         u'GGGAGAATGAAGCTAAAAGTAGACAAGACA', u'ATTCCACTTTTCCACAGGTGTACAAAAGGA',
         u'AACAGCATAAAGTGAGGTCTTGGGTTTTTC', u'GTTGCTACATGGGAGAGATTTGCCATTATT',
         u'ACTTTCTACAAGCTGATATTCCGGATAGCA', u'GCCCTTGGCTGACCATAAATAAAATTCCTT',
         u'CTGGGTCTGCTTTTTAGGTTCTGAAAACAA', u'TTATGCCTTTTAAAACAGACTCCTACCTCT',
         u'TGGGCTTACATTCAACAATAAAGCTTAAAG', u'ACTCTAGGTAATAGTTGCATTAACATTCAC',
         u'TACTGTGTCTTTCATACACAAAAGCTGATT', u'GACTTAGGAAAGCTTTTTAAGACAAAACCT',
         u'GAATCACCAGCATCTTGATTTAAGAAAGAA', u'TTGTGTAATCAAAGTGCACATCATTTTTCT',
         u'TCCTTTAAAGATTTATAGCTATCGACTTTC', u'TATGCATGAGAAAACAAGATTTTTCTTCAG',
         u'ACTAGGAAAAAAAAAAACCAGCTTTATCCA', u'CAAACATCTAAGTTGTGCACAATAAACTTA',
         u'TCCACATGACGTGTTTTATCCTTACAGCTG']

color_1 = pcdha[::3]
color_2 = pcdha[1::2]
color_3 = pcdha[2::2]

f_primer, r_primer = prime_iter.next()
init_3c = [initiators[1], initiators[2], initiators[3]]

pcdh_template = "{f_primer} {f_cutting} TAG {hcr_seq} TATA {probe} GAT {r_cutting} {r_primer}"
for n, (init, color) in enumerate(zip(init_3c, (color_1, color_2, color_3))):
    for pn, probe in enumerate(color):
        vals = {
            'f_primer': f_primer,
            'f_cutting': f_cutting,
            'r_cutting': r_cutting,
            'hcr_seq': init,
            'probe': probe,
            'r_primer': r_primer,
            'filler': 'AGTC' * 5
        }
        used_primers.append(((f_primer, r_primer), "PCDH Colocalizing"))
        all_probes.append((pcdh_template.format(**vals),
                           "PCDHaConst-{}-{}".format(n + 1, pn + 1)))

# Setup split PGK1
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

# Make split PGK1 probes except for template
split_pgk1_template = []
for cn, color_block in enumerate(split_pgk1):
    pgk1_bridge = pgk1_bridges[cn]
    for n, probe in enumerate(color_block):
        split_pgk1_template.append(({
            'f_cutting': f_cutting,
            'bridge': reverse_complement(pgk1_bridge),
            'probe': probe,
            'r_cutting': r_cutting,
        }, 'Split PGK1-Color{}-#{}'.format(cn, n)))

n_used_primers = {
    (oligo.split(' ')[0], oligo.split(' ')[-1])
    for oligo, name in all_probes
}
for f_primer, r_primer in n_used_primers:
    primer_dict = {'f_primer': f_primer, 'r_primer': r_primer}
    for template, name in split_pgk1_template:
        template.update(primer_dict)
        all_probes.append((seq_template.format(**template), name))

control_seq = seq_template.format(**{
    'f_primer': 'GACGCACATATGCGGGCAAG',
    'f_cutting': f_cutting,
    'r_cutting': r_cutting,
    'bridge': 'ATGCATGCATGCATGCATGC',
    'probe': 'GCTTGCAAGCTTGCAAGCTTGCAAGCTTGCAAGCA',
    'r_primer': reverse_complement('TCCGCAGTCACGAAGATGCC')
})

used_primers.append(
    (('GACGCACATATGCGGGCAAG', 'GGCATCTTCGTGACTGCGGA'), "Control"))
all_probes += [(control_seq, 'Control-Exclude')]


def primer_getter(seq):
    vals = seq[0].split(' ')
    return (vals[0], vals[-1])


discovered_primers = []
for primers, probe in groupby(sorted(all_probes), primer_getter):
    discovered_primers.append(primers)
discovered_primers = set(discovered_primers)
existing_primers = {p for p, n in set(used_primers)}

assert (existing_primers == discovered_primers)


def sort_key(i):
    # return i.strip(range(10))
    return i[1].split('hyb')[0].strip()

    # return " ".join(i[1].split()[:2])


cut_groups = ['Brain', 'Zak', 'C Bridge Secondary', 'Colocalizing Introns'
              'PCDH Colocalizing']
primer_lookup = set(used_primers)
t7_seq = "TAATACGACTCACTATAGGG"
with open('temp/primers_8-27-15.csv', 'wb') as f_out:
    cv = csv.writer(f_out)
    cv.writerow(['Group', 'Name', 'Forward', 'Reverse', 'T7', 'Forward Cutting'
                 'Reverse Cutting'])
    for group, primer_pairs in groupby(sorted(primer_lookup,
                                              key=sort_key), sort_key):
        for (f_primer, r_primer_rc), name in primer_pairs:
            r_primer = reverse_complement(r_primer_rc)
            t7_primer = t7_seq + r_primer
            f_cutting_seq = 'AGTACTTAG'  # ScaI +TAG
            r_cutting_seq = 'GATGAATTC'  # GAT + EcoRI
            forward_cutting = f_primer + f_cutting_seq
            reverse_cutting = r_primer + reverse_complement(r_cutting_seq)
            if group in cut_groups:
                cv.writerow([group, name, f_primer, r_primer, t7_primer,
                             forward_cutting, reverse_cutting])
            else:
                cv.writerow([group, name, f_primer, r_primer, t7_primer])

assert (len(all_probes) <= 92918)


def primer_getter(seq):
    vals = seq[0].split(' ')
    return (vals[0], vals[-1])


discovered_primers = []
for primers, probe in groupby(sorted(all_probes), primer_getter):
    discovered_primers.append(primers)
discovered_primers = set(discovered_primers)
existing_primers = {p for p, n in set(used_primers)}

with open('temp/8-27-15.csv', 'wb') as f_out:
    f_out2 = csv.writer(f_out)
    for line in all_probes:
        f_out2.writerow(line)

# Make an archive storing all files
import zipfile
from tempfile import NamedTemporaryFile
for file in ['temp/8-27-15.zip', 'temp/8-27-15_ordering.zip']:
    with zipfile.ZipFile(file, 'w',
                         compression=zipfile.ZIP_DEFLATED) as archive:
        with NamedTemporaryFile() as f_temp:
            cv = csv.writer(f_temp)
            if 'ordering' in file:
                t_chip = [(''.join(seq.split(' ')), namel)
                          for seq, namel in all_probes]
                cv.writerows(t_chip)
            else:
                cv.writerows(all_probes)
            f_temp.flush()
            archive.write(f_temp.name, arcname='oligos' + '.csv')
        archive.write("temp/barcodes_8-27-15_brain.csv",
                      arcname='brain_barcodes.csv')
        archive.write("temp/barcodes_8-27-15_zak.csv",
                      arcname='zak_barcodes.csv')
        archive.write('temp/quats_8hybs.csv', arcname='quats_8hybs.csv')
        archive.write('temp/tertiaries_2k_genes.csv',
                      arcname='tertiaries_2k_genes.csv')
        archive.write('temp/quats_8hybs.csv', arcname='quats_8hybs.csv')
        archive.write('temp/primers_8-27-15.csv', arcname='primers.csv')
