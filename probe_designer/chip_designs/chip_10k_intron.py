"""
This file was used to design probes for 10k genes. 
Both the old database from the previous intron design project and a new database designing more liberal intron probes were used.  These databases were generated on whole genome runs using the intron_designer2.design_introns() function. 
All subsequent steps are documented in this file.
"""

from __future__ import print_function, division
import dataset
import probe_designer.intron_designer2
from collections import defaultdict
from itertools import groupby
import probe_designer.blaster2
from progressbar import ProgressBar
import probe_designer.mRNA_designer
import probe_designer.blaster2
import csv
import json
from probe_designer.get_seq import reverse_complement
from collections import Counter

filterer = probe_designer.intron_designer2.ProbeFilter()
###### Blast Adapters ###########
adapters = []
with open("db/adapter10k.csv", "r") as f_in:
    for line in f_in:
        adapters.append(line.strip("\n"))

fasta_adapter = "\n".join(">{}\n{}".format(n, probe)
                          for n, probe in enumerate(adapters))
blast_hits_adapter = probe_designer.blaster2.local_blast_query(fasta_adapter,
                                                db='gencode_tracks_reversed')
blast_res_adapter = probe_designer.blaster2.parse_hits(blast_hits_adapter,
                                        match_thresh=10,
                                        strand=1)
off_target = {
    bridge: filterer.get_copynum(hits)
    for bridge, hits in blast_res_adapter.iteritems()
}

####### Filter New and old Probes  #############
dbs = ["sqlite:///db/intron_probes3.db", "sqlite:///db/intron_probes2.db.bk"]
filterer = probe_designer.intron_designer2.ProbeFilter()
intron_db_filtered = dataset.connect(
    "sqlite:///db/intron_probes_filtered_10k.db")
filtered_probe_table = intron_db_filtered['mouse']

# Commented out to save runtime
# for db in dbs:
#     intron_db = dataset.connect(db)
#     probes = intron_db['mouse']
#     genes = [g['Name'] for g in probes.distinct('Name')]
#     used_probes = {
#         p['target']
#         for p in filtered_probe_table.distinct('target')
#     }
#     remaining_genes = set(genes).difference(used_probes)
#     p_bar = ProgressBar(maxval=len(remaining_genes))
#     for n, gene in enumerate(remaining_genes):
#         if gene in used_probes: continue
#         res = [r["Probe (5'-> 3')"] for r in probes.find(Name=gene)]
#         if len(res) >= 20:
#             f_probes = [{'target': gene,
#                          'seq': p} for p in filterer.run(res, gene,
#                                                          n_probes=25)]
#             filtered_probe_table.insert_many(f_probes)
#         p_bar.update(n)
#     p_bar.finish()

# Get Genes with >=24 probes
all_probes = {}
for gene_d in filtered_probe_table.distinct('target'):
    gene = gene_d['target']
    probes = [p['seq'] for p in filtered_probe_table.find(target=gene)]
    if len(probes) >= 24:
        all_probes[gene] = probes

# Search for redundant nested sequences
p_set2 = probe_designer.mRNA_designer.probe_set_refiner(all_probes)

# Blast Everything to check that all probes together don't cause big problems
flat_probes = [
    {"gene": gene,
     "seq": probe,
     "number": n,
     "name": "{}-{}".format(gene, n)}
    for gene, probes in p_set2.iteritems() for n, probe in enumerate(probes)
]
# Commented out to save runtime
# fasta_str = "\n".join(">{name}\n{seq}".format(**probe)
#                       for probe in flat_probes)
# blast_hits2 = blaster2.local_blast_query(fasta_str,
#                                          db='gencode_tracks_reversed')
# blast_res2 = blaster2.parse_hits(blast_hits2, match_thresh=20, strand=1)

# with open('temp/chip_10k_intron.json', 'w') as f_out:
#     json.dump(dict(blast_res2), f_out)

# Load from file to save time
with open('temp/chip_10k_intron.json', 'r') as f_out:
    blast_res2 = json.load(f_out)

# For every gene hit in blast results see what probes hit it
hit_to_probe = defaultdict(list)
for probe, hits in blast_res2.iteritems():
    for hit in hits:
        hit_to_probe[hit].append(probe)

# Get probes that hit unintended target
bad_hits = [(hit, [probe for probe in probes
                   if probe.split('-')[0] not in hit])
            for hit, probes in hit_to_probe.iteritems()]

bad_hits = sorted([(k, hits) for k, hits in bad_hits if any(hits)],
                  key=lambda x: len(x[1]),
                  reverse=True)

# Drop probes if more than thresh probes hit x-target or copy number is greater than c_thresh
thresh = 5
c_thresh = 50
drop_list = []
for hit, probes in bad_hits:
    ensembl_name = hit.split(',')[0]
    if len(probes) >= thresh or filterer.get_copynum(
        [ensembl_name]) > c_thresh:
        for probe in probes:
            drop_list.append(probe)

drop_probes = set(drop_list)

# Make final list of good probes
flat_probes_filtered = [probe for probe in flat_probes
                        if probe['name'] not in drop_probes]

# Count remaining probes
grouped_probes = [(group, list(probes))
                  for group, probes in groupby(flat_probes_filtered,
                                               key=lambda x: x['gene'])]
passed = [g for g, probes in grouped_probes if len(probes) == 25]

# Make x-refrence between refseq -> to encode names
import json
with open('db/encode_to_refseq.json', 'r') as f_in:
    enc_ref = json.load(f_in)
ref_enc = defaultdict(list)
for k, v in enc_ref.iteritems():
    ref_enc[v].append(k)

# Sort Genes by Copy Number
c_num = {
    gene: filterer.get_copynum(ref_enc[gene])
    for gene in passed if gene in ref_enc
}

c_sort = sorted(c_num.items(), key=lambda x: x[1], reverse=True)

###### Make Probes #########
# Start running form here to shorten evaluation time on already run file
seq_template = "{f_primer} {f_cutting} TAG {bridge} TATA {probe} GAT {r_cutting} {r_primer}"
f_cutting = 'AGTACT'  # ScaI
r_cutting = 'GAATTC'  # EcoRI

primers = [
    ['AAGCGCCACGAGTTGTCACG', 'GCATCACTTTGGGCTCGGCT'],
    ['GCCCCATCATGTGCCTTTCC', 'GCCTGCGAGAGGGAAATCCA'],
    ['TGTGCGCTCCGATTGTCCTC', 'GGCCAACAGACCCCATTTGC'],
    ['ATTGAGGGTCTTCGCGTGCC', 'GGTTGCAAAGCGCCGGTTAC'],
    ['AATTGAGCAGCTCGGGCCAC', 'AGTTGCAGGCTTCCATCGCC'],
]  # yapf: disable

used_bridges = []

# Already picked probes from last probe set
existing_probes = defaultdict(list)
with open('temp/5-27-15_intron_probes2.csv', 'r') as f_in:
    for probe, gene in csv.reader(f_in, delimiter='\t'):
        probe_seq = probe[53:88]
        bridge = probe[29:49]
        used_bridges.append((reverse_complement(bridge), gene))
        existing_probes[gene].append({
            'probe': probe_seq,
            'bridge': bridge,
            'f_primer': primers[0][0],
            'r_primer': reverse_complement(primers[0][1]),
            'f_cutting': f_cutting,
            'r_cutting': r_cutting,
        })

# Remove Already picked probes from remaining choices
remaining_choices = set(passed).difference(existing_probes.keys())

# Get highly expressed probes
remaining_c_sort = [k for k, v in c_sort if k in remaining_choices]
high_expression = [k for k in remaining_c_sort[:1000]]

# Get remaining probe sets
from random import shuffle
remaining_choices = list(set(remaining_c_sort).difference(high_expression))
shuffle(remaining_choices)
set_size = int(len(remaining_choices) / 3)
probe_blocks = [[
    g for g in remaining_choices[i * set_size:(i + 1) * set_size]
] for i in range(3)]

######### Sort and Choose Adapters #################

# Sort adapters by copy number of hits and make iterable
ub = set(map(lambda x: x[0], used_bridges))
br2 = sorted([(k, v) for k, v in off_target.items()
              if adapters[int(k)] not in ub],
             key=lambda x: x[1], )
iter_bridge = iter([adapters[int(k)] for k, v in br2])

############# Make Chips ####################

probe_dict = {g: probes for g, probes in grouped_probes if len(probes) == 25}

# First make high expression set
high_exp = defaultdict(list)
for gene in high_expression:
    bridge = iter_bridge.next()
    used_bridges.append((bridge, gene))
    for probe in probe_dict[gene]:
        high_exp[gene].append({
            'probe': probe['seq'],
            'bridge': reverse_complement(bridge),
            'f_primer': primers[3][0],
            'r_primer': reverse_complement(primers[3][1]),
            'f_cutting': f_cutting,
            'r_cutting': r_cutting,
        })

# Make Remaining Sets
meta_p_block = []
for probe_set in probe_blocks:
    p_block = defaultdict(list)
    for gene in probe_set:
        bridge = iter_bridge.next()
        used_bridges.append((bridge, gene))
        for probe in probe_dict[gene]:
            p_block[gene].append({
                'probe': probe['seq'],
                'bridge': reverse_complement(bridge),
                'f_primer': primers[0][0],
                'r_primer': reverse_complement(primers[0][1]),
                'f_cutting': f_cutting,
                'r_cutting': r_cutting,
            })
    meta_p_block.append(p_block)

bridge, gene = zip(*used_bridges)

# Check for bugs
assert (len(set(used_bridges)) == len(set(bridge)) == len(set(gene)))

###### Barcode the Bridges ########

# Get integer barcodes
with open('db/barcode10k.csv', 'r') as f_in:
    replace = {1: 4, 2: 2, 3: 3, 4: 1}  # Swap 1 and 4
    barcodes = {
        tuple(map(lambda x: replace[int(x)], line))
        for line in csv.reader(f_in)
    }

bridge_template = "{f_primer} GGTACC TCC {initiator} TATA {adapter} CTA AGTACT {r_primer}"

initiators = ["gCATTCTTTCTTgAggAgggCAgCAAACgggAAgAg",
              "AgCTCAgTCCATCCTCgTAAATCCTCATCAATCATC",
              "AAAgTCTAATCCgTCCCTgCCTCTATATCTCCACTC",
              "CACATTTACAgACCTCAACCTACCTCCAACTCTCAC"]
hyb_barcodes = zip(*barcodes)

all_bridges = []
bridge_primers = [primers[1], primers[2], primers[1], primers[2], primers[1],
                  primers[2], primers[1], primers[2]]

# Make the bridges
for n, barcode_round in enumerate(hyb_barcodes):
    bridge_set = []
    for (bridge, gene), init_n in zip(set(used_bridges), barcode_round):
        f_primer, r_primer = bridge_primers[n]
        params = {
            "adapter": bridge,
            "f_primer": f_primer,
            "r_primer": reverse_complement(r_primer),
            "initiator": initiators[init_n - 1]
        }
        seq = bridge_template.format(**params)
        name = "{}-{}-{}-Bridge".format(gene, n, init_n)
        bridge_set.append((seq, name))
    all_bridges.append(bridge_set)

# Write barcodes to file
with open("temp/barcodes_6-24-15.csv", "w") as f_out:
    cv = csv.writer(f_out)
    cv.writerow(['name', 'bridge', 'Hyb1', 'Hyb2', 'Hyb3', 'Hyb4', 'Hyb5',
                 'Hyb6', 'Hyb7', 'Hyb8'])
    for (bridge, name), barcode in zip(set(used_bridges), barcodes):
        cv.writerow([name, bridge] + list(barcode))

############## Sort Probes into Chips ################
seq_template = "{f_primer} {f_cutting} TAG {bridge} TATA {probe} GAT {r_cutting} {r_primer}"
f_cutting = 'AGTACT'  # ScaI
r_cutting = 'GAATTC'  # EcoRI

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
    pgk1_bridge = iter_bridge.next()
    used_bridges.append((pgk1_bridge, "Split PGK1 {}-Bridge".format(cn)))
    for n, probe in enumerate(color_block):
        split_pgk1_template.append(({
            'f_cutting': f_cutting,
            'bridge': reverse_complement(pgk1_bridge),
            'probe': probe,
            'r_cutting': r_cutting,
        }, 'Split PGK1-Color{}-#{}'.format(cn, n)))

pgk1_sub_set = ['Split PGK1-Color0-#1', 'Split PGK1-Color0-#7',
                'Split PGK1-Color0-#12', 'Split PGK1-Color0-#18',
                'Split PGK1-Color0-#23', 'Split PGK1-Color1-#2',
                'Split PGK1-Color1-#6', 'Split PGK1-Color1-#11',
                'Split PGK1-Color1-#15', 'Split PGK1-Color1-#17']

pgk1_short = [(d, name) for d, name in split_pgk1_template
              if name in pgk1_sub_set]


def append_pgk1(chip, pgk1_template, ):
    """
    Adds PGK1 oligos to chip.  One set for every primer pair.
    """
    # Get primers for PGK1
    f_primers, r_primers = [], []
    for oligo, name in chip:
        sv = oligo.split(' ')
        f_primers.append(sv[0])
        r_primers.append(sv[-1])

    primer_pairs = set(zip(f_primers, r_primers))
    for f_primer, r_primer, in primer_pairs:
        primer_dict = {'f_primer': f_primer, 'r_primer': r_primer}
        for template, name in pgk1_template:
            template.update(primer_dict)
            chip.append((seq_template.format(**template), name))
    return chip

# Control Oligo
control_seq = seq_template.format(**{
    'f_primer': 'GACGCACATATGCGGGCAAG',
    'f_cutting': f_cutting,
    'r_cutting': r_cutting,
    'bridge': 'ATGCATGCATGCATGCATGC',
    'probe': 'GCTTGCAAGCTTGCAAGCTTGCAAGCTTGCAAGCA',
    'r_primer': reverse_complement('TCCGCAGTCACGAAGATGCC')
})

# Chip 1
chip_primers = []
chip_1 = []
from random import sample
for gene, probes in existing_probes.iteritems():
    if len(probes) >= 47:  # 48 probes filled the chip
        samp_size = 47
    else:
        samp_size = len(probes)
    for n, probe in enumerate(sample(probes, samp_size)):
        chip_1.append((seq_template.format(**probe), "{}-{}".format(gene, n)))

chip_primers.append(
    (chip_1[-1][0].split(' ')[0], chip_1[-1][0].split(' ')[-1]))

for gene, probes in high_exp.iteritems():
    for n, probe in enumerate(probes):
        chip_1.append((seq_template.format(**probe), "{}-{}".format(gene, n)))

chip_primers.append(
    (chip_1[-1][0].split(' ')[0], chip_1[-1][0].split(' ')[-1]))

chip_1 += all_bridges[0]

chip_primers.append(
    (chip_1[-1][0].split(' ')[0], chip_1[-1][0].split(' ')[-1]))

chip_1 += all_bridges[1]

chip_primers.append(
    (chip_1[-1][0].split(' ')[0], chip_1[-1][0].split(' ')[-1]))

chip_1 = append_pgk1(chip_1, pgk1_short)
chip_1 += [(control_seq, 'Control-Exclude')]

chip_primers.append(
    (chip_1[-1][0].split(' ')[0], chip_1[-1][0].split(' ')[-1]))
assert (len(chip_1) <= 92918)
assert (all(v == 1 for k, v in Counter(chip_primers).iteritems()))

# Chip 2
chip_primers = []
chip_2 = []
for gene, probes in meta_p_block[0].iteritems():
    for n, probe in enumerate(probes):
        chip_2.append((seq_template.format(**probe), "{}-{}".format(gene, n)))
chip_primers.append(
    (chip_2[-1][0].split(' ')[0], chip_2[-1][0].split(' ')[-1]))
chip_2 += all_bridges[2]
chip_primers.append(
    (chip_2[-1][0].split(' ')[0], chip_2[-1][0].split(' ')[-1]))
chip_2 += all_bridges[3]
chip_primers.append(
    (chip_2[-1][0].split(' ')[0], chip_2[-1][0].split(' ')[-1]))

# Get Ahmet's Probes
import random
ahmet_primers = ['TGCAGCTCCGCGAAATGAAG',
                 reverse_complement('AATGGCACAGACAGGCAGCG')]
ahmet_probes = []
ahmet_template = "{f_primer} {f_cutting} TAG {seq} GAT {r_cutting} {r_primer}"
f_cutting = 'AGTACT'  # ScaI
r_cutting = 'GAATTC'  # EcoRI

with open('temp/ahmet_riboprobes.csv', 'r') as f_in:
    for line in csv.reader(f_in):
        seq = line[1]
        full_seq = ahmet_template.format(**{
            'f_primer': ahmet_primers[0],
            'f_cutting': f_cutting,
            'r_primer': ahmet_primers[1],
            'r_cutting': r_cutting,
            'seq': seq
        })
        ahmet_probes.append((full_seq, line[0] + '-Ahmet'))

remaining_space = 92918 - len(chip_2) - 48 * 4 - 1
chip_2 += sample(ahmet_probes, remaining_space)
chip_primers.append(
    (chip_2[-1][0].split(' ')[0], chip_2[-1][0].split(' ')[-1]))

chip_2 = append_pgk1(chip_2, split_pgk1_template)
chip_2 += [(control_seq, 'Control-Exclude')]
chip_primers.append(
    (chip_2[-1][0].split(' ')[0], chip_2[-1][0].split(' ')[-1]))

assert (len(chip_2) <= 92918)
assert (all(v == 1 for k, v in Counter(chip_primers).iteritems()))

# Chip 3
chip_primers = []
# Get 100 bridges for Zak's Genes
zak_bridges = []
for bridge_file in ('temp/1k Bridges #10.csv', 'temp/1k Bridges #11.csv', ):
    with open(bridge_file, 'r') as f_in:
        for n, line in enumerate(csv.reader(f_in)):
            if n == 0:
                continue
            zak_bridges.append(reverse_complement(line[2][5:]))

iter_zak_bridge = iter(zak_bridges[:100])

zak_oligos = []
with open("temp/zak_priority_genes.csv", "r") as f_in:
    for line in csv.reader(f_in):
        zak_oligos.append(line[1:])

# Get Zak's Probes
final_zak_oligos = []
used_zak_bridges = []
for gene, oligos in groupby(zak_oligos, lambda x: x[0]):
    adapter = iter_zak_bridge.next()
    used_zak_bridges.append((reverse_complement(adapter), gene))
    for n, (g_name, oligo) in enumerate(oligos):
        params = {
            'probe': oligo,
            'bridge': adapter,
            'f_primer': primers[4][0],
            'r_primer': reverse_complement(primers[4][1]),
            'f_cutting': f_cutting,
            'r_cutting': r_cutting,
        }
        final_zak_oligos.append(
            (seq_template.format(**params), "{}-{}-ZAK".format(g_name, n)))

chip_3 = []
for gene, probes in meta_p_block[1].iteritems():
    for n, probe in enumerate(probes):
        chip_3.append((seq_template.format(**probe), "{}-{}".format(gene, n)))
chip_primers.append(
    (chip_3[-1][0].split(' ')[0], chip_3[-1][0].split(' ')[-1]))
chip_3 += all_bridges[4]
chip_primers.append(
    (chip_3[-1][0].split(' ')[0], chip_3[-1][0].split(' ')[-1]))
chip_3 += all_bridges[5]
chip_primers.append(
    (chip_3[-1][0].split(' ')[0], chip_3[-1][0].split(' ')[-1]))

# Add Zak Oligos
remaining_space = 92918 - len(chip_3) - 48 * 4 - 1
chip_3 += final_zak_oligos[:remaining_space]
chip_primers.append(
    (chip_3[-1][0].split(' ')[0], chip_3[-1][0].split(' ')[-1]))

remaining_zak = final_zak_oligos[remaining_space:
                                 ]  # Remaining oligos to add to chip 4

chip_3 = append_pgk1(chip_3, split_pgk1_template)
chip_3 += [(control_seq, 'Control-Exclude')]
chip_primers.append(
    (chip_3[-1][0].split(' ')[0], chip_3[-1][0].split(' ')[-1]))

remaining_space = 92918 - len(chip_3)
assert (len(chip_3) <= 92918)
assert (all(v == 1 for k, v in Counter(chip_primers).iteritems()))

# Chip 4
chip_4 = []
chip_primers = []
for gene, probes in meta_p_block[2].iteritems():
    for n, probe in enumerate(probes):
        chip_4.append((seq_template.format(**probe), "{}-{}".format(gene, n)))
chip_primers.append(
    (chip_4[-1][0].split(' ')[0], chip_4[-1][0].split(' ')[-1]))
chip_4 += all_bridges[6]
chip_primers.append(
    (chip_4[-1][0].split(' ')[0], chip_4[-1][0].split(' ')[-1]))
chip_4 += all_bridges[7]
chip_primers.append(
    (chip_4[-1][0].split(' ')[0], chip_4[-1][0].split(' ')[-1]))

chip_4 += remaining_zak  # Add Zak Oligos
chip_primers.append(
    (chip_4[-1][0].split(' ')[0], chip_4[-1][0].split(' ')[-1]))

# Add failed bridges from last chip
old_bridges_4 = []
with open('temp/bridge_set4.csv', 'r') as f_in:
    for n, bridge in enumerate(csv.reader(f_in)):
        seq = "{} {} {}".format(ahmet_primers[0], bridge[0][20:-20],
                                ahmet_primers[1])
        old_bridges_4.append((seq, "OldBridge4-{}".format(n)))

chip_4 += old_bridges_4
chip_primers.append(
    (chip_4[-1][0].split(' ')[0], chip_4[-1][0].split(' ')[-1]))
remaining_space = 92918 - len(chip_4)

chip_4 = append_pgk1(chip_4, pgk1_short)

chip_4 += [(control_seq, 'Control-Exclude')]
chip_primers.append(
    (chip_4[-1][0].split(' ')[0], chip_4[-1][0].split(' ')[-1]))
remaining_space = 92918 - len(chip_4)
assert (len(chip_4) <= 92918)
assert (all(v == 1 for k, v in Counter(chip_primers).iteritems()))

# Make an archive storing all files
import zipfile
from tempfile import NamedTemporaryFile
for file in ['temp/chip_10k_intron.zip', 'temp/chip_10k_intron_ordering.zip']:
    with zipfile.ZipFile(file, 'w',
                         compression=zipfile.ZIP_DEFLATED) as archive:
        for name, chip in {
            'chip_1': chip_1,
            'chip_2': chip_2,
            'chip_3': chip_3,
            'chip_4': chip_4
        }.iteritems():
            with NamedTemporaryFile() as f_temp:
                cv = csv.writer(f_temp)
                if 'ordering' in file:
                    t_chip = [(''.join(seq.split(' ')), namel)
                              for seq, namel in chip]
                    cv.writerows(t_chip)
                else:
                    cv.writerows(chip)
                f_temp.flush()
                archive.write(f_temp.name, arcname=name + '.csv')
        archive.write("temp/barcodes_6-24-15.csv", arcname='barcodes.csv')
        with NamedTemporaryFile() as f_temp:
            cv = csv.writer(f_temp)
            cv.writerow(['bridge', 'name'])
            cv.writerows(set(used_bridges))
            f_temp.flush()
            archive.write(f_temp.name, arcname='bridges.csv')
        with NamedTemporaryFile() as f_temp:
            cv = csv.writer(f_temp)
            cv.writerow(['bridge', 'name'])
            cv.writerows(set(used_zak_bridges))
            f_temp.flush()
            archive.write(f_temp.name, arcname='zak_bridges.csv')

######### BLAT ALL PROBES ###############
import probe_designer.align_probes
intron_db_filtered = dataset.connect(
    "sqlite:///db/intron_probes_filtered_10k.db")
filtered_probe_table = intron_db_filtered['mouse']

probe_alignments = dataset.connect("sqlite:///db/intron_probes_aligned_10k.db")
mouse_alignments = probe_alignments['mouse']
n_genes = len(list(filtered_probe_table.distinct('target')))
p_bar = ProgressBar(maxval=n_genes)
for n, target_d in enumerate(filtered_probe_table.distinct('target')):
    gene = target_d['target']
    passed = [probe['seq'] for probe in filtered_probe_table.find(target=gene)]
    chr = probe_designer.align_probes.find_chromosome(gene)
    chr_positions = probe_designer.align_probes.blat_probes(passed, chr)
    for pos in chr_positions:
        pos.update({'target': gene})
    mouse_alignments.insert_many(chr_positions)
    p_bar.update(n)
p_bar.finish()

with open('temp/intron_10k.bed', 'w') as f_out:
    bed = csv.writer(f_out, delimiter='\t')
    bed.writerow(['chrom', 'chromStart', 'chromEnd', 'strand', 'name'])
    for row in mouse_alignments:
        bed.writerow([row['tName'], row['tStart'], row['tEnd'], row['strand'],
                      row['target'] + '-' + row['qName']])

###### PRIMER LIST #######

primers = [
    ['AAGCGCCACGAGTTGTCACG', 'GCATCACTTTGGGCTCGGCT'],
    ['GCCCCATCATGTGCCTTTCC', 'GCCTGCGAGAGGGAAATCCA'],
    ['TGTGCGCTCCGATTGTCCTC', 'GGCCAACAGACCCCATTTGC'],
    ['ATTGAGGGTCTTCGCGTGCC', 'GGTTGCAAAGCGCCGGTTAC'],
    ['AATTGAGCAGCTCGGGCCAC', 'AGTTGCAGGCTTCCATCGCC'],
]  # yapf: disable

ahmet_primers = ['TGCAGCTCCGCGAAATGAAG', 'AATGGCACAGACAGGCAGCG']
control_primer = ['GACGCACATATGCGGGCAAG', 'TCCGCAGTCACGAAGATGCC']
