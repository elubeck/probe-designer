import dataset
from probe_designer.probe_refiner import probe_set_refiner
import csv
import os

print(os.getcwd())
print("BOBDOLE")

purchased_10k = []
with open("db/10k_order_8-15.csv", 'r') as f_in:
    cv = csv.reader(f_in)
    for line in cv:
        purchased_10k.append(line)

bought_genes = set(k.split("-")[0] for k,v in purchased_10k)
bought_probes = set(v for k,v in purchased_10k)
print("BOD")



intron_db_filtered = dataset.connect("sqlite:///db/intron_probes_4_filtered.db")
filtered_probe_table = intron_db_filtered['mouse']
print("BOB")

all_probes = dict()
for gene_d in filtered_probe_table.distinct('target'):
    gene = gene_d['target']
    probes = [p['seq'] for p in filtered_probe_table.find(target=gene)]
    if len(probes) >= 24 and gene in bought_genes:
        all_probes[gene] = probes


all_p2 = probe_set_refiner(all_probes, block_size=18)
all_p3 = dict()
for gene, probes in all_p2.items():
    all_p3[gene] = set(probes).difference(bought_probes)

with open("temp/accessory10k_probes_8-15.csv", "w") as f_out:
    cv = csv.writer(f_out)
    for gene, probes in all_p3.items():
        for probe in probes:
            cv.writerow((gene, probe))

flat_1 = [seq for k,v in all_probes.items() for seq in v]
flat_2 = [seq for k,v in all_p2.items() for seq in v]
flat_3 = [seq for k,v in all_p3.items() for seq in v]
print("BOBDOLEO")
print("BOBDOLEO")
print("D!")
