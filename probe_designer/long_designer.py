from __future__ import print_function
from probe_designer import get_probes
from itertools import groupby
import pandas as pd
from random import getrandbits
from probe_designer import batch_blast_probes

with open("/home/eric/cds.fasta", 'r') as f1:
    names = []
    for line in f1:
        if line.startswith(">"):
            gene_name = line.strip(">").split("$")[0]
            names.append(gene_name)

probes = []
with open("/home/eric/chipdesign0415.txt", 'r') as f:
    for n, line in enumerate(f):
        sub_line = line.strip("{}").split(',')[1:]
        for line in sub_line:
            cds = line.strip(' ""')
            cds = cds.strip('}\r\n"')
            probes.append((names[n], cds))

longs_probes = {}
for name, group in groupby(probes, lambda x: x[0]):
    longs_probes[name] = [v for k,v in group]

probe_df = pd.DataFrame(probes, columns=['Name', "Probe (5'-> 3')"])
probe_df['Percent GC'] = ["NA" for i in range(len(probe_df))]
probe_df['Probe Position*']= [getrandbits(10) for i in range(len(probe_df))]
probe_df['Masking'] = [632 for i in range(len(probe_df))]
# probe_df['CDS Region #'] = [i for i in range(len(probe_df))]
    

cds_org = '"Mus musculus"[porgn:__txid10090]'
organism = "mouse"
debug = False
debug = True
max_probes = 200
min_probes = 12
timeout = 1800
bob = [[{'Name':k, 'Masking': 632, 'Probes':v, 'CDS Region #': getrandbits(5)}]
       for k, v in probe_df.groupby("Name")]
# g_set = batch_blast_probes(cds_org, debug, max_probes, min_probes,
#                            organism, bob, timeout, local=True)
passed_probes = []
passed_names = set()
for probe in bob:
    probe = probe[0]
    probes = get_probes(probe['Name'], "mouse", probe['Masking'])
    if any(probes):
        print(probe['Name'], len(probes))
        passed_probes.append(pd.DataFrame(probes))
        passed_names.add(probe['Name'])
big_df = pd.concat(passed_probes)
big_df.to_csv("35mer_probes.csv")
i = ['Loxl1', 'Plagl1', 'Col5a1', 'Plcxd2', 'Htr2c', 'Matn2', 'Ctgf',
     'Marcksl1', 'Ppp1r1b', 'Opn3', 'Sema3e', 'Lhx5', 'Cyp39a1', 'Rph3a',
     'Bhlhe22', 'Calb1', 'Pvrl3', 'Dtx4', 'Unc5d', 'Itga7', 'Rorb',
     'Coch', 'Gpr101', 'Thsd7a', 'Rasgrf2', 'Reln', 'Nefh', 'Cplx3',
     'Lhx2', 'Foxo1', 'Sulf2', 'Tbr1', 'Pou3f3', 'DCN', 'Pou3f2', 'Cux1',
     'Sox5', 'Nr4a2', 'Cux2', 'Ldb2', 'Gpc3', 'Prss12', 'Pou3f1',
     'Kcnk2', 'Wfs1', 'Rprm', 'Mdga1', 'Syt6', 'Kcnip2', 'Tle3',
     'Plxnd1', 'Foxp2', 'Gpr6', 'Synpr', 'Cntn6', 'Fign', 'Atp10a',
     'Esyt3', 'Esr2', 'Map3k15', 'Ccdc3', 'Sema7a', 'Ppfibp1', 'Lrmp',
     'Rcn1', 'Ano2', 'Gpr155', 'Me2', 'Ntsr1', 'Ntn1', 'Gucy2c',
     'Plxdc2', 'Fnbp1l', 'Slc17a6', 'Ankrd34b', 'Panx2', 'Vav2', 'Myo16',
     'Chrna3', 'B3galt2', 'Neto2', 'Thbs2', 'Kcnj5', 'Ptpru', 'Satb2',
     'Ltk', 'Zfpm2', 'Rreb1', 'Isoc1', 'Sipa1l2', 'Prox1', 'Fam163b',
     'Glra1', 'Lypd6b', 'Enox2', 'Iqcj', 'Fgd5', 'Rnf144b', 'Zfhx3',
            'Prkcd']
for gene in i:
    if gene not in passed_names:
        print("SAD", gene)
