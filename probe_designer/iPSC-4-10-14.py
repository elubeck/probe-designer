__author__ = 'Eric'

import os

import pandas as pd

from probe_designer import del_probes, get_probes, main
import get_seq


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

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[c] for c in seq][::-1])

import oligoarray_designer as od
import traceback
probe_designer = od.OligoarrayDesigner(blast_db="/home/eric/blastdb/old_format/mouse_refseq_rnaDB")
target_genes = [x.strip() for x in zak_genes]
cds_org = '"Mus musculus"[porgn:__txid10090]'
min_probes = 12
passed_genes = {}
# Get CDS and make reverse complement for oligoarray
for gene in target_genes:
    try:
        cds = get_seq.CDS(gene, cds_org).run()
        if len(''.join(cds['CDS List'])) > (35 * min_probes):
            rc_cds = [reverse_complement(c)
                                for c in cds['CDS List']]
            cds['CDS List'] = rc_cds
            passed_genes[gene] = cds
        else:
            print("CDS too short of %s at %int" % (gene, len(''.join(cds['CDS List']))))
    except get_seq.CDSError as e:
        print(traceback.format_exc())
        print(e)
pd = dict(passed_genes.items()[:10])
f_loc = "temp/temp_rna.fasta"

with open(f_loc, "w") as f:
    for gene, probes in pd.iteritems():
        for chunk_n, seq_chunk in enumerate(probes['CDS List']):
            fasta_str = ">{}=Chunk:{}\n{}\n".format(gene, chunk_n,
                                                    seq_chunk)
            f.write(fasta_str)

probes = probe_designer.run(f_loc,
            max_dist=8000,
            min_length=35,
            max_length=35,
            max_tm=100,
            cross_hyb_temp=72,
            min_tm=74,
            secondary_struct_temp=76,
            max_oligos=100,
            timeout=10,
            num_processors=1)
probes = [[probe] for probe in probes]
print("DONE OLIGOAARAY")
        
