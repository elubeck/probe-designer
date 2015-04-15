__author__ = 'Eric'
genes = """
Foxa2
,Gata4
,Ptf1a
,Cdx2
,Eomes
,Gcm1
,Krt1
,Afp
,Serpina1a
,Fn1
,Lama1
,Lamb1
,Lamc1
,Sox17
,Wt1
,Des
,Myf5
,Myod1
,T
,Hba-x
,Hbb-y
,Col1a1
,Runx2
,Nes
,Neurod1
,Pax6
,Cd34
,Cdh5
,Flt1
,Pecam1
,Ddx4
,Sycp3
,Gcg
,Iapp
,Ins2
,Pax4
,Pdx1
,Sst
,Olig2
,Tat
,Foxd3
,Gata6
,Gbx2
,Nanog
,Nr5a2
,Nr6a1
,Pou5f1
,Sox2
,Tcfcp2l1
,Utf1
,Zfp42
,Commd3
,Crabp2
,Ednrb
,Fgf4
,Fgf5
,Gabrb3
,Gal
,Grb7
,Hck
,Ifitm1
,Il6st
,Kit
,Lefty1
,Lefty2
,Lifr
,Nodal
,Nog
,Numb
,Pten
,Sfrp2
,Tdgf1
,Fgf4
,Fgf5
,Gdf3
,Lefty1
,Lefty2
,Nodal
,Brix1
,Cd9
,Diap2
,Ifitm2
,Igf2bp2
,Lin28a
,Podxl
,Rest
,Sema3a
,Tert
,Il6st
,Lifr
,Socs3
,Stat3
,Fgfr1
,Fgfr2
,Fgfr3
,Fgfr4
,Ptch1
,Ptchd2
,Smo
,Gli1
,Gli2
,Gli3
,Sufu
,Ncstn
,Notch1
,Notch2
,Notch3
,Notch4
,Psen1
,Psen2
,Psenen
,Rbpjl
,Hes1
,Hey1
,Acvr1
,Acvr1b
,Acvr1c
,Acvr2a
,Acvr2b
,Acvrl1
,Amhr2
,Bmpr1a
,Bmpr1b
,Bmpr2
,Eng
,Ltbp1
,Ltbp2
,Ltbp3
,Ltbp4
,Rgma
,Tgfbr1
,Tgfbr2
,Tgfbr3
,Tgfbrap1
,Crebbp
,E2f5
,Ep300
,Rbl1
,Rbl2
,Smad1
,Smad2
,Smad3
,Smad4
,Smad5
,Smad6
,Smad7
,Smad9
,Sp1
,Zeb2
,Fzd1
,Fzd2
,Fzd3
,Fzd4
,Fzd5
,Fzd6
,Fzd7
,Fzd8
,Fzd9
,Lrp5
,Lrp6
,Vangl2
,Bcl9
,Bcl9l
,Ctnnb1
,Lef1
,Nfat5
,Nfatc1
,Nfatc2
,Nfatc3
,Nfatc4
,Pygo2
,Tcf7l1
,Tcf7
,Tcf7l2
"""

from probe_designer import main, del_probes, get_probes
import os
import pandas as pd

target_genes = [x.strip() for x in genes.split(",")]
del_probes(target_genes, organism='mouse')
probes = main(target_genes, min_probes=16, max_probes=200,
              timeout=60, debug=True, organism='mouse',
              probe_design='oligoarray')
dir_name = "passed_probes_iPSCc_oligoarray"
n_probes = 0
p = [(gene, masking, len(get_probes(gene, "mouse", masking)))
     for gene in target_genes for masking in [6]]
print(p)
for gene in target_genes:
    try:
        os.mkdir(dir_name)
    except:
        pass
    passed = False
    for masking in [6]:
        probes = get_probes(gene, "mouse", masking)
        print(gene, masking, len(probes))
        if len(probes) > 16:
            out_path = os.path.join(dir_name, gene + '_Mask={}_N={}.csv'.format(masking, len(probes)))
            pd.DataFrame(probes).to_csv(out_path)
            passed = True
    if passed:
        n_probes += 1
    if passed == False:
        print("Couldn't Find probes for {}".format(gene))
print(n_probes)
