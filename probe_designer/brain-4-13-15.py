__author__ = 'Eric'
from probe_designer import main, del_probes, get_probes
import os
import pandas as pd
genes = """
Tbr1
,Rasgrf2
,Pvrl3
,Cux2
,Rorb
,Plcxd2
,Thsd7a
,Kcnk2
,Cplx3
,Sulf2
,Foxp2
,Syt6
,Rprm
,Nr4a2
,Synpr
,Lhx5
,Reln
,Rgs8
,Cartpt
,Tle1
,Tle3
,Mdga1
,Kcnip2
,Rorb
,Cyp39a1
,Lhx2
,Unc5d
,Gpr6
,Dtx4
,Cux1
,Cux2
,Bhlhe22
,Plxnd1
,Pou3f1
,Marcksl1
,Pou3f3
,Pou3f2
,Nefh
,Cntn6
,Foxo1
,Opn3
,Lix1
,Ldb2
,Sox5
,Sema3e
,Foxp2
,Ppp1r1b
,Ctgf
,Nr4a2
,Wfs1
,DCN
,Htr2c
,Grp
,Gpr101
,Col5a1
,Gpc3
,Prss12
,Ndst4
,Calb1
,Matn2
,Rph3a
,Loxl1
,Plagl1
,Coch
,Itga7
,Lyt
,Gprin3
,Adora2a
,Cd4
,Rgs9
,Drd1
,Serpina9
,Ngfr
,Ntrk1
,Rarb
,Nts
,Pcdh11x
,Slc6a3
,Magell2
,Iyd,
Lyd,
Lm1,
Tiam2,
Robo1,
Cadps2
"""
more_genes = [u'Phactr2',
              u'Atp10a',
              u'Ndst4',
              u'Fign',
              u'Lrmp',
              u'Nts',
              u'Hpcal1',
              u'Dcc',
              u'Col23a1',
              u'Fam19a2',
              u'Esr2',
              u'Adcyap1',
              u'Esyt3',
              u'Crim1',
              u'Glra3',
              u'Itih3',
              u'Slc17a6',
              u'Gpx3',
              u'Casr',
              u'Pcdh20',
              u'Ccdc3',
              u'Cacng5',
              u'Map3k15',
              u'Sh3kbp1',
              u'Fam40b',
              u'Amigo2',
              u'Arhgef26',
              u'Pqlc1',
              u'Stard5',
              u'Sms',
              u'Coch',
              u'Ace',
              u'Acvr1c',
              u'Me2',
              u'Gpr155',
              u'Ddit4l',
              u'Dach1',
              u'Gm261',
              u'Rgs7bp',
              u'Ier5',
              u'Myo16',
              u'Vav2',
              u'Syt6',
              u'Chrna3',
              u'AF529169',
              u'Lrrc55',
              u'Tac2',
              u'Slc18a3',
              u'Ano1',
              u'Scube1',
              u'Ndst4',
              u'Nov',
              u'Ptpru',
              u'Satb2',
              u'Ltk',
              u'Mm.26272',
              u'Kcnh7',
              u'Tyro3',
              u'Cpne8',
              u'Gpr161',
              u'Loxl2',
              u'Arc',
              u'Pln',
              u'Rnf144b',
              u'Fgd5',
              u'Grm3',
              u'Esrrg',
              u'Ubash3b',
              u'Cacna1i',
              u'Trh',
              u'Plxdc2',
              u'Lypd1',
              u'Fnbp1l',
              u'Slc8a2',
              u'Cdh24',
              u'Prkar2b',
              u'Calb1',
              u'Rasal1',
              u'Col15a1',
              u'Fam131b',
              u'Fam69c',
              u'Lypd6b',
              u'Glra1',
              u'Pstpip1',
              u'Lhfpl2',
              u'C130018E23Rik*',
              u'Rgs2',
              u'Lin7a',
              u'Glra4',
              u'Etv1',
              u'Ttc6',
              u'Neto2',
              u'Kcng4',
              u'Zfp618',
              u'B3galt2',
              u'Ntng2',
              u'Cdh23',
              u'Alox8',
              u'L3mbtl4',
              u'A930038C07Rik',
              u'Iqcj',
              u'Prss35',
              u'1700019N12Rik',
              u'Iyd',
              u'Enox2',
              u'Rps6ka3',
              u'Prkcd',
              u'Tmem215',
              u'Gls',
              u'Pvrl3',
              u'Fam19a2',
              u'Lrmp',
              u'Sh3kbp1',
              u'Ndst4',
              u'Ppfibp1',
              u'Phactr2',
              u'Grem2',
              u'Sun2',
              u'Lmo3',
              u'Filip1',
              u'Gpr101',
              u'Ano2',
              u'Rcn1',
              u'Smoc1',
              u'Spint2',
              u'Ecel1',
              u'Trh',
              u'Gm10413',
              u'Glra3',
              u'Lypd1',
              u'Kcng4',
              u'Panx2',
              u'Pax6',
              u'Limk1',
              u'Ankrd34b',
              u'Ubash3b',
              u'Drd2',
              u'Cdh24',
              u'Spp1',
              u'Slc32a1',
              u'Prox1',
              u'Sipa1l2',
              u'Tdo2',
              u'Slc26a10',
              u'Fam163b',
              u'Slc39a6',
              u'Csf2rb2',
              u'LOC432928',
              u'Npnt',
              u'Btg1',
              u'Dap',
              u'Prcp',
              u'Isoc1',
              u'Rreb1',
              u'Zfpm2',
              u'D0H4S114',
              u'Btg1',
              u'LOC432748',
              u'Igfbp5',
              u'Abhd3',
              u'Cyp39a1',
              u'Sema7a',
              u'Tmem215',
              u'Rorb',
              u'Kcng4',
              u'Col15a1',
              u'Steap2',
              u'Hsd11b1',
              u'Pde5a',
              u'Glrx',
              u'Gpr101',
              u'2310042E22Rik',
              u'Dlk1',
              u'Nts',
              u'Slc17a6',
              u'Psrc1',
              u'Grp',
              u'Dcn',
              u'Prss12',
              u'Sik3',
              u'Gucy2c',
              u'Ntsr1',
              u'Clptm1l',
              u'Ntn1',
              u'Chrna6',
              u'Foxa1',
              u'Gsg1l',
              u'Bend5',
              u'Zdhhc2',
              u'Fam102b',
              u'Tac2',
              u'Gpr101',
              u'Zfhx3',
              u'Gabrg1',
              u'Prkcd',
              u'Gabrq',
              u'LOC433088',
              u'LOC381076',
              u'Plcxd2',
              u'Doc2a',
              u'Ntrk1',
              u'Kcnj5',
              u'Gpx3',
              u'Thbs2',
              u'Slc41a3',
              u'Insrr',
              u'Fam65b',
              u'Fam102b',
              u'Cntnap3',
              u'Ccdc164']


target_genes = set([x.strip() for x in genes.split(",")] + map(str, more_genes))
target_genes = list(target_genes)
print(len(target_genes))
del_probes(target_genes, organism='mouse')
probes = main(target_genes, min_probes=16, max_probes=200,
              timeout=60, debug=True, organism='mouse',
              probe_design='oligoarray')
dir_name = "passed_probes_Brain-4-13-15"
n_probes = 0
p = [(gene, masking, len(get_probes(gene, "mouse", masking)))
     for gene in target_genes for masking in [6,]]
print(p)
for gene in target_genes:
    try:
        os.mkdir(dir_name)
    except:
        pass
    passed = False
    for masking in [6,]:
        probes = get_probes(gene, "mouse", masking)
        print(gene, masking, len(probes))
        if len(probes) > 50:
            out_path = os.path.join(dir_name, gene + '_Mask={}_N={}.csv'.format(masking, len(probes)))
            pd.DataFrame(probes).to_csv(out_path)
            passed = True
    if passed:
        n_probes += 1
    if passed == False:
        print("Couldn't Find probes for {}".format(gene))
print(n_probes)

# from pathlib import Path
# from itertools import groupby
# count = 0
# max_probe = {}
# for gene, csv in groupby(Path(dir_name).glob("*.csv"),
#                          lambda x: x.name.split("_")[0]):
#     probe_files = list(csv)
#     n_probes = [int(p.name.split("N=")[1].strip(".csv")) for p in probe_files]
#     max_probes = max(n_probes)
#     max_probe[gene] = max_probes

#     if max_probes >= 50:
#         count += 1
#         print(gene,count, max_probes)
#         probes = ['Loxl1', 'Plagl1', 'Col5a1', 'Plcxd2', 'Htr2c', 'Matn2', 'Ctgf',
#                   'Marcksl1', 'Ppp1r1b', 'Opn3', 'Sema3e', 'Lhx5', 'Cyp39a1', 'Rph3a',
#                   'Bhlhe22', 'Calb1', 'Pvrl3', 'Dtx4', 'Unc5d', 'Itga7', 'Rorb',
#                   'Coch', 'Gpr101', 'Thsd7a', 'Rasgrf2', 'Reln', 'Nefh', 'Cplx3',
#                   'Lhx2', 'Foxo1', 'Sulf2', 'Tbr1', 'Pou3f3', 'DCN', 'Pou3f2', 'Cux1',
#                   'Sox5', 'Nr4a2', 'Cux2', 'Ldb2', 'Gpc3', 'Prss12', 'Pou3f1',
#                   'Kcnk2', 'Wfs1', 'Rprm', 'Mdga1', 'Syt6', 'Kcnip2', 'Tle3',
#                   'Plxnd1', 'Foxp2', 'Gpr6', 'Synpr', 'Cntn6', 'Fign', 'Atp10a',
#                   'Esyt3', 'Esr2', 'Map3k15', 'Ccdc3', 'Sema7a', 'Ppfibp1', 'Lrmp',
#                   'Rcn1', 'Ano2', 'Gpr155', 'Me2', 'Ntsr1', 'Ntn1', 'Gucy2c',
#                   'Plxdc2', 'Fnbp1l', 'Slc17a6', 'Ankrd34b', 'Panx2', 'Vav2', 'Myo16',
#                   'Chrna3', 'B3galt2', 'Neto2', 'Thbs2', 'Kcnj5', 'Ptpru', 'Satb2',
#                   'Ltk', 'Zfpm2', 'Rreb1', 'Isoc1', 'Sipa1l2', 'Prox1', 'Fam163b',
#                   'Glra1', 'Lypd6b', 'Enox2', 'Iqcj', 'Fgd5', 'Rnf144b', 'Zfhx3',
#                   'Prkcd']

#         first_round = [probe for probe in probes if probe in max_probe]