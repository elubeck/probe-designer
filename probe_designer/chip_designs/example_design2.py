gene_str = """
ABCG2; NeuroD1; ASCL1/Mash1; Noggin; Beta-catenin; Notch-1; Notch-2; Brg1  ; Nrf2  ; N-Cadherin; Nucleostemin; Calcitonin R; Numb; CD15/Lewis X; Otx2; CDCP1; Pax3; COUP-TF I/NR2F1; Pax6; CXCR4; PDGF R alpha; FABP7/B-FABP; PKC zeta; FABP 8/M-FABP; Prominin-2; FGFR2; ROR2; FGFR4; RUNX1/CBFA2; FoxD3; RXR alpha/NR2B1; Frizzled-9; sFRP-2; GATA-2; SLAIN 1; GCNF/NR6A1; SOX1; GFAP; SOX2; Glut1; SOX9; HOXB1; SOX11; ID2; SOX21; Meteorin; SSEA-1; MSX1; TRAF-4; Musashi-1; Vimentin; Musashi-2; ZIC1; Nestin

A2B5; AP-2 Alpha; ATPase Na+/K+ transporting alpha 1; Activin RIIA; Brg1; CD168/RHAMM; CD4; Doublecortin/DCX; Frizzled 4/CD344; GAP43; Jagged1; Laminin; MSX1/HOX7; Mash1; Musashi-1; Nestin; Netrin-1; Netrin-4; Neuritin; NeuroD1; Neurofilament alpha-internexin/NF66; Notch1; Notch2; Notch3; Nucleostemin; Otx2; PAX3; S100B; SOX2; Semaphorin 3C; Semaphorin 6A; Semaphorin 6B; Semaphorin 7A; TROY/TNFRSF19; Tubulin BII; Tuj 1; Vimentin

ATH1/MATH1; ASH1/MASH1; HES5; HuC/Hu; HuD; Internexin alpha; L1 neural adhesion molecule; MAP1B/MAP5; MAP2A; MAP2B; Nerve Growth Factor Rec/NGFR); Nestin; NeuroD; Neurofilament L 68 kDa; Neuron Specific Enolase/NSE; NeuN; Nkx-2.2/NK-2; Noggin; Pax-6; PSA-NCAM; Tbr1; Tbr2; Tubulin beta III; TUC-4; Tyrosine hydroxylase/TH

Collapsin Response Mediated Protein 1 /CRMP1; Collapsin Response Mediated Protein 2 /CRMP2; Collapsin Response Mediated Protein 5 /CRMP5; Contactin-1; Cysteine-rich motor neuron 1/CRIM1; c-Ret phosphor Serine 696; Doublecortin/DCX; Ephrin A2; Ephrin A4; Ephrin A5; Ephrin B1; Ephrin B2; GAP-43; HuC; HuD; Internexin alpha; Laminin-1; LINGO-1; MAP1B/MAP5; Mical-3; NAP-22; NGFR; Nestin; Netrin-1; Neuropilin; Plexin-A1; RanBPM; Semaphorin 3A; Semaphorin 3F; Semaphorin 4D; Slit2; Slit3; Staufen; Tbr 1; Tbr 2; Trk A; Tubulin BIII; TUC-4

NeuN; NF-L; NF-M; GAD; TH; PSD-95; Synaptophysin; VAMP; ZENON

4.1G; Acetylcholinesterase; Ack1; AMPA Receptor Binding Protein/ABP; ARG3.1; Arp2; E-Cadherin; N-Cadherin; Calcyon; Catenin alpha and beta; Caveolin; CHAPSYN-110/PSD93; Chromogranin A; Clathrin light chain; Cofilin; Complexin 1/CPLX1/Synaphin 2; Contactin-1; CRIPT; Cysteine String Protein/CSP; Dynamin 1; Dymanin 2; Flotillin-1; Fodrin; GRASP; GRIP1; Homer; Mint-1; Munc-18; NSF; PICK1; PSD-95; RAB4; Rabphillin 3A; SAD A; SAD B; SAP-102; SHANK1a; SNAP-25; Snapin; Spinophilin/Neurabin-1; Stargazin; Striatin; SYG-1; Synaptic Vesicle Protein 2A; Synaptic Vesicle Protein 2B; Synapsin 1; Synaptobrevin/VAMP; Synaptojanin 1; Synaptophysin; Synaptotagmin; synGAP; Synphilin-1; Syntaxin 1; Syntaxin 2; Syntaxin 3; Syntaxin 4; Synuclein alpha; VAMP-2; Vesicular Acetylcholine Transporter/VAChT; Vesicular GABA transporter/VGAT/VIAAT; Vesicular Glutamate Transporter 1, 2, 3/VGLUT; Vesicular monoamine transporter 1, 2

Acetylcholine/ACh; Acetylcholinesterase; Choline Acetyltransferase/ChAT; Choline transporter; Vesicular Acetylcholine Transporter/VAChT

Adrenaline; Dopamine; Dopamine Beta Hydroxylase/DBH; Dopamine Transporter/DAT; L-DOPA; Nitric Oxide-Dopamine; Norepinephrine; Norepinephrine Transporter/NET; Parkin; Tyrosine Hydroxylase/TH; TorsinA


DL-5-Hydroxytryptophan; Serotonin; Serotonin Transporter/SERT; Tryptophan Hydroxylase

DARPP-32; GABA; GABA Transporters 1; GABA Transporters 2; GABA Transporters 3; Glutamate Decarboxylase/GAD; Vesicular GABA transporter/VGAT/VIAAT

Glutamate; Glutamate Transporter; Glutamine; Glutamine Synthetase; Vesicular Glutamate Transporter 1; Vesicular Glutamate Transporter 2; Vesicular Glutamate Transporter 3
"""

from Bio import Entrez
from probe_designer.mRNA_designer import RNARetriever2, design_step
from probe_designer.probe_refiner import ProbeFilter, probe_set_refiner
import json
retriever = RNARetriever2()
filterer = ProbeFilter(db='gencode_tracks_reversed_introns+mRNA', copy_num='brain')
Entrez.email = 'elubeck@caltech.edu'
passed = {}
try:
    with open("brain_type_genes.json".format(), 'r') as f_out:
        used = json.load(f_out).keys()
except FileNotFoundError:
    print("No json")
    used = set()
for frag in gene_str.split(";"):
    for sub_frag in frag.split('/'):
        name = sub_frag.strip(" \n")
        name1 = name[0].upper() + name[1:].lower()
        # Make best set of probes
        if name in used:
            continue
        name, probes, seq = design_step(name1, cds_only=True, length=26, spacing=0,
                                        gc_target=0.55, gc_min=0.35, gc_max=0.75)
        if name == "FAILED":
            esearch = Entrez.read(Entrez.esearch(db='gene',
                                                term='{}[gene] AND "Mus musculus"[orgn]'.format(
                                                    name1),
                                                retmode='xml'))
            esearch2 = Entrez.read(Entrez.esearch(db='gene',
                                                 term='{} AND "Mus musculus"[orgn]'.format(
                                                     name1),
                                                 retmode='xml'))
            if len(esearch['IdList']) == 1:
                result = Entrez.read(Entrez.efetch(db='gene',
                                                id=esearch['IdList'][0],
                                                retmode='xml'))
            elif len(esearch2['IdList']) == 1:
                result = Entrez.read(Entrez.efetch(db='gene',
                                                   id=esearch2['IdList'][0],
                                                   retmode='xml'))
            else:
                continue
            name1 = result[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_locus']
            name, probes, seq = design_step(name1, cds_only=True, length=26, spacing=0,
                                            gc_target=0.55, gc_min=0.35, gc_max=0.75)
        if name in used:
            continue
        probes = set(probes)
        if not probes:
            continue
        probes = filterer.run(probes, name, match_thresh=14, n_probes=48,
                              max_off_target=2000, off_target_hits=6)
        print(name, len(probes), len(''.join(seq)))
        passed[name] = probes
        with open("brain_type_genes.json".format(), 'w') as f_out:
            json.dump(passed, f_out)


import arrow
import csv
with open("brain_genes_{}.csv".format(arrow.now().timestamp), 'w') as f_out:
    cv = csv.writer(f_out)
    for name, genes in passed.items():
        for gene in genes:
            cv.writerow([name, gene])


print("DONE")

