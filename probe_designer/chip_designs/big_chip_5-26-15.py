from __future__ import division, print_function, unicode_literals
import probe_designer.mRNA_designer
from progressbar import ProgressBar
import csv
import dataset
import random

####### Brain Genes
# # First Get Long's 100 Genes from List
# gene_list = []
# with open("/home/eric/brain100list2.csv", 'rU') as fin:
#     for line in csv.reader(fin):
#         gene_list.append(line)

# # Design All Genes
# flat_gene_list = [gene for gene, exp_level in gene_list]
# gene_lookup = mRNA_designer.batch_design(flat_gene_list)

# flat_gene2 = ['Lhx1','Lhx3','Lhx4','Lhx6','Lhx2','Neurod4','Neurog1','Foxa1','Foxa2','Foxc1','Foxd1','Foxd3','Foxd4','Foxg1','Foxh1','Foxi1','Foxm1','Foxn1','Foxo1','Gata6','Barhl1','Bcl3','Bhlhe41','Blzf1','Cdc5l','Cdc6','Cdx2','Hoxb3','Hoxb8','Creb1','Creb3l1','Crx','Dbx1','Ddb2','Nfe2l2','Nfkb2','Nfkbiz','Nfya','Nhlh1','Nkx6-1']

# gene_lookup2 =  mRNA_designer.batch_design(flat_gene2)
# # Filter genes on blast
# from intron_designer2 import ProbeFilter
# probe_filterer = ProbeFilter(copy_num='brain')
# db = dataset.connect("sqlite:///db/tf_mrna_probes.db")
# p_table = db['unfiltered_probes']
# finished_probes = db['filtered_probes']
# flat_probes = []
# p_bar = ProgressBar(maxval=len(gene_list)).start()

# debug = {gene: len([p['seq'] for p in p_table.find(target=gene)]) for gene in enumerate(gene_lookup2.keys())}
# debug2 = {k:v for k, v in debug.iteritems() if v>=24}
# for n, gene in enumerate(gene_lookup.keys()):
#     p_bar.update(n)
#     hits = [{'target':probe['target'], 'seq':probe['seq']}
#             for probe in finished_probes.find(target=gene)]
#     if hits:
#         flat_probes.append(hits)
#     else:
#         p_set = [p['seq'] for p in p_table.find(target=gene)]
#         if len(p_set) < 24: continue
#         iters = 0 
#         f_probes = []
#         while iters < 5:
#             max_off = 1 + iters * 2000
#             f_probes = probe_filterer.run(p_set, gene, n_probes=30, max_off_target=max_off)
#             if len(f_probes) >=30: break
#             iters += 1
#         f_probes = [{'target':gene, 'seq':probe} for probe in f_probes]
#         finished_probes.insert_many(f_probes)
#         flat_probes.append(f_probes)
# p_bar.finish()

# flat_probes2 = []
# for n, gene in enumerate(gene_lookup2.keys()):
#     p_bar.update(n)
#     hits = [{'target':probe['target'], 'seq':probe['seq']}
#             for probe in finished_probes.find(target=gene)]
#     if hits:
#         flat_probes2.append(hits)
#     else:
#         p_set = [p['seq'] for p in p_table.find(target=gene)]
#         if len(p_set) < 24: continue
#         iters = 0 
#         f_probes = []
#         while iters < 5:
#             max_off = 1 + iters * 2000
#             f_probes = probe_filterer.run(p_set, gene, n_probes=30, max_off_target=max_off)
#             if len(f_probes) >=30: break
#             iters += 1
#         f_probes = [{'target':gene, 'seq':probe} for probe in f_probes]
#         finished_probes.insert_many(f_probes)
#         flat_probes2.append(f_probes)
# p_bar.finish()
# high={'Acap1','Vangl','Ddx60','Fancc','Ak4','Tdo2','Fat4','Klk8','Nr3c'}

# flat2 = {gene[0]['target']:[p['seq'] for p in gene] for gene in flat_probes2 if gene if len(gene) >= 24}
# flat = {gene[0]['target']:[p['seq'] for p in gene] for gene in flat_probes if gene if len(gene) >= 24}
# flat = {k: v for k,v in flat.iteritems() if k not in high}
# total_probes = 100
# remaining_choices = total_probes - len(flat)

# while remaining_choices != 0:
#     flat2_choices = {k: flat2[k] for k in random.sample(flat2, remaining_choices)}
#     flat.update(flat2_choices)
#     refined = mRNA_designer.probe_set_refiner(flat, block_size=13)
#     refined2 = {}
#     for k,v in refined.iteritems():
#         if len(v) < 24:
#             pass
#         elif len(v)<= 30:
#             refined2[k] = v
#         else:
#             refined2[k] = random.sample(v, 30)
#     remaining_choices = total_probes - len(refined2)

# refined = {k: refined2[k] for k in random.sample(refined2, 100) if len(refined2[k])>=24}
# refined = [{'target': gene, 'seq': p} for gene, seqs in refined.iteritems() if len(seqs) >= 24 for p in seqs]

# with open('5-26-15_brain.csv', 'w') as f_out:
#     fwriter = csv.writer(f_out)
#     fwriter.writerows(seq.values() for seq in refined)
#     # fwriter.writerows(seq.values() for seq in random.sample(list(groupby(refined, lambda x: x['target'])), 100))

############### INTRONS ###################
import probe_designer.probe_refiner
from probe_designer.utils.misc import reverse_complement
import dataset
import probe_designer.utils.misc

with open("/home/eric/tflistall.txt", "r") as fin:
    tfs = [line.strip('\r\n').lower() for line in fin]

with open('/home/eric/signalinglist.txt', 'r') as fin:
    signaling = [line.strip('\r\n').lower() for line in fin]

intron_db = dataset.connect("sqlite:///db/intron_probes2.db")
intron_db_filtered = dataset.connect("sqlite:///db/intron_probes_filtered.db")
probe_db = intron_db['mouse']
filtered_probe_table = intron_db_filtered['mouse']
all_targets = set(tfs + signaling)

n_probes = 0
all_probes = {}
for gene_d in filtered_probe_table.distinct('target'):
    gene = gene_d['target']
    probes = [p['seq'] for p in filtered_probe_table.find(target=gene)]
    if len(probes) >= 48 and gene.lower() in all_targets:
        all_probes[gene] = probes

# Introns - I was placing this in the intron_designer2 terminal
adapter_seq = []
with open('adapters.txt', 'r') as fin:
    for line in fin:
        seq = line.strip('\n')
        r_probe = reverse_complement(seq)
        adapter_seq.append(r_probe)

r_adapter = list(reversed(adapter_seq))

import probe_designer.intron_designer as id2
import random
from progressbar import ProgressBar
import probe_designer.mRNA_designer

p_set2 = probe_designer.probe_refiner.probe_set_refiner(all_probes)

# Rock1 is short
p_keys = set(p_set2.keys()) - {'Rock1'}
f_probes = {
    k: ['TATA' + seq for seq in p_set2[k]]
    for n, k in enumerate(random.sample(p_keys, 1000))
}

probe_filterer_intron = probe_designer.probe_refiner.ProbeFilter()
intron_probes2 = []
p_bar = ProgressBar(maxval=1000).start()
# # used_probes = {s['target'] for probe in intron_probes for s in probe}
for pn, (gene, probes) in enumerate(f_probes.iteritems()):
    # for pn, (gene, probes) in enumerate(p_set.iteritems()):
    p_bar.update(pn)
    # if gene in used_probes: continue
    iters = 0
    f_probes = []
    while iters < 5:
        max_off = 1 + iters * 2000
        f_probes = probe_filterer_intron.run(probes, gene,
                                             n_probes=48,
                                             max_off_target=max_off)
        if len(f_probes) >= 48: break
        iters += 1
    f_probes = [{'target': gene, 'seq': probe} for probe in f_probes]
    intron_probes2.append(f_probes)
    print(pn, len(f_probes), )
p_bar.finish()

import csv
flat_probes = [p for gene in intron_probes2 for p in gene]

with open('intron_probes_5-27-15.csv', 'w') as f_out:
    cw = csv.writer(f_out)
    cw.writerows(p.values() for p in flat_probes)

primers = [['ACAGTGGAACGCGCCATGAG', 'CGCTTGTTGGCTGGTATGCG'], [
    'AATTGAGCAGCTCGGGCCAC', 'AGTTGCAGGCTTCCATCGCC'
], ['GACGCACATATGCGGGCAAG', 'TCCGCAGTCACGAAGATGCC'], ['ATTGAGGGTCTTCGCGTGCC',
                                                      'GGTTGCAAAGCGCCGGTTAC'],
           ['GATGGACTGGCAATTCCGGG',
            'GAATATGCCCGGCGAGAACG'], ['TGTGCGCTCCGATTGTCCTC',
                                      'GGCCAACAGACCCCATTTGC'],
           ['CCGCACGCCGTCCTTAAATC', 'AGATCCGGCAGCACGGAAAG'], [
               'TGCAGCTCCGCGAAATGAAG', 'AATGGCACAGACAGGCAGCG'
           ], ['CATGCCAGATCCCAATCCCC',
               'CAATCTGCGGGCATGTTTCG'], ['TCCATGAGCTCTGTTGCGGG',
                                         'CCGGCCTGAAAAACGTAGCG']]

p2 = [(forward, probe_designer.utils.misc.reverse_complement(reverse))
      for forward, reverse in primers]
random.shuffle(intron_probes2)
final_intron = []
for n, start in enumerate(range(0, len(intron_probes2), 100)):
    fp, rp = p2[n]
    for gene in intron_probes2[start:start + 100]:
        for probe in gene:
            ScaI = 'AGTACT'
            Spacer1 = 'TAG'
            Spacer2 = 'GAT'
            EcoRI = 'GAATTC'
            schematic = " ".join([fp, ScaI, Spacer1, probe['seq'], Spacer2,
                                  EcoRI, rp])
            final_intron.append({'target': probe['target'], 'seq': schematic})

with open('intron_probes+primers_5-27-15.csv', 'w') as f_out:
    cw = csv.writer(f_out)
    cw.writerows(p.values() for p in final_intron)

    # passed_probes = [pset for pset in intron_probes if len(pset)==48]
    # passed_adapters = {p['seq'][:20] for pset in passed_probes for p in pset}

    # failed_probes = [pset for pset in intron_probes if len(pset) != 48]
    # failed_adapters = {p['seq'][:20] for pset in failed_probes for p in pset}

    # used_adapters = passed_adapters.union(failed_adapters)
    # unused_adapters = [(a, c_num) for a, c_num in adapter_seq if a not in used_adapters]

    # import intron_designer2
    # pf_i = intron_designer2.ProbeFilter()
    # passed2 = []
    # for n2, probe in enumerate(failed_probes):
    #     target = probe[0]['target']
    #     iters = 0 
    #     f_probes = []
    #     chosen = unused_adapters[n2][0]
    #     probes = [chosen + p for p in p_set[target]]
    #     while iters < 5:
    #         max_off = 1 + iters * 2000
    #         f_probes = pf_i.run(probes, gene, n_probes=48, max_off_target=max_off)
    #         if len(f_probes) >=48: break
    #         iters += 1
    #     f_probes = [{'target':gene, 'seq':probe} for probe in f_probes]
    #     print(target, len(f_probes), len(used_adapters), chosen)
    #     passed2.append(f_probes)

    # passed2 = []
    # c_ind = 0
    # for n2, probe in enumerate(reversed(failed_probes)):
    #     target = probe[0]['target']
    #     while True:
    #         iters = 0 
    #         f_probes = []
    #         unused_adapters = [(a, c_num) for a, c_num in adapter_seq if a not in used_adapters]
    #         chosen = unused_adapters[0][0]
    #         probes = [chosen + p for p in p_set[target]]
    #         while iters < 5:
    #             max_off = 1 + iters * 2000
    #             f_probes = probe_filterer_intron.run(probes, gene, n_probes=48, max_off_target=max_off)
    #             if len(f_probes) >=48: break
    #             iters += 1
    #         f_probes = [{'target':gene, 'seq':probe} for probe in f_probes]
    #         print(n2, len(f_probes), c_ind, len(used_adapters))
    #         if len(f_probes) == 48:
    #             break
    #         else:
    #             used_adapters = used_adapters.union(chosen)
    #             c_ind + 1
    #     c_ind += 1
    #     passed2.append(f_probes)

    #     # files = ['tf-midbrainlist317.txt', 'candidategenes-all4000.txt', 'candidategenes1000-midbrain.txt', 'candidategenes173.txt']
    #     # file_path = ["/home/eric/Downloads/{}".format(f) for f in files]
    #     # gene_set = []
    #     # for file_name in file_path:
    #     #     with open(file_name, 'r') as f_in:
    #     #         for line in f_in:
    #     #             for gene in line.strip('\r\n').split('\t'):
    #     #                 if gene:
    #     #                     gene_set.append(gene)
    #     # gene_set = set(gene_set)
    #     # probes = batch_design(gene_set)
    #     # from intron_designer2 import ProbeFilter
    #     # from progressbar import ProgressBar
    #     # probe_filterer = ProbeFilter(copy_num='brain')
    #     # db = dataset.connect("sqlite:///db/tf_mrna_probes.db")
    #     # p_table = db['unfiltered_probes']
    #     # filtered_probes = db['filtered_probes']
    #     # used_probes = [target['target'] for target in filtered_probes.distinct('target')]
    #     # n_probes = len(list(p_table.distinct('target')))
    #     # p_bar = ProgressBar(maxval=n_probes).start()
    #     # to_use = [target['target'] for target in p_table.distinct('target') if target['target'] not in used_probes]
    #     # flat_use = ([p['seq'] for p in p_table.find(target=gene)] for gene in to_use)

    #     # def do_it(gene, p_set):
    #     #     if len(p_set) < 20: return
    #     #     f_probes= ProbeFilter(copy_num='brain').run(p_set, gene, max_off_target=5000)
    #     #     f_probes = [{'target':gene, 'seq':probe} for probe in f_probes]
    #     #     return f_probes

    #     # from joblib import Parallel, delayed
    #     # from itertools import izip
    #     # p = Parallel(n_jobs=-1, pre_dispatch=50)
    #     # res_dump = p(delayed(do_it)(gene, p_set) for gene, p_set in izip(to_use, flat_use))

    #     # for n, gene in enumerate(reversed(list(p_table.distinct('target')))):
    #     #     p_bar.update(n)
    #     #     gene = gene['target']
    #     #     if gene in used_probes:
    #     #         print("Skipping {}".format(gene))
    #     #         continue
    #     #     p_set = [p['seq'] for p in p_table.find(target=gene)]
    #     #     if len(p_set) < 20: continue
    #     #     f_probes= probe_filterer.run(p_set, gene, max_off_target=5000)
    #     #     f_probes = [{'target':gene, 'seq':probe} for probe in f_probes]
    #     #     if f_probes:
    #     #         tries = 0
    #     #         while tries != 10:
    #     #             try:
    #     #                 filtered_probes.insert_many(f_probes)
    #     #                 break
    #     #             except:
    #     #                 tries += 1
    #     # p_bar.finish()

    #     # n_probes = 0
    #     # with open("/home/eric/easy_passed.csv", 'w') as f_out:
    #     #     for gene_d in filtered_probes.distinct('target'):
    #     #         gene = gene_d['target']
    #     #         probes = [p['seq'] for p in filtered_probes.find(target=gene)]
    #     #         if len(probes) > 24:
    #     #             f_out.write("{},{}\n".format(gene, len(probes)))
    #     #             n_probes += 1

    #     # probes = [probe + 'CCAGC' for probe in gene_probes[gene]]
    #     # if len(probes) < 20:
    #     #     continue
    #     # passed_set[gene] = probe_filterer.run(probes, gene)

    # # import csv
    # # gene_list = []
    # # with open("/home/eric/brain100list2.csv", 'rU') as fin:
    # #     for line in csv.reader(fin,) :
    # #         gene_list.append(line)

    # # from intron_designer2 import ProbeFilter
    # # from progressbar import ProgressBar
    # # probe_filterer = ProbeFilter(copy_num='brain')
    # # db = dataset.connect("sqlite:///db/tf_mrna_probes.db")
    # # p_table = db['unfiltered_probes']
    # # flat_probes = []
    # # p_bar = ProgressBar(maxval=len(gene_list)).start()
    # # for n, (gene, exp_lvl) in enumerate(gene_list):
    # #     p_bar.update(n)
    # #     p_set = [p['seq'] for p in p_table.find(target=gene)]
    # #     if len(p_set) < 20: continue
    # #     iters = 0 
    # #     f_probes = []
    # #     while iters < 5:
    # #         max_off = 1 + iters * 2000
    # #         f_probes = probe_filterer.run(p_set, gene, n_probes=30, max_off_target=max_off)
    # #         if len(f_probes) >=30: break
    # #         iters += 1
    # #     f_probes = [{'target':gene, 'seq':probe} for probe in f_probes]
    # #     flat_probes.append(f_probes)
    # # p_bar.finish()
    # # used_probes = {k[0]['target'] for k in flat_probes if len(k) != 0}

    # # with open('/home/eric/brain_redundant_probes_filtered.json', 'r') as f_in:
    # #     passed_set =  json.load(f_in)
    # #     dog_cat = probe_set_refiner(passed_set, 18)
    # #     count = 0 
    # #     with open('/home/eric/brain_greater_than_20.csv', 'w') as f_ount:
    # #         for k, v in dog_cat.iteritems():
    # #             if len(v) >= 20:
    # #                 f_ount.write("{},{}\n".format(k, len(v)))
    # #     # with open('/home/eric/brain_redundant_probes_filtered_redundant.json', 'w') as f_out:
    # #     #     json.dump(dog_cat, f_out)

    # # with open('/home/eric/5-20-15_brain4genes.csv', 'w') as f_out:
    # #     for gene in ['Npy2r', 'Rgs14', 'Arhgef26', 'Gdf5']:
    # #         probes = dog_cat[gene]
    # #         if len(probes) > 24:
    # #             # Get GC counts of every probe
    # #             probe_gc = [(probe, gc_count(probe))
    # #                         for probe in probes]
    # #             gc_target = 0.55
    # #             gc_range = 0.01
    # #             multiplier = 1
    # #             chosen_gc = []
    # #             # Get closest probe set to 0.55 gc
    # #             while len(chosen_gc) < 24:
    # #                 gc_min = gc_target - gc_range * multiplier
    # #                 gc_max = gc_target + gc_range * multiplier
    # #                 chosen_gc = [probe for probe, gc in probe_gc
    # #                                 if gc_max >= gc >= gc_min]
    # #                 if len(chosen_gc) >= 24:
    # #                     chosen_gc = random.sample(chosen_gc, 24)
    # #                     break
    # #                 multiplier += 1
    # #             print(len(chosen_gc))
    # #             probes = chosen_gc
    # #         for n, line in enumerate(probes):
    # #             f_out.write("{},{},{}\n".format(gene, n+1, line))
