import intron_designer2

import blaster2

filterer = intron_designer2.ProbeFilter()
###### Blast Adapters ###########
adapters = []
with open("db/70k_bridges.txt", "r") as f_in:
    for line in f_in:
        adapters.append(line.strip("\n"))
adapters = list(set(adapters))
print(len(adapters))

fasta_adapter = "\n".join(">{}\n{}".format(n, probe)
                          for n, probe in enumerate(adapters))
blast_hits_adapter = blaster2.local_blast_query(
    fasta_adapter,
    db='gencode_tracks_reversed_introns+mRNA',
    strand='plus')
blast_res_adapter = blaster2.parse_hits(blast_hits_adapter,
                                        match_thresh=10,
                                        strand=1)
best_oligos = [bridge for bridge, hits in blast_res_adapter.iteritems()
               if len(hits) == 0]
oligo_ind = sorted(best_oligos)

passed_adapters = [adapters[int(ind)] for ind in oligo_ind]

bad_seqs = [4 * 'a', 4 * 't', 4 * 'c', 4 * 'g']
passed_2 = [adapter for adapter in passed_adapters
            if not any(seq in adapter for seq in bad_seqs) if adapter]

with open("db/passed_bridges_mouse.txt", 'w') as f_out:
    for oligo in passed_2:
        f_out.write("{}\n".format(oligo))

# ###### Blast Adapters ###########
# yandong_adapters = []
# with open("temp/yandong_adapters.txt", "r") as f_in:
#     for line in f_in:
#         yandong_adapters.append(line.strip("\n"))
# adapters = list(set(yandong_adapters))

# fasta_adapter = "\n".join(">{}\n{}".format(n, probe)
#                           for n, probe in enumerate(adapters))
# blast_hits_adapter = blaster2.local_blast_query(
#     fasta_adapter,
#     db='gencode_tracks_reversed_introns+mRNA',
#     strand='plus')
# blast_res_adapter = blaster2.parse_hits(blast_hits_adapter,
#                                         match_thresh=10,
#                                         strand=1)

# not_pgk1 = [k for k, v in blast_res_adapter.items()
#             if not any(['Pgk1' in e for e in v])]
# off_target = {
#     bridge: filterer.get_copynum(hits)
#     for bridge, hits in blast_res_adapter.iteritems()
# }

# best_oligos = [bridge for bridge, hits in blast_res_adapter.iteritems()
#                if len(hits) == 0]
# oligo_ind = sorted(best_oligos)

# passed_adapters = [adapters[int(ind)] for ind in oligo_ind]

# bad_seqs = [4 * 'a', 4 * 't', 4 * 'c', 4 * 'g']
# passed_2 = [adapter for adapter in passed_adapters
#             if not any(seq in adapter for seq in bad_seqs) if adapter]

# with open("db/passed_bridges_mouse.txt", 'w') as f_out:
#     for oligo in passed_2:
#         f_out.write("{}\n".format(oligo))

# #### CHECK BLAST DB
# pgk1 = ['gcttttgaagcgtgcagaat',
#         'aacaccgtgaggtcgaaagg',
#         'tcagcttgttggaaagcgac',
#         'taggaacgttgaagtccacc',
#         'gtttgttatctggttgttct',
#         'ttggaacagcagccttgatc',
#         'ccattgtccaagcagaattt',
#         'ggctcataaggacaacggac',
#         'ctggctctaaggagtacttg',
#         'gcagagatttgagttcagca',
#         'cccacacaatccttcaagaa',
#         'acaggcattctcgacttctg',
#         'catgaaagcggaggttttcc',
#         'ctcagctttaaccttgtttc',
#         'tgaggctcggaaagcatcaa',
#         'agacatctcctagtttggac',
#         'gtcccaaaagcatcattgac',
#         'ccaaagccttggcaaagtag',
#         'caagatagccaggaagggtc',
#         'attgatcagctggatcttgt',
#         'ccttaaggaaggtaaaggcc',
#         'ccaatctccatgttgttgag',
#         'ctccttcttcatcatacaga',
#         'tttctcagctttggacatga', ]

# fasta_adapter = "\n".join(">{}\n{}".format(n, probe)
#                           for n, probe in enumerate(pgk1))
# bra = {}
# dbl = ['gencode_tracks_reversed', 'gencode_tracks_reversed_introns+mRNA']
# for db in dbl:
#     blast_hits_adapter = blaster2.local_blast_query(fasta_adapter,
#                                                     db=db,
#                                                     strand='plus')
#     bra[db] = blaster2.parse_hits(blast_hits_adapter, match_thresh=10, )
# for p_num in range(20):
#     d = [bra[db][str(p_num)] for db in dbl]

# # import numpy as np
# # from scipy.spatial.distance import pdist
# # char_map = {'a':0, 't':1, 'g': 2, 'c':3}
# # p2= np.array([ list(char_map[e] for e in v) for v in passed_2 ])
# # d = pdist(p2, metric='hamming') * 20
# # bob = np.min(d)
