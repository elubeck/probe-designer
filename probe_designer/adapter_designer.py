from __future__ import print_function, division, unicode_literals
import csv
from intron_designer2 import ProbeFilter

pf = ProbeFilter()
with open('/home/eric/20mer_adaptors.csv', 'rU') as fin:
    flat_probes = [line[0] for line in csv.reader(fin)]
with open('20mer_adapter2.txt', 'r') as fin:
    flat_probes2 = [line[0] for line in csv.reader(fin)]
flat_probes = flat_probes + flat_probes2        

probe_lookup = {"{},{}".format(n*3, n*2): probe for n, probe in enumerate(flat_probes)}
res = pf.run_blast(probe_lookup, match_thresh=12)
hit_vals = pf.blast2copynum(res, drop_self=False)
hit_vals = sorted(hit_vals, key=lambda x: x[1])
hit_vals = [(probe_lookup[hit], n) for hit, n in hit_vals]
h2 = [k for k,c in hit_vals if c == 0 ]

from get_seq import reverse_complement
probe_lookup = {"{},{}".format(n*3, n*2): reverse_complement(probe) for n, probe in enumerate(h2)}
res = pf.run_blast(probe_lookup, match_thresh=12)
hit_vals = pf.blast2copynum(res, drop_self=False)
hit_vals = sorted(hit_vals, key=lambda x: x[1])
hit_vals = [(probe_lookup[hit], n) for hit, n in hit_vals]
h3 = [reverse_complement(k) for k,c in hit_vals if c == 0 ][:1000]
with open('adapters.txt', 'w') as f_out:
    f_out.write('\n'.join(h3))
    



# with open('adapter_bridge_hits.txt', 'w') as f_out:
#     for line in hit_vals:
#         f_out.write(",".join(map(str, line))+'\n')

