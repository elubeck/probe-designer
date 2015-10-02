from __future__ import print_function, division, unicode_literals
import csv
from probe_designer.intron_designer import ProbeFilter

pf = ProbeFilter()
with open('/home/eric/20mer_adaptors.csv', 'rU') as fin:
    flat_probes = [line[0] for line in csv.reader(fin)]
with open('20mer_adapter2.txt', 'r') as fin:
    flat_probes2 = [line[0] for line in csv.reader(fin)]
flat_probes = set(flat_probes + flat_probes2 )

probe_lookup = {"{},{}".format(n*3, n*2): probe for n, probe in enumerate(flat_probes)}
res = pf.run_blast(probe_lookup, match_thresh=12)
hit_vals = pf.blast2copynum(res, drop_self=False)
hit_vals = sorted(hit_vals, key=lambda x: x[1])
hit_vals = [(probe_lookup[hit], n) for hit, n in hit_vals]
h2 = {k for k,c in hit_vals if c == 0 }

from probe_designer.get_seq import reverse_complement
probe_lookup = {"{},{}".format(n*3, n*2): reverse_complement(probe) for n, probe in enumerate(h2)}
res = pf.run_blast(probe_lookup, match_thresh=12)
hit_vals = pf.blast2copynum(res, drop_self=False)
hit_vals = list(set(hit_vals))
hit_vals = sorted(hit_vals, key=lambda x: x[1])
hit_vals2 = [(probe_lookup[hit], n) for hit, n in hit_vals]
# h3 = set()
# i = 0 
# i2 = []
# while len(h3) < 1000:
#     k, c = hit_vals[i]
#     h3.add(reverse_complement(k))
#     i2.append(reverse_complement(k))
#     i+=1
# with open('replacement_adapters.txt', 'w') as f_out:
#     f_out.write('\n'.join(i2[-90:]))
h3 = [reverse_complement(k) for k,c in hit_vals2 ]
        
with open('adapters.txt', 'w') as f_out:
    f_out.write('\n'.join(h3))

# with open('adapter_bridge_hits.txt', 'w') as f_out:
#     for line in hit_vals:
#         f_out.write(",".join(map(str, line))+'\n')


with open('/home/eric/old_adapters.txt', 'r') as fin:
    old_adapters = {line.strip('\n') for line in fin}

non_overlapping = [p for p in h3 if p not in old_adapters][:101]
with open('/home/eric/fixed_adapters.txt', 'w') as fout:
    fout.write('\n'.join(non_overlapping))
        