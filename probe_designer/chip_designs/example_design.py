import os
import sys

from probe_designer.mRNA_designer import RNARetriever2, design_step
from probe_designer.probe_refiner import ProbeFilter

retriever = RNARetriever2()
filterer = ProbeFilter(db='gencode_tracks_reversed_introns+mRNA', copy_num='brain')

name, probes, sequence = design_step('Pgk1', cds_only=True, length=35,
                                     spacing=0, gc_target=0.55, gc_min=None, gc_max=None)
probes = set(probes)
final_probes = filterer.run(probes, name, match_thresh=18, n_probes=48, max_off_target=2000)

print(len(final_probes))

