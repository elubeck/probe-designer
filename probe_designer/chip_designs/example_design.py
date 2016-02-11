import os
import sys

from probe_designer.mRNA_designer import RNARetriever2, design_step
from probe_designer.probe_refiner import ProbeFilter

retriever = RNARetriever2()
filterer = ProbeFilter(db='gencode_tracks_reversed', copy_num='brain')

name, probes, sequence = design_step('Pgk1', cds_only=True, length=35)
probes = set(probes)
final_probes = filterer.run(probes, name, match_thresh=18, n_probes=48, max_off_target=2000)


print("BOB")