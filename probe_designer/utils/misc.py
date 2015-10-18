from __future__ import division
from builtins import len, sum

__author__ = 'Eric'


def gc_count(probe):
    """
    Returns the % GC for a string
    """
    return len([1 for c in probe.lower() if c in ['c', 'g']]) / len(probe)


def n_probes(chunk_list, probe_size=35):
    # Get # of probes that could be made from set of gene chunks
    return sum(len(chunk) // probe_size for chunk in chunk_list)


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement[c] for c in seq.upper())[::-1]


def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement[c] for c in seq.upper())
