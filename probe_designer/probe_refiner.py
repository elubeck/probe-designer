from builtins import list, range, len, set, enumerate
from collections import defaultdict, Counter
import copy
import csv
from itertools import groupby
import random
import blaster
from utils.misc import gc_count

__author__ = 'Eric'


def probe_set_refiner(pset_i, block_size=18):
    """
    Checks that probes don't hit the same target of block_size.  Drops probes if they do.  Returns a filtered probeset.
    Makes the assumption that all probes from different sets that hit the same target must be dropped.  This assumption can be faulty in the case of non-mRNA targeting sequences such as adapters added to probes.  Prefered behavior in this case would be to drop nothing or only mRNA binding sequences.
    """
    pset = copy.deepcopy(pset_i)
    p_lookup = defaultdict(list)
    flat_pfrags = []
    for gene, probes in pset.items():
        for probe in probes:
            for i in range(0, len(probe) - block_size):
                p_lookup[probe[i:i + block_size]].append(gene)
                flat_pfrags.append(probe[i:i + block_size])
    # Get probes that hit multiple targets
    counts = [(k, v) for k, v in Counter(flat_pfrags).items() if v > 1]
    counts = [(bad_probe, n_hits) for bad_probe, n_hits in counts
              if len(set(p_lookup[bad_probe])) > 1]
    for probe, n_hits in counts:
        for target in set(p_lookup[probe]):
            for n, probe_seq in enumerate(pset[target]):
                if probe in probe_seq:
                    del pset[target][n]
    return pset


class ProbeFilter(object):
    def run_batch(self):
        hit_val2 = sorted(hit_vals, key=lambda x: x[0].split(',')[0])
        finished_probes = {}
        # Iterate through every probeset
        for target, probe_set in groupby(hit_val2,
                                         key=lambda x: x[0].split(',')[0]):
            pass

    def filter_probes(self, probe_set, target, hits, probe_lookup,
                      max_hits=10000,
                      off_target_thresh=7,
                      probe_num=48, ):
        # Filter out probes with greater than max_hits off target hits
        probe_set = [probe for probe in probe_set if probe[1] < max_hits]
        # First pick probes with no off target hits
        passed_probes = [probe for probe in probe_set if probe[1] == 0]
        remaining_probes = set(probe_set)
        remaining_probes.difference_update(passed_probes)

        # Sort remaining probes by off target hits
        choices = sorted(list(remaining_probes), key=lambda x: x[1])
        choice_index = 0
        off_target = Counter([])

        # Iterate until at least probe_num probes are picked or no more probes are available
        while len(passed_probes) < probe_num and choice_index != len(
            remaining_probes):
            selected = choices[choice_index]
            # Get off target hits for chosen probe
            off_target_hit = [k for k in hits[selected[0]]
                              if target != k.split(',')[1]]
            test_counter = off_target + Counter(off_target_hit)

            # Check if adding this probe makes anything go over off target threshold
            over_thresh = [True for k in test_counter.values()
                           if k >= off_target_thresh]
            if not any(over_thresh):
                passed_probes.append(selected)
                off_target = test_counter
            choice_index += 1

        # If more than 24 probes chosen optimize on GC
        if len(passed_probes) > probe_num:
            # Get GC counts of every probe
            probe_gc = [(probe[0], gc_count(probe_lookup[probe[0]]))
                        for probe in passed_probes]
            gc_target = 0.55
            gc_range = 0.01
            multiplier = 1
            chosen_gc = []
            # Get closest probe set to 0.55 gc
            while len(chosen_gc) < probe_num:
                gc_min = gc_target - gc_range * multiplier
                gc_max = gc_target + gc_range * multiplier
                chosen_gc = [probe for probe, gc in probe_gc
                             if gc_max >= gc >= gc_min]
                multiplier += 1
            # If too many probes still exists choose a random subset
            if len(chosen_gc) != probe_num:
                chosen_gc = random.sample(chosen_gc, probe_num)
            passed_probes = chosen_gc
        else:
            passed_probes = [probe for probe, n in passed_probes]
        return passed_probes

    def blast2copynum(self, hits, drop_self=True):
        """
        drop_self = True: make sure probes that probes that hit themselves are not counted in blast tally
        """
        hit_vals = []
        for probe_name, matches in hits.items():
            gene_name = probe_name.split(",")[0].lower()
            if any(matches):  # Added incase no matches
                gencode_id, refseq = map(list, zip(*[match.lower().split(',')
                                                     for match in matches]))
            else:
                gencode_id, refseq = [], []
            if drop_self:
                if not gene_name in refseq:
                    continue
                    bad_count += 1
                    # raise Exception("Query not found in hits")
                del gencode_id[refseq.index(gene_name)]
            try:
                hit_vals.append((probe_name, self.get_copynum(gencode_id)))
            except:
                print("Failed @ {}".format(gencode_id))
        hit_vals = sorted(hit_vals, key=lambda x: x[1], reverse=True)
        return hit_vals

    def run_blast(self, probe_lookup, match_thresh,
                  strand=None,
                  db='gencode_tracks_reversed'):
        if not strand:
            strand = self.strand
        fasta_str = "\n".join(">{}\n{}".format(*items)
                              for items in probe_lookup.items())
        res = blaster.local_blast_query(fasta_str, db=db)
        hits = blaster.parse_hits(res,
                                   strand=strand,
                                   match_thresh=match_thresh)
        return hits

    def run(self, flat_probes, gene_name,
            match_thresh=18,
            n_probes=48,
            max_off_target=10000):
        """
        flat_probes is a list of probes.  Sequence only
        """
        probe_lookup = {
            "{},{}".format(gene_name, n): probe
            for n, probe in enumerate(flat_probes)
        }
        res = self.run_blast(probe_lookup, match_thresh, db=self.db)
        hit_vals = self.blast2copynum(res)
        filtered_probes = self.filter_probes(hit_vals, gene_name, res,
                                             probe_lookup,
                                             probe_num=n_probes,
                                             max_hits=max_off_target)
        finished_probes = [probe_lookup[probe] for probe in filtered_probes]
        return finished_probes

    def get_copynum(self, hits):
        split_hits = (hit.split(',')[0] for hit in hits)
        false_hits = sum(self.counts[name] for name in split_hits
                         if name in self.counts.keys())
        return false_hits

    def __init__(self,
                 db='gencode_tracks_reversed',
                 strand='+',
                 copy_num='embryonic11.5'):
        self.db = db
        if strand == '+':
            self.strand = 1
        else:
            self.strand = -1
        if copy_num == 'embryonic11.5':
            # Get merged embryonic 11.5 encode data
            with open('db/encode_counts.csv', 'r') as f:
                self.counts = {
                    line[0]: float(line[1])
                    for line in csv.reader(f)
                }
        elif copy_num == 'brain':
            with open('db/brain_counts.csv', 'r') as f:
                self.counts = {
                    line[0]: float(line[1])
                    for line in csv.reader(f)
                }
        else:
            raise Exception('need a valid copy_num database')

