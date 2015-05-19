from __future__ import division, print_function, unicode_literals

import csv
import gzip
import random
import sys
from collections import Counter, defaultdict
from math import ceil
from tempfile import NamedTemporaryFile

import dataset
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path

import blaster2
import get_seq
from oligoarray_designer import Oligoarray, OligoarrayDesigner

csv.field_size_limit(sys.maxsize)  # Prevent field size overflow.


def gc_count(probe):
    return len([1 for c in probe.lower() if c in ['c', 'g']]) / len(probe)




# TODO: FINISH ME
class mRNARetriever(object):
    def get_mRNA(self, target, chunk_size=1000):
        hits = []
        for file_name in Path("db/").glob("mouse.*.rna.fna"):
            with open(str(file_name), 'r') as fin:
                for record in SeqIO.parse(fin, 'fasta'):
                    if "({})".format(
                        target.lower()) in record.description.lower():
                        hits.append(record)
        for sub_seq in get_seq.CDS(target).align_seqs(hits):
            if len(sub_seq) >= 20:
                record = SeqRecord(Seq(sub_seq),
                                   name=target,
                                   id=target,
                                   description='')
                for start in range(0, len(record), chunk_size):
                    yield record[start:start + chunk_size]

    def __iter__(self):
        with open("db/Refseq2Gene.tsv", "r") as f_in:
            ref_names = {
                gene_name
                for rnum, gene_name in csv.reader(f_in,
                                                  delimiter='\t')
            }
            for gene in ref_names:
                for record in self.get_mRNA(gene):
                    yield record

    def __init__(self, organism='mouse'):
        self.organism = organism


def n_probes(chunk_list, probe_size=35):
    # Get # of probes that could be made from set of gene chunks
    return sum(len(chunk) // probe_size for chunk in chunk_list)


class ProbeDesigner(object):
    def run_oligoarray(self, chunks):
        with NamedTemporaryFile("w") as fasta_input:
            # Break records into 1000nt chunks
            for gene_chunks in chunks:
                if n_probes(gene_chunks, self.probe_size) > 150:
                    # Reduce potential probeset size to drop design time
                    while n_probes(gene_chunks, self.probe_size) > 150:
                        gene_chunks = random.sample(gene_chunks,
                                                    len(gene_chunks) - 1)
                SeqIO.write(gene_chunks, fasta_input, 'fasta')
            fasta_input.flush()
            res = self.o.run(fasta_input.name,
                             max_dist=self.max_dist,
                             min_length=self.probe_size,
                             max_length=self.probe_size,
                             max_tm=self.max_tm,
                             cross_hyb_temp=self.cross_hyb_temp,
                             min_tm=self.min_tm,
                             secondary_struct_temp=self.secondary_struct_temp,
                             max_oligos=100,
                             timeout=15)
        for name, probes in res:
            probes['Percent GC'] = 100 * probes["Probe (5'-> 3')"].map(gc_count)
            yield probes

    def design(self, target):
        gene = [list(mRNARetriever().get_mRNA(target))]
        for probe in self.run_oligoarray(gene):
            print("TODO: ADD to DB")

    def batch_design(self, target_list, chunksize=3):
        from progressbar import ProgressBar
        p_bar = ProgressBar(maxval=len(target_list)).start()
        chunks = []
        n = 0
        probe_list = []
        while True:
            if len(chunks) == chunksize or n == len(target_list):
                for name, probes in self.run_oligoarray(chunks):
                    for probe in probes.T.to_dict().values():
                        probe_list.append(probe)
                    # TODO: ADD to DB
                chunks = []
                p_bar.update(n)
                if n == len(target_list):
                    break
            gene_name = target_list[n]
            chunks.append(list(mRNARetriever().get_mRNA(gene_name)))
            n += 1
        p_bar.finish()
        return probe_list

    def __init__(self,
                 db='/home/eric/blastdb/old_format/transcribed_mouse',
                 max_dist=1000,
                 probe_size=35,
                 max_tm=100,
                 cross_hyb_temp=72,
                 min_tm=74,
                 secondary_struct_temp=76, ):
        self.o = OligoarrayDesigner(blast_db=db)
        self.probe_size = probe_size
        self.max_dist = max_dist
        self.cross_hyb_temp = cross_hyb_temp
        self.min_tm = min_tm
        self.max_tm = max_tm
        self.secondary_struct_temp = secondary_struct_temp


def filter_probes_by_blast(probes):

    def probe_to_fasta(probe):
        return ">{Name},{Probe #}\n{Probe (5'-> 3')}\n".format(**probe)

    iv = [probe
          for probe_set in probes for probe in probe_set.T.to_dict().values()]
    probe_lookup = {
        "{Name},{Probe #}".format(**probe): probe["Probe (5'-> 3')"]
        for probe in probes
    }
    fasta_str = "".join(probe_to_fasta(probe) for probe in probes)
    # for n, probe in enumerate(probes.values):
    #     probe_lookup.append(("{},{}".format(gene, n), probe))
    # probe_lookup = dict(probe_lookup )
    # fasta_str = "\n".join(">{}\n{}".format(k, v) for k,v in probe_lookup.iteritems())

    strand = 1
    res = blaster2.local_blast_query(fasta_str, db='gencode_tracks_reversed')
    hits = blaster2.parse_hits(res, strand=strand, match_thresh=18)

    import pdb
    pdb.set_trace()
    hit_vals = []
    for probe_name, matches in hits.iteritems():
        gene_name = probe_name.split(",")[0]
        gencode_id, refseq = map(list, zip(*[match.split(',')
                                             for match in matches]))
        if not gene_name in refseq:
            continue
            bad_count += 1
            # raise Exception("Query not found in hits")
        del gencode_id[refseq.index(gene_name)]
        try:
            hit_vals.append((probe_name, blaster2.get_copynum(gencode_id)))
        except:
            print("Failed @ {}".format(gencode_id))
    hit_vals = sorted(hit_vals, key=lambda x: x[1], reverse=True)

    def gc_count(probe):
        return len([1 for c in probe.lower() if c in ['c', 'g']]) / len(probe)

    # Iterate everyg ene
    from itertools import groupby
    from collections import Counter
    hit_val2 = sorted(hit_vals, key=lambda x: x[0].split(',')[0])
    off_target_thresh = 7
    finished_probes = {}
    # Iterate through every probeset
    for target, probe_set in groupby(hit_val2,
                                     key=lambda x: x[0].split(',')[0]):
        probe_set = list(probe_set)

        # First pick probes with no off target hits
        passed_probes = [probe for probe in probe_set if probe[1] == 0]
        remaining_probes = set(probe_set)
        remaining_probes.difference_update(passed_probes)

        # Sort remaining probes by off target hits
        choices = sorted(list(remaining_probes), key=lambda x: x[1])
        choice_index = 0
        off_target = Counter([])

        # Iterate until at least 24 probes are picked or no more probes are available
        while len(passed_probes) < 24 and choice_index != len(remaining_probes):
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
        if len(passed_probes) > 24:
            # Get GC counts of every probe
            probe_gc = [(probe[0], gc_count(probe_lookup[probe[0]]))
                        for probe in passed_probes]
            gc_target = 0.55
            gc_range = 0.01
            multiplier = 1
            chosen_gc = []
            # Get closest probe set to 0.55 gc
            while len(chosen_gc) < 24:
                gc_min = gc_target - gc_range * multiplier
                gc_max = gc_target + gc_range * multiplier
                chosen_gc = [probe for probe, gc in probe_gc
                             if gc_max >= gc >= gc_min]
                multiplier += 1

            # If too many probes still exists choose a random subset
            if len(chosen_gc) != 24:
                chosen_gc = random.sample(chosen_gc, 24)
            passed_probes = chosen_gc
        else:
            passed_probes = [probe for probe, n in passed_probes]
        finished_probes[target] = [probe_lookup[probe]
                                   for probe in passed_probes]
        return finished_probes



def sub_seq_splitter(seq, size, gc_target=0.55, stride=1, spacing=2, gc_min=None, gc_max=None, o_mode=1):
    assert(isinstance(seq, unicode))
    # Generate all probes
    from itertools import groupby
    probes = []
    for seed in range(0, len(seq)-size, stride):
        p_seq = seq[seed:seed+size]
        p_gc = gc_count(p_seq)
        if gc_min and p_gc < gc_min:
            continue
        if gc_max and p_gc > gc_max:
            continue
        probes.append({'start': seed-spacing, 'end': seed+size+spacing,
                       'gc': p_gc, 'seq': p_seq})
    probes = sorted(probes, key=lambda x: abs(gc_target - x['gc']))

    # Build list of overlapping probes
    probe_indices = [set(range(probe['start'], probe['end'])) for probe in probes]
    if o_mode == 1:
        overlapping = []
        for probe_n, probe_ind in enumerate(probe_indices):
            for i in range(len(probe_indices)):
                if i != probe_n:
                    if any(probe_ind & probe_indices[i]):
                        overlapping.append((probes[probe_n], probes[i]))
        overlapping_2 = {probe['seq']: [hit for query, hit in over_probes]
                            for probe, over_probes in groupby(overlapping, key=lambda x: x[0])}
    if o_mode == 2:
        seq_lookup = {probe['seq']: probe for probe in probes}
        from collections import defaultdict
        overlapping = defaultdict(list)
        for probe_n, probe_ind in enumerate(probe_indices):
            for i in range(probe_n, len(probe_indices)):
                if i != probe_n:
                    if any(probe_ind & probe_indices[i]):
                        overlapping[probes[probe_n]['seq']].append(probes[i])
        import copy
        overlapping_b = copy.deepcopy(overlapping)
        for k, v in overlapping.iteritems():
            pl = seq_lookup[k]
            for i2 in v:
                overlapping_b[i2['seq']].append(pl)
        # Get probes that don't overlap anything!
        for seq in set(seq_lookup.keys()).difference(overlapping_b.keys()):
            overlapping_b[seq] = []
        overlapping_2 = {seq: [hit for query, hit in over_probes][0]
                         for seq, over_probes in groupby(overlapping_b.items(), key=lambda x: x[0])}
            
    # Group probes by difference from target gc
    probe_groups = [(k, list(v)) for k,v in groupby(probes, key= lambda x:abs(gc_target - x['gc']))]
    probe_index = 0
    passed_probes = []
    disallowed_probes = []

    # Iterate by probe GC and choose non-overlapping probes
    while probe_index != len(probe_groups):
        gc, probe_group = probe_groups[probe_index]
        # Filter already overlapping probes and convert list of dict to list of seqs
        probe_group = [probe['seq'] for probe in probe_group
                       if probe['seq'] not in disallowed_probes]
        # Iterate until no probes overlap in probe_group
        if len(probe_group) > 1:
            n_iters = 0
            while True:
                om = [(probe_seq, [overlapped for overlapped in overlapping_2[probe_seq] if overlapped['seq'] in probe_group])
                            for probe_seq in probe_group]
                om = sorted(om, key=lambda x: len(x[1]), reverse=True)
                # If worst probe doesn't overlap any other probes
                if len(om[0][1]) == 0:
                    break
                # Otherwise throw out worst probe
                else:
                    probe_group = [probe_seq for probe_seq, overind in om[1:]]
                n_iters += 1
            probe_group = [probe_seq for probe_seq, overind in om[1:]]
        for probe_seq in probe_group:
            passed_probes.append(probe_seq)
            for overlapped in overlapping_2[probe_seq]:
                disallowed_probes.append(overlapped['seq'])
        probe_index += 1

    # DEBUG_START
    from collections import defaultdict
    passed_or = defaultdict(list)
    for probe in passed_probes:
        for probe_m in probes:
            if probe_m['seq'] == probe:
                passed_or[probe].append(probe_m)
    dog_head = sorted([i['start'] for p in passed_or.values() for i in p])
    dog_poop = [i2-i1 for i2, i1 in zip(dog_head[1:], dog_head[:-1])]
    for v in dog_poop:
        if v <37:
            raise Exception("Probe design fail")
    # DEBUG_END
    return passed_probes

if __name__ == '__main__':
    # genes = []
    # with open("/home/eric/Downloads/candidategenes173.txt", "r") as f_in:
    #     for line in f_in:
    #         line = line.strip('\r\n')
    #         if line:
    #             genes.append(line.split('\t'))
    # flat_genes = [gene for line in genes for gene in line]
    # mr = mRNARetriever()
    # cds_records = {}
    # from collections import defaultdict
    # gene_probes = defaultdict(list)
    # from progressbar import ProgressBar
    # pbar = ProgressBar(maxval=len(flat_genes)).start()
    # failed = []
    # for n, gene in enumerate(flat_genes):
    #     try:
    #         gene_records = [str(record.seq.reverse_complement())
    #                             for record in mr.get_mRNA(gene, chunk_size=100000)]
    #     except:
    #         failed.append(gene)
    #     cds_records[gene] = gene_records
    #     for chunk in cds_records[gene]:
    #         for probe in sub_seq_splitter(chunk, 35, gc_min=0.35, o_mode=2):
    #             gene_probes[gene].append(probe)
    #     pbar.update(n)
    # pbar.finish()
    # import json
    # with open('brain_redundant.json', 'w') as f_out:
    #     json.dump(cds_records, f_out)
    # with open('brain_redundant_probes.json', 'w') as f_out:
    #     json.dump(gene_probes, f_out)
    # chunks = list(mRNARetriever().get_mRNA(flat_genes[0], chunk_size=10000))
    import json
    gene_probes = defaultdict(list)
    from progressbar import ProgressBar
    pbar = ProgressBar(maxval=173).start()
    n = 0
    with open('/home/eric/brain_redundant.json', 'r') as f_in:
        for gene, probes in json.load(f_in).iteritems():
            for chunk in probes:
                for probe in sub_seq_splitter(unicode(chunk), 35, gc_min=0.35, o_mode=2):
                    gene_probes[gene].append(probe)
            n += 1
            pbar.update(n)
        pbar.finish()
    with open('/home/eric/brain_redundant_probes.json', 'w') as f_out:
        json.dump(gene_probes, f_out)