from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tempfile import NamedTemporaryFile
from probe_designer.probe_refiner import ProbeFilter


class DNARetriever(object):
    def get_single_dna(self, coord, chrom, strand, chunk_size=1000):
        """
        :param coord: (start, end)
        :param chunk_size:
        :return:
        """
        start_coord, end_coord = coords

        chrom_path = self.chromo_folder.joinpath(
            "{}.fa.masked".format(chrom))
        chrom_seq = SeqIO.read(str(chrom_path), 'fasta')
        return self.get_dna(start_coord, end_coord, strand, chrom_seq, chunk_size)

    def get_dna(self, start, end, strand, chrom_seq, chunk_size=1000):
        seq = chrom_seq[start:end]
        if strand == "+": # TODO: YODAI Should this be positive or negative?
            seq = seq.reverse_complement()
        # Drop masked sequences
        for sub_seq in seq.seq.split("N"):
            if len(sub_seq) >= 20:
                name = "start:{}, end:{}, strand:{}".format(start,end,strand)
                record = SeqRecord(sub_seq,
                                   name=name,
                                   id=name,
                                   description='')
                for start in range(0, len(record), chunk_size):
                    yield record[start:start + chunk_size]

    def __init__(self, organism='mouse', skip_probes=None,):
        self.organism = organism
        self.chromo_folder = Path("db/{}/chromFaMasked/".format(organism))

def design(genes, probe_length=50, chunksize=1):
    """

    :param genes: dictionary of {name: {start, end, strand, chrom}}
    :param probe_length: integer of probe length
    :param chunksize: integer amount of genes to run oligoarray at once.
    :return:
    """
    dna_retriever = DNARetriever()
    chunks = []
    results = []
    for name, g_vals in genes.items():
        # TODO: Yodai you really need to make a blast database for DNA
        blast_path = os.path.join(os.path.expanduser("~"),
                                  "blastdb", "old_format",
                                  "transcribed_mouse2")
        o = OligoarrayDesigner(
            blast_db=blast_path)
        chunksize = 12
        coords = (g_vals['start'], g_vals['end'])
        gene_chunks = list(dna_retriever.get_single_dna(coords, g_vals['chrom'], g_vals['strand']))
        chunks.append(gene_chunks)
        if len(chunks) == chunksize:
            with NamedTemporaryFile("w") as fasta_input:
                # Break records into 1000nt chunks
                for gene_chunks in chunks:
                    SeqIO.write(gene_chunks, fasta_input, 'fasta')  # write all sequences to file
                fasta_input.flush()
                res = o.run(fasta_input.name,
                            max_dist=1000,
                            min_length=probe_length,
                            max_length=probe_length,
                            max_tm=100,
                            cross_hyb_temp=72,
                            min_tm=74,
                            secondary_struct_temp=76,
                            max_oligos=100,
                            timeout=15)
                for gene in res:
                    results.append(gene)
    return results

initial_probes = design(genes)

probe_filterer = ProbeFilter(db='gencode_tracks_reversed_introns+mrna')
passed = []
for gene in initial_probes:
    probes_that_are_good_for_one_gene = probe_filterer.run(flat_probes, )
    passed.append(probes_that_are_good_for_one_gene)


import probe_refiner
p_set2 = probe_refiner.probe_set_refiner(all_p2)

p_set2 = {k: v for k, v in p_set2.items() if len(v) >= 24}

# blast everything to check that all probes together don't cause big problems
flat_probes = [
    {"gene": gene,
     "seq": probe,
     "number": n,
     "name": "{}-{}".format(gene, n)}
    for gene, probes in p_set2.items() for n, probe in enumerate(probes)
    ]

import blaster
fasta_str = "\n".join(">{name}\n{seq}".format(**probe)
                      for probe in flat_probes)
blast_hits2 = blaster.local_blast_query(
    fasta_str,
    db='gencode_tracks_reversed_introns+mrna')
blast_res2 = blaster.parse_hits(blast_hits2, match_thresh=18, strand=1)

from collections import defaultdict
# for every gene hit in blast results see what probes hit it
hit_to_probe = defaultdict(list)
for probe, hits in blast_res2.items():
    for hit in hits:
        hit_to_probe[hit].append(probe)

# get probes that hit unintended target
bad_hits = [(hit, [probe for probe in probes
                   if probe.split('-')[0].lower() not in hit.lower()])
            for hit, probes in hit_to_probe.items()]

bad_hits = sorted([(k, hits) for k, hits in bad_hits if any(hits)],
                  key=lambda x: len(x[1]),
                  reverse=true)

# drop probes if more than thresh probes hit x-target or copy number is greater than c_thresh
thresh = 5
c_thresh = 50
drop_list = []
for hit, probes in bad_hits:
    ensembl_name = hit.split(',')[0]
    if len(probes) >= thresh or probe_filterer.get_copynum(
            [ensembl_name]) > c_thresh:
        for probe in probes:
            drop_list.append(probe)

drop_probes = set(drop_list)

# make final list of good probes
flat_probes_filtered = [probe for probe in flat_probes
                        if probe['name'] not in drop_probes]

from itertools import groupby
# count remaining probes
grouped_probes = [(group, list(probes))
                  for group, probes in groupby(flat_probes_filtered,
                                               key=lambda x: x['gene'])]
passed = {g for g, probes in grouped_probes if len(probes) >= 24}
grouped_probes2 = [(g, p) for g, p in grouped_probes if g in passed]
p_set2 = {g: {p['seq'] for p in probes} for g, probes in grouped_probes2}

n_oligos = len([g for k, v in p_set2.items() for g in v])

with open("temp/intron_11-15.csv", "wb") as f_out:
    c = csv.writer(f_out)
    for gene, probes in p_set2.items():
        for n, probe in enumerate(probes):
            c.writerow([gene, n, probe])
