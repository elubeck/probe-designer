import csv
from collections import defaultdict
import os

from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline as blastn
from tempfile import NamedTemporaryFile

os.environ['BLASTDB'] = os.path.join(os.path.expanduser("~"), "blastdb")

os.environ['$BLASTDB'] = os.path.join(os.path.expanduser("~"), "blastdb")


# Get merged embryonic 11.5 encode data
with open('db/encode_counts.csv', 'r') as f:
    counts = {line[0]: float(line[1]) for line in csv.reader(f)}


def local_blast_query(query, db='gencode_tracks_reversed', strand='both'):
    """
    Blast fasta query via NCBI blastn.
    :param query: Fasta string
    :param strand: Strand to search.  Default is both.  Options: plus, minus, both
    :return: a handle to an xml formatted blast result
    """
    output_handle = NamedTemporaryFile(mode='w+',
                                       suffix='.xml',
                                       encoding='utf-8')
    with NamedTemporaryFile("w") as f:
        f.write(query)
        f.flush()
        res = blastn(query=f.name,
                     db=db,
                     out=output_handle.name,
                     outfmt=5,
                     word_size=11,
                     num_threads=8,
                     strand=strand)
        res()
    return output_handle


def parse_hits(handle, strand=1, match_thresh=13):
    """
    Parses the results from a blast search.  Returns all blast hits on a given strand
    :param handle: results xml handle from blast search
    :param strand: which strand to search for hits.  By default searches complementary strand (-1).
    :return: dictionary of query:[target...] hits for each query
    """
    #Parse blast hits
    gene_hits = defaultdict(list)
    for record in NCBIXML.parse(handle):
        gene_hits[record.query] = [
        ]  # This was added so probes that couldn't match anything weren't dropped
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                #Check that hit is on opposite strand
                if hsp.frame[1] == strand:
                    block_len = [len(block) for block in hsp.match.split(" ")]
                    if any(True for b in block_len if b >= match_thresh):
                        gene_hits[record.query].append(alignment.hit_def)
    handle.close()
    return gene_hits


def blast_probes(probes, timeout=120, debug=False, ):
    """
    Runs entire blast routine for a given gene returning
    only the complementary hits for a probeset.  Stores result in a db.
    :param debug: Whether to print program output
    :param local: Run local or global blast
    :return:
    """
    res = local_blast_query(fasta_query, )
    gene_hits = parse_hits(res)
    res.close()
    return gene_hits


def probe_to_fasta(probe):
    return ">{Name},{Probe #}\n{Probe (5'-> 3')}\n".format(**probe)


def get_copynum(hits):
    off_target = {name: counts[name] for name in hits if name in counts.keys()}
    false_hits = sum(off_target.values())
    return false_hits
