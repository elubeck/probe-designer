import signal
from collections import Counter, defaultdict
import sys
from time import sleep
import collections

from Bio.Blast import NCBIWWW, NCBIXML
import pandas as pd


class TimeoutError(Exception):
    pass


class timeout(object):
    def __init__(self, seconds=1, error_message='Timeout'):
        self.seconds = seconds
        self.error_message = error_message
    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)
    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)
    def __exit__(self, type, value, traceback):
        signal.alarm(0)


def probe_df_to_fasta(probe_df):
    """
    Converts a Dataframe of probes to a fasta string.
    Dataframe must have two collumns, 'Probe Name' and 'Probe (5'-> 3')'.
    :param probe_df: Pandas DataFrame of Probes.
    :return: a fasta string
    """
    p = []
    for k, v in probe_df.iterrows():
        p.append('>%s\n%s'
                 % (v['Probe Name'], v["Probe (5'-> 3')"]))
    fasta = '\n'.join(p)
    return fasta


def blast_query(query, max_time=120, max_iterations=10):
    """
    Blast fasta query via NCBI qblast.
    Requires web access.  Qblast often times out so request is killed after max_time
    and repeated for max_iterations.
    :param query: Fasta string
    :param max_time: Maximum amount of time in seconds to wait for a response from NCBI
    :param max_iterations: Maximum times to request from NCBI
    :return: a handle to an xml formatted blast result
    """
    res = None
    iterations = 0
    while True:
        with timeout(seconds=max_time):
            try:
                res = NCBIWWW.qblast('blastn', 'refseq_rna', query,
                                   entrez_query='"Mus musculus"[porgn:__txid10090]',
                                   word_size=7)
            except Exception as inst:
                print(type(inst))     # the exception instance
                print(inst.args)      # arguments stored in .args
                print(inst)          # __str__ allows args to printed directly
                print("Failed on iteration %i" % iterations)
                if iterations >= max_iterations:
                    raise Exception("Timeout while blasting query")
                sys.stdout.flush()
                sleep(5)
        if res:
            return res
        iterations += 1


def parse_hits(handle, strand=-1):
    """
    Parses the results from a blast search.  Returns all blast hits on a given strand
    :param handle: results xml handle from blast search
    :param strand: which strand to search for hits.  By default searches complementary strand (-1).
    :return: dictionary of query:[target...] hits for each query
    """
    #Parse blast hits
    gene_hits = defaultdict(list)
    for record in NCBIXML.parse(handle):
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                #Check that hit is on opposite strand
                if hsp.frame[1] == strand:
                    gene_hits[record.query].append(alignment.hit_def)
    return gene_hits


def blast_probes(gene, probe_df, timeout=120, debug=False):
    """
    Runs entire blast routine for a given gene returning
    only the complementary hits for a probeset.
    :param gene: Name of gene
    :param probe_df: Pandas DataFrame of probes
    :param debug: Whether to print program output
    :return: dictionary of probe:[target...] hits for each probe
    """
    fasta_query = probe_df_to_fasta(probe_df)
    res = blast_query(fasta_query, max_time=timeout, )
    gene_hits = parse_hits(res)
    if debug:
        #Print the false hits for an entire gene
        i = Counter([hit for k, v in gene_hits.iteritems() for hit in v])
        s = pd.Series({k:v for k,v in i.iteritems() if v>1})
        s.sort(ascending=False)
        print(gene, s)
        print('...\n'*3)
    return gene_hits


def flatten(l):
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el


def filter_probes_based_on_blast(gene, blast_hits, probe_df, n_probes=24, max_false_hits=7, debug=False):
    """
    Attempts to pick the n_probes best probes for a given set of probes.
    :param gene: String - name of gene
    :param blast_hits: dict probe:[target...] for each probe in probe_df
    :param probe_df: Pandas DataFrame of Probes
    :param n_probes: # of probes to design for gene
    :param max_false_hits: Maximum # of off target hits for a given probeset
    :param debug: Whether to print filtering
    :return:
    """
    gene = gene.strip()
    #Make an index of bad blast hits
    i = Counter([hit for k, v in blast_hits.iteritems() for hit in v])
    s = pd.Series({k:v for k,v in i.iteritems()})
    s.sort(ascending=False)
    p = []
    for desc in s.index:
        for name in flatten(map(lambda x:x.split(")"), desc.split("("))):
            if name.lower() == gene.lower():
                p.append(False)
                break
        else:
            p.append(True)
    bad_genes = s.index[p]
    bad = s[s.index[p]]
    extra_bad = bad[bad > 2].index
    bad_hits = [(probe, sum(bad_genes.isin(probe_hits)), sum(extra_bad.isin(probe_hits)))
                 for probe, probe_hits in blast_hits.iteritems()]
    bad_hits = sorted(bad_hits, key=lambda x: (x[1], x[2]))
    select_probes = map(lambda x:x[0], bad_hits[:n_probes])
    i2 = Counter([hit for probe in select_probes for hit in blast_hits[probe]])
    s2 = pd.Series({k:v for k,v in i2.iteritems()})
    s2.sort(ascending=False)
    if debug:
        print(s2[s2>1])
    if s2[bad_genes].max() > max_false_hits:
        print("Bad probeset for %s" % gene)
    passed_df = pd.concat([probe_df.ix[probe_df['Probe Name']==probe]
                           for probe in select_probes])
    return passed_df
