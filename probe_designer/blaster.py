import collections
import signal
import sys
from collections import Counter, defaultdict
from time import sleep

import dataset
import pandas as pd
from Bio.Blast import NCBIWWW, NCBIXML
from future.builtins import object

blast_db = dataset.connect("sqlite:///db/blast.db")

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


def local_blast_query(query, max_time=120, max_iterations=10, organism='"Mus musculus"[porgn:__txid10090]'):
    """
    Blast fasta query via NCBI qblast.
    Requires web access.  Qblast often times out so request is killed after max_time
    and repeated for max_iterations.
    :param query: Fasta string
    :param max_time: Maximum amount of time in seconds to wait for a response from NCBI
    :param max_iterations: Maximum times to request from NCBI
    :return: a handle to an xml formatted blast result
    """
    from Bio.Blast.Applications import NcbiblastnCommandline as blastn
    query_file = "temp/blast_query.fasta"
    output_file = "temp/blast_res.xml"
    with open(query_file, "w") as f:
        f.write(query)
    res = blastn(query=query_file, db="refseq_rna_mouse", out=output_file, outfmt=5,
                 word_size=11, num_threads=8)
    res()
    print("DONE BLAST INNER LOOP")
    output_handle = open(output_file, 'r')
    return output_handle

def test():
    fasta_query = ">BOB1\naaatccgtgtttatccgatatgttgttggtgagtttc"
    organism='"Mus musculus"[porgn:__txid10090]'
    res = local_blast_query(fasta_query, max_time=timeout, organism=organism)

def blast_query(query, max_time=120, max_iterations=10, organism='"Mus musculus"[porgn:__txid10090]'):
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
                                   entrez_query=organism,
                                   word_size=11)
            except Exception as inst:
                print(type(inst))     # the exception instance
                print(inst.args)      # arguments stored in .args
                print(inst)          # __str__ allows args to printed directly
                print("Failed on iteration %i" % iterations)
                if iterations >= max_iterations:
                    raise BlastError("Timeout while blasting query")
                sys.stdout.flush()
                sleep(5)
        if res:
            return res
        iterations += 1


def parse_hits(handle, strand=-1, match_thresh=13):
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
                    block_len = [len([1 for digit in block if digit == "|"])
                                    for block in hsp.match.split(" ")]
                    if any([True for b in block_len if b >= match_thresh]):
                        gene_hits[record.query].append(alignment.hit_def)
    return gene_hits


def check_db(probe_df, organism):
    blast_table = blast_db[organism]
    res = {}
    for k, v in probe_df.iterrows():
        query = v['Probe Name']
        hits = list(blast_table.find(probe_name=query))
        if any(hits):
            res[query] = [match['hit'] for match in hits]
    return res

def write_db(hits, organism):
    blast_table = blast_db[organism]
    insert = [{'probe_name': query, 'hit': hit}
              for query, v in hits.iteritems()
              for hit in v]
    blast_table.insert_many(insert)

def blast_probes(gene, probe_df, timeout=120, debug=False, organism='"Mus musculus"[porgn:__txid10090]',
                 local=False):
    """
    Runs entire blast routine for a given gene returning
    only the complementary hits for a probeset.  Stores result in a db.
    :param gene: Name of gene
    :param probe_df: Pandas DataFrame of probes
    :param debug: Whether to print program output
    :param local: Run local or global blast
    :return:
    """
    # db_gene_hits = check_db(probe_df, organism)
    db_gene_hits = {}
    exists = probe_df['Probe Name'].isin(db_gene_hits.keys())
    blast_these = exists != True
    blast_df = probe_df.ix[blast_these]
    if not blast_df.empty:
        print("Blasting {} of {} oligos".format(len(blast_df), len(probe_df)))
        fasta_query = probe_df_to_fasta(blast_df)
        if local == True:
            if '10090' not in organism:
                raise Exception("Organism not supported for local blast")
            print("DOING LOCAL BLAST")
            res = local_blast_query(fasta_query, max_time=timeout, organism=organism)
            print("DONE LOCAL BLAST, PARSE AWAY")
            gene_hits = parse_hits(res)
            print("DONE PARSING")
            res.close()
        else:
            res = blast_query(fasta_query, max_time=timeout, organism=organism)
            gene_hits = parse_hits(res)
        write_db(gene_hits, organism)


def flatten(l):
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el

class BlastError(Exception):
    def __init__(self, exception):
        Exception.__init__(self, exception)

def filter_probes_2(gene, blast_hits, probe_df, max_probes=24, min_probes=16, n_probes=24, max_false_hits=7, debug=False):
    gene = gene.strip()
    #Count all genes that have off target hits
    d = pd.DataFrame({k:Counter(v) for k, v in blast_res.iteritems()})
    d = d.fillna(0)
    tot_hits = d.sum(1)
    tot_hits.sort(ascending=False)
    tot_hits


def filter_probes_based_on_blast(gene, blast_hits, probe_df, max_probes=24, min_probes=16, n_probes=24, max_false_hits=7, debug=False):
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
    i = Counter([hit for k, v in blast_hits.items() for hit in v])
    s = pd.Series({k:v for k,v in i.items()})
    s.sort(ascending=False)
    p = []
    #Determines a gene is a match if its abbreviation is found in description and is the same as the given name
    for desc in s.index:
        for name in flatten(map(lambda x:x.split(")"), desc.split("("))):
            if name.lower() == gene.lower():
                p.append(False)
                break
        else:
            p.append(True)
    bad_genes = s.index[p]
    bad = s[bad_genes]
    #Tallies up in order of how many incorrect genes are hit all the bad hits
    eb = [bad[bad > i].index for i in range(1, 10)]
    bad_hits = [ [p] + [sum(badi.isin(hits)) for badi in eb]
           for p, hits in blast_hits.items()]
    #Creates tuple=(probe_name, # of off target hits, # of 2 off target hits )
    bad_hits = sorted(bad_hits, key=lambda x: [x[i] for i in range(1, 10)])
    # number of probes to design
    n_probes = min((len(blast_hits), max_probes))
    while n_probes >= min_probes:
        select_probes = [x[0] for x in bad_hits[:n_probes]]
        i2 = Counter([hit for probe in select_probes for hit in blast_hits[probe]])
        s2 = pd.Series({k:v for k,v in i2.items()})
        s2.sort(ascending=False)
        if debug:
            print(s2[s2>1], "\n")
        if s2[bad_genes].max() > max_false_hits:
            print("Bad probeset for %s at %i probes\n" % (gene, n_probes))
        else:
            print("Good probeset designed for %s at %i probes\n" % (gene, n_probes))
            break
        n_probes -= 1
    else:
        raise BlastError("Couldn't Design at least {} probes for {}".format(min_probes, gene))
    passed_df = pd.concat([probe_df.ix[probe_df['Probe Name']==probe]
                           for probe in select_probes])
    return passed_df


def filter_probes2(gene, blast_hits, probe_df, max_probes=24, min_probes=16, n_probes=24, max_false_hits=7, debug=False):
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
    i = Counter([hit for k, v in blast_hits.items() for hit in v])
    s = pd.Series({k:v for k,v in i.items()})
    s.sort(ascending=False)
    p = []
    #Determines a gene is a match if its abbreviation is found in description and is the same as the given name
    for desc in s.index:
        for name in flatten(map(lambda x:x.split(")"), desc.split("("))):
            if name.lower() == gene.lower():
                p.append(False)
                break
        else:
            p.append(True)
