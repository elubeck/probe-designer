from functools import wraps
import errno
import os
import signal
from collections import Counter, defaultdict
from Bio.Blast import NCBIWWW, NCBIXML
import sys
from time import sleep
import pandas as pd

class TimeoutError(Exception):
    pass

class timeout:
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
    p = []
    for k, v in probe_df.iterrows(): 
        p.append('>%s\n%s' 
             %(v['Probe Name'], v["Probe (5'-> 3')"]))
    fasta = '\n'.join(p)
    return fasta

def blast_query(query, max_time=120, max_iterations=10):
    # Loop with timeout while blasting incase server doesn't respond
    q = None
    iterations = 0
    while True:
        with timeout(seconds=max_time):
            try:
                q = NCBIWWW.qblast('blastn', 'refseq_rna', query,
                               entrez_query='"Mus musculus"[porgn:__txid10090]',
                               word_size=7)
            except:
                print("Failed on iteration %i" % iterations)
                if iterations >= max_iterations:
                    raise Exception("Timeout while blasting query")
                sys.stdout.flush()
                sleep(5)
        if q:
            return q
        iterations+=1


def blast_probes(gene, probe_df, debug=False):
    fasta_query = probe_df_to_fasta(probe_df)
    q = blast_query(fasta_query)
    #Parse blast hits
    gene_hits = defaultdict(list)
    for record in NCBIXML.parse(q):
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                #Check that hit is on opposite strand
                if hsp.frame[1] == -1:
                    gene_hits[record.query].append(alignment.hit_def)
    if debug:
        #Print the false hits for an entire gene
        i = Counter([hit for k, v in gene_hits.iteritems() for hit in v])
        s = pd.Series({k:v for k,v in i.iteritems() if v>1})
        s.sort(ascending=False)
        print(gene, s)
        print('...\n'*3)
    return gene_hits
        ### Naievely reduces probes to minimize global false hits for a given probeset

import collections
def flatten(l):
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el


def filter_probes_based_on_blast(gene, blast_hits, probe_df, max_false_hits=8):
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
    select_probes = map(lambda x:x[0], bad_hits[:24])
    i2 = Counter([hit for probe in select_probes for hit in blast_hits[probe]])
    s2 = pd.Series({k:v for k,v in i2.iteritems()})
    s2.sort(ascending=False)
    if s2[bad_genes].max()>8:
        raise Exception("Bad probeset for %s" % gene)
    passed_df = pd.concat([probe_df.ix[probe_df['Probe Name']==probe]
                           for probe in select_probes])
    return passed_df
    
