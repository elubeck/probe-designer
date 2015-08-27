from __future__ import print_function, division
import dataset
from subprocess import call
from tempfile import NamedTemporaryFile
from multiprocessing import Pool, cpu_count
from progressbar import ProgressBar
from itertools import imap
import csv


def find_chromosome(gene):
    db = dataset.connect("sqlite:///db/refGene.db")
    table = db['mouse']
    # Some genes may hit annotated chromosome fragments.  Ignore these fragments to prevent errors
    hits = {r['chrom'].split('_')[0] for r in table.find(name2=gene)}
    if len(hits) != 1:
        raise Exception()
    return list(hits)[0]


def find_txStart(gene):
    db = dataset.connect("sqlite:///db/refGene.db")
    table = db['mouse']
    hits = {r['txStart'] for r in table.find(name2=gene)}
    if len(hits) != 1:
        raise Exception()
    return list(hits)[0]


def blat_probes(probe_list, chromosome):
    psl_header = ['matches', 'misMatches', 'repMatches', 'nCount',
                  'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert',
                  'strand', 'qName', 'qSize', 'qStart', 'qEnd', 'tName',
                  'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes',
                  'qStarts', 'tStarts', 'queryAligns', 'subjectAligns']
    db_path = "db/mouse/chromFaMasked/{}.fa.masked".format(chromosome)
    res = []
    with NamedTemporaryFile() as fasta_file:
        fasta_str = "\n".join(">{}\n{}".format(n, items)
                              for n, items in enumerate(probe_list))
        fasta_file.write(fasta_str)
        fasta_file.flush()
        with NamedTemporaryFile() as f:
            p_c = call(['./blat', db_path, fasta_file.name, f.name, '-q=dna',
                        '-out=pslx'])
            for n, line in enumerate(csv.reader(f, delimiter='\t')):
                if n < 5:
                    continue
                vals = dict(zip(psl_header, line))
                query = {'query': probe_list[int(vals['qName'])]}
                vals.update(query)
                res.append(vals)
    return res


def blat_wrapper(d):
    try:
        gene, passed = d
        gene = gene[0].upper() + gene[1:]
        if not passed:
            return {}
        chr = find_chromosome(gene)
        chr_positions = blat_probes(passed, chr)
        for pos in chr_positions:
            pos.update({'target': gene})
        return chr_positions
    except:
        return {}


def blat_pset(p_set, f_name):
    # TODO: ORGANIZE BY CHROMOSOME
    n_genes = len(p_set)
    p_bar = ProgressBar(maxval=n_genes)
    p = Pool(processes=cpu_count())
    with open('temp/{}'.format(f_name), 'wb') as f_out:
        bed = csv.writer(f_out, delimiter='\t')
        for n, chr_positions in enumerate(p.imap(blat_wrapper, p_set.items())):
            for align in chr_positions:
                bed.writerow([align['tName'], align['tStart'], align['tEnd'],
                              align['strand'],
                              align['target'] + '-' + align['qName']])
            p_bar.update(n)
