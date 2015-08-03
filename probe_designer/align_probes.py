from __future__ import print_function, division
import dataset
import intron_designer2
from subprocess import call
from tempfile import NamedTemporaryFile
import csv


def find_chromosome(gene):
    db = dataset.connect("sqlite:///db/refGene.db")
    table = db['mouse']
    hits = {r['chrom'] for r in table.find(name2=gene)}
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
    db_path = "db/chromFaMasked/{}.fa.masked".format(chromosome)
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

# #TODO:  Blast all probes and search for off target hits
# gene = "Pgk1"
# import dataset
# intron_db_filtered = dataset.connect(
#     "sqlite:///db/intron_probes_filtered_10k.db")
# filtered_probe_table = intron_db_filtered['mouse']
# passed = [probe['seq'] for probe in filtered_probe_table.find(target=gene)]

# chr = find_chromosome(gene)
# txStart = find_txStart(gene)
# chr_positions = blat_probes(passed, chr)
# iStart = int(txStart)
# for probe in chr_positions:
#     print(iStart - int(probe['tStart']))