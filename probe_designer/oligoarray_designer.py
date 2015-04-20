from __future__ import print_function, with_statement, division
import dataset
from multiprocessing import cpu_count
from sh import oligoarray_cl as oligoarray
import arrow 
import csv
from itertools import groupby
import pandas as pd
from random import getrandbits

class Oligoarray(object):

    def check_db(self, gene, masking):
        res = self.table.find(Name=gene, Variants=self.variants, Masking=masking)
        passed = []
        for item in res:
            if arrow.get(item['Date']) >= self.date.replace(years = -1):
                item['Probes'] = pd.DataFrame(list(self.p_table.find(ProbeID=item['ProbeID'])))
                passed.append(item)
        return passed


    def write_db(self, probes):
        probes['ProbeID'] = getrandbits(50)
        probes['Date'] = arrow.now().datetime
        probes['Variants'] = self.variants
        probe_df = probes.pop('Probes')
        self.table.insert(probes)
        probe_df['ProbeID'] = probes['ProbeID']
        for k, probe in probe_df.T.to_dict().items():
            self.p_table.insert(probe)


    def check_database(self, gene, mask_range=(5,)):
        probe_set =  [p for p in [self.check_db(gene, mask) for mask in mask_range]
                        if p]
        if probe_set:
            for probe in probe_set: 
                probes.append(probe)
                pl = [pd.concat([p['Probes'] for p in probe], ignore_index=True).shape[0]
                        or probe in probe_set]
            print("Found unblasted probeset for {} with {} probes".format(gene,pl))
        return probe_set


    def design(self, cds, ml=30, Ml=30, max_tm=100, x_hyb=72):
        probes = []
        fasta = FastaFile()
        fake_mask = 11
        for gene, data in cds.items():
            probe_set = self.check_db(gene, fake_mask)
            if probe_set:
                probes.extend(probe_set)
                print("Found {}".format(gene))
                continue
            for n, seq in enumerate(data['CDS List']):
                if len(seq) < 20:
                    continue
                print("Adding {} cds #{} to design list".format(gene, n))
                #TODO: Put probe design here
                temp_name = "{}${}${}".format(gene, n, data['# Isoforms'])
                fasta.add_seq(temp_name, seq)
        od = OligoarrayDesigner()
        results = od.run(fasta.location, min_length=ml, max_length=Ml, min_tm=70, max_tm=100, cross_hyb_temp=x_hyb,secondary_struct_temp=76,
                min_gc=35, max_gc=70)
        for seq_code, table in results:
            seq = fasta.seqs[seq_code]
            codename = fasta.codes[seq_code]
            name, cds_region, num_isoform = codename.split('$')
            print(seq_code, codename, name)
            table['Name'] = table['Name'].map(lambda x: fasta.codes[x])
            probes.append( {'Name': name, 'Masking': fake_mask,
                            'Target Seq': seq, 'Probes': table,
                            'CDS Region #': n, '# isoforms': data['# Isoforms'],
                            })
            self.write_db(probes[-1])
            print("Designed %i probes for %s cds#:%i with masking %i"
                    % (len(table), name, n, fake_mask))
        fasta.close()
        return probes


    def __init__(self, organism="mouse", variants=True):
        self.organism = organism
        self.variants = variants
        self.db = dataset.connect("sqlite:///db/oligoarray.db")
        self.date = arrow.now()
        self.table = self.db[organism]
        self.p_table = self.db['probes']


class OligoarrayDesigner(object):

    def run(self, input_file, output_file="temp/output.txt", rejected_file="temp/rejected.txt",
            log_file="temp/log.txt", max_oligos=10000,
            min_length=30, max_length=30, max_dist=10000, min_tm=85, max_tm=90,
            secondary_struct_temp=65.0, cross_hyb_temp=65, min_gc=40, max_gc=60,
            prohibited_seqs=[6*'C', 6*'G', 6*'T', 6*'A'],
            num_processors=-1, min_dist=2):
        min_spacing = max_length + min_dist
        prohibited_seqs = ";".join(prohibited_seqs)
        if num_processors == -1:
            num_processors = cpu_count() * 3
        oligoarray(i=input_file, d=self.blast_db, o=output_file, r=rejected_file,
                   R=log_file, n=max_oligos, l=min_length, L=max_length, D=max_dist,
                   t=min_tm, T=max_tm, s=secondary_struct_temp, x=cross_hyb_temp,
                   p=min_gc, P=max_gc, m=prohibited_seqs, N=num_processors,
                   g=min_spacing)
        results = OligoArrayResults(output_file)
        p = []
        for name, passed_probe in groupby(results.parse(), lambda x:x['Name']):
            so_probes = sorted(passed_probe, key=lambda x: x['Probe Position*'])
            probes_table = pd.DataFrame(so_probes)
            probes_table['Probe #'] = probes_table.index
            # probes_table['"Probe (5'-> 3')"] =probes_table['"Probe (5'-> 3')"].map(reverse_complement) 
            p.append((name, probes_table))
        results.close()
        return p
            

    def __init__(self, blast_db="/home/eric/blastdb/old_format/mouse_refseq_rnaDB"):
        self.blast_db = blast_db


class OligoArrayResults(object):

    def close(self):
        import os
        os.remove(self.location)

    def parse(self):
        results = []
        try:
            with open(self.location, "r") as f:
                tsvin = csv.reader(f, delimiter='\t')
                for row in tsvin:
                    off_targets = row[7].count(' ;')
                    if off_targets == 0:
                        results.append({"Name":row[0], 'Probe Position*': int(row[1]),
                                        "Probe (5'-> 3')": row[-1],
                                        "Percent GC":"NA"})
        except:
            pass
        return results
    
    def __init__(self, location):
        self.location = location

            
            
class FastaFile(object):
    def to_fasta(self, gene, seq):
        return ">{}\n{}\n".format(gene, seq)

    def add_seq(self, gene, seq):
        seq_code = str(getrandbits(50))
        self.codes[seq_code] = gene
        self.seqs[seq_code] = seq
        fasta_seq = self.to_fasta(seq_code, seq)
        with open(self.location, "a") as f:
            f.write(fasta_seq)

    def close(self):
        import os
        os.remove(self.location)

    def __init__(self, location="temp/temp_fasta.fasta"):
        self.location = location
        self.seqs = {}
        self.codes = {}
        try:
            self.close()
        except:
            pass
        with open(self.location, "w") as f:
            f.write("")
                
        
seq = """
ATGGCCGAGAATGTGGTGGAACCCGGGCCGCCTTCAGCCAAGCGGCCTAAACTCTCATCTCCGGCCCTCT
CGGCGTCCGCCAGCGATGGCACAGATTTTGGTTCACTGTTTGACCTGGAACATGACTTACCAGATGAATT
AATCAACTCTACAGAATTGGGACTAACCAATGGTGGCGATATCAGTCAGCTTCAGACAAGTCTTGGCATA
GTACAAGATGCAGCCTCGAAACATAAACAGCTGTCAGAACTGCTGAGGTCTGGTAGCTCCCCAAACCTCA
ACATGGGAGTCGGTGGCCCAGGCCAAGCGATGGCCAGCCAGGCCCAACAGAACAGCCCTGGATTAAGTTT
GATAAATAGCATGGTCAAAAGCCCAATGGCACAGACAGGCTTGACTTCTCCAAACATGGGGATTGGCAGT
AGTGGACCAAATCAGGGTCCTACTCAGTCCCCAGCAGGTATGATGAACAGTCCAGTGAACCAGCCTGCCA
TGGGAATGAACACAGGGATGAATGCTGGCATGAATCCTGGAATGTTGGCTGCAGGCAATGGACAAGGGAT
AATGCCCAATCAAGTCATGAACGGTTCCATTGGAGCAGGCCGGGGACGGCCAAACATGCAGTACCCAAAT
GCAGGCATGGGCAATGCTGGCAGTTTATTGACTGAGCCACTACAGCAGGGCTCTCCTCAGATGGGAGGAC
AGCCAGGATTGAGAGGCCCCCAACCACTTAAGATGGGAATGATGAACAATCCCAGTCCTTATGGTTCACC
ATACACTCAGAATTCTGGACAGCAGATTGGAGCAAGTGGCCTTGGTCTCCAAATTCAGACAAAGACTGTT
CTACCAAATAACTTATCTCCATTTGCAATGGACAAAAAGGCAGTTCCTGGTGGGGGAATGCCCAGTATGG
GCCAGCAGCCTACCCCATCGGTCCAGCAGCCAGGCCTGGTGACTCCAGTTGCCGCAGGAATGGGTTCTGG
AGCACACACAGCTGATCCAGAGAAGCGCAAGCTCATCCAGCAGCAGCTTGTTCTCCTTTTACATGCTCAC
AAGTGCCAGCGCCGGGAGCAAGCTAATGGGGAAGTGAGGCAGTGCAACCTTCCTCACTGTCGTACCATGA
"""            

def test(seq):
    fasta = FastaFile()
    fasta.add_seq("gi|123173876:416-7654",seq)
    fasta.add_seq("gi|123173876:416-7654",seq)
    fasta.add_seq("gi|123173876:416-7654",seq)
    fasta.add_seq("gi|123173876:416-7654",seq)
    fasta.add_seq("gi|123173876:416-7654",seq)
    od = OligoarrayDesigner()
    # od.run(fasta.location)
    od.run(fasta.location, min_tm=70, max_tm=100, cross_hyb_temp=72,secondary_struct_temp=76,
            min_gc=35, max_gc=70)
    results = OligoArrayResults("temp/output.txt")

# test(seq)
