from __future__ import print_function
from Bio import SeqIO
from Bio.Seq import Seq
import gzip
import csv
from oligoarray_designer import Oligoarray

class Introns(object):
    """
    Iterator that iterates over introns retrieved from UCSC.
    """

    def get_ids(self):
        # Table of every RefSeqID -> GeneName
        with open("db/Refseq2Gene.tsv", "r") as tsvin:
            return dict(csv.reader(tsvin, delimiter='\t'))

    def get_introns(self, repeat_masking=True):
        # Refseq intron fasta indexed by ID downloaded from UCSC genome browser
        # Repetive elements are masked in lowercase
        # By default repeats are masked
        used_names = []
        from collections import Counter
        with gzip.open("db/refseq_intron_maskedLower.gzip.gz", 'rb') as f:
            for fasta in SeqIO.parse(f, 'fasta'):
                id = fasta.id.split("refGene_")[1]
                ref_name = self.refseq_id[id]
                used_names.append(ref_name)
                intron_num = used_names.count(ref_name)
                fasta.name = "{}_Intron#{}".format(ref_name, intron_num)
                if repeat_masking:
                    fasta.seq = Seq("".join([c for c in fasta.seq if c.isupper()]))
                yield fasta
                
    def __iter__(self):
        for intron in self.get_introns(self.repeat_masking):
            yield intron
    
    def __init__(self):
        self.refseq_id = self.get_ids()
        self.repeat_masking = True
        self.organism = "mouse"

            
class Intron(object):
    def to_fasta(self, chunk_size=1000):
        rev_comp = self.record.seq.reverse_complement()
        for chunk_n, chunk_start in enumerate(range(0, len(rev_comp), 1000)):
            seq_chunk = rev_comp[chunk_start: chunk_start+1000]
            fasta_str = ">{}=Chunk:{}\n{}\n".format(self.name, chunk_n, seq_chunk)
            yield fasta_str

    def __init__(self, record):
        self.record = record
        # TODO: Embed Intron number in name
        self.name = "{}_{}".format(self.record.name, "intron")
    

if __name__ == "__main__":
    introns = Introns()
    from oligoarray_designer import OligoarrayDesigner
    queue = []
    for n, i in enumerate(introns):
        print(i.name)
        i.name = "{}_{}_{}".format(i.name, n, "intron")
        print(len(i.format("fasta")))
        queue.append(i)
        if (n+1)%10 == 0:
            f_loc = "temp/temp.fasta"
            with open(f_loc, "w") as f:
                for intron in queue:
                    rev_comp = intron.seq.reverse_complement()
                    for chunk_n, chunk_start in enumerate(range(0, len(rev_comp), 1000)):
                        seq_chunk = rev_comp[chunk_start: chunk_start+1000]
                        fasta_str = ">{}=Chunk:{}\n{}\n".format(intron.name, chunk_n, seq_chunk)
                        f.write(fasta_str)
            o = OligoarrayDesigner()
            res = o.run(f_loc, max_dist=1000, min_length=35, max_length=35, max_tm=100,
                        cross_hyb_temp=72, min_tm=74, secondary_struct_temp=76)
            for name, probes in res:
                print(name, len(probes))
                probes.to_csv("oligoarray_probes/{}.csv".format(name))
            queue = []
    print(n)
            