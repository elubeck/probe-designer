import os
import sys

from probe_designer.mRNA_designer import RNARetriever2

retriever = RNARetriever2()

pgk1 = retriever.get_single_rna("Pgk1", cds_only=True)
print("BOB")
print(len(pgk1),)


