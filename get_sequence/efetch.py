import sys
target=sys.argv[1]

from Bio import Entrez
Entrez.email = "quyixiang666@gmail.com"
handle = Entrez.efetch(db="nucleotide", id=target, rettype="fasta", retmode="text")
print (handle.read())