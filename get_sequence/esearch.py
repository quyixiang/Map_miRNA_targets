import sys
target=sys.argv[1:]
target=" ".join(target)
from Bio import Entrez
Entrez.email = "quyixiang666@gmail.com"
result=Entrez.esearch(db="Nucleotide", term=target)
record = Entrez.read(result)
final=record["IdList"][0]
print(final) 
