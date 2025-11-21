import sys
from screed import fasta

in_path1 = sys.argv[1]

with open(in_path1, "r") as fin1:
    record_count=0
    for record1 in fasta.fasta_iter(fin1):
        if record_count % 2 == 0:
            record_count+=1
        else:
            record_count+=1
            print(f">{record1.name}")
            print(f"{record1.sequence}")

