import sys
from screed import fasta

in_path1 = sys.argv[1]
in_path2 = sys.argv[2]

def count_snps(string1, string2):
    count=0
    length=0
    for letter1, letter2 in zip(string1, string2):
        #print(letter1, letter2)
        length += 1
        if letter1 != letter2:
            count+=1
    #print(count/length)i
    return(count/length)

with open(in_path1, "r") as fin1:
    record_count=0
    total_distance=0
    record_1_1=""
    record_1_2=""
    with open(in_path2, "r") as fin2:
        record_2_1=""
        record_2_2=""
        for record1, record2 in zip(fasta.fasta_iter(fin1), fasta.fasta_iter(fin2)):
            if record_count % 2 == 0:
                record_1_1=record1
                record_2_1=record2
                record_count+=1
            else:
                distance=0
                record_1_2=record1
                record_2_2=record2
                distance += count_snps(record_1_1.sequence, record_2_1.sequence)
                distance += count_snps(record_1_2.sequence, record_2_1.sequence)
                distance += count_snps(record_1_1.sequence, record_2_2.sequence)
                distance += count_snps(record_1_2.sequence, record_2_2.sequence)
                distance = distance / 4
                #print(distance)
                record_count+=1
                total_distance+=distance
    print(total_distance/(record_count/2))

