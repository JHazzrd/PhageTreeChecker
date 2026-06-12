import os, sys
#from tqdm import tqdm
from Bio import SeqIO
if False:
    with open ("rep_phage_outliers.tsv","r") as file:
        lines= file.readlines()
    for i in tqdm(range(0,len(lines))):
        os.system(f"esearch -db nuccore -query {lines[i].strip()} | efetch -format fasta | tee -a output.txt")
else:
    c = 0
    with open("Tree_check_metadata/rep_phage_gen_complete.fasta") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id == "AF234172":
                print(record.id)
                print(record.seq)
                print(record)
    print(c)
