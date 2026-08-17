import os, sys
#from tqdm import tqdm
from Bio import SeqIO
with open ("rep_phage_outliers.tsv","r") as file:
    ines= file.readlines()
for i in tqdm(range(0,len(lines))):
    os.system(f"esearch -db nuccore -query {lines[i].strip()} | efetch -format fasta | tee -a output.txt")
