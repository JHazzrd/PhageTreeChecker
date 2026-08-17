import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
phage_meta = pd.read_csv("representative_phages.tsv",sep="\t")
accessions = list(phage_meta["Accession"].unique())
for i in tqdm(range(0,len(accessions))): #Clean up accessions for finding.
    accessions[i] = str(accessions[i])
    accessions[i] = accessions[i].split(' ',1)[0] # Remove spaces after the accesion
    accessions[i] = accessions[i].split('.',1)[0] # Remove any .1s; not in the fasta file
outfile = open("representative_phage_genomes.fasta","w")
out_text = ""
c = 0

for record in SeqIO.parse("Phage_data/7Apr2026_genomes_excluding_refseq.fa","fasta"):
    if record.id in accessions:
            out_text = out_text + ">"+record.id+"\n"
            out_text = out_text + str(record.seq)+"\n"
            if c % 100 == 0:
                print(c) # Counter for sanity.
            c+=1
            accessions.remove(record.id)
    else:
        pass
    
outfile.write(out_text)
outfile.close()
print(accessions)
print(f"count:{c}")
print(len(accessions))
            
    
    


