import pandas as pd

PATH = "Phage_data/"

df_meta = pd.read_excel(PATH+"VMR_MSL41.v1.20260320.xlsx", index_col = 23)
output_df = []

subset = df_meta[df_meta["Host source"] == "bacteria"]

for genus in subset["Genus"].unique():
    example_phage = subset.loc[subset["Genus"] == genus].head(2) # Find two representive phages in genus
    for index, row in example_phage.iterrows():
        if "RNA" not in row["Genome"]: #quick fix
            d = {"Genus":genus,
                 "Subfamily":row["Subfamily"],
                 "Family":row["Family"],
                 "Accession":str(index),
                 "Phage":row["Virus name(s)"],
                 "Genome":row["Genome"]}
            output_df.append(d)
output_df = pd.DataFrame(output_df)
output_df.to_csv("representative_phages.tsv",sep="\t",index=False)

# Does lead to a single error result: A Punavirus with Accession Number NaN. Made a note to
# manually replace this with the next Punavirus: FO818745 ; Escherichia phage RCS47




