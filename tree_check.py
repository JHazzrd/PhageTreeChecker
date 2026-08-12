'''
Script to determine if provided label of new_genus for a phage is true or false. Most logic is based on calculated median genome length of approx. 1,500
phage genuses. If genus is unknown, calculated median of neighbouring genomes is used. Generally, the number of neighbours is not significant for calcuation;
only the neighbouring genus has to be calcuateed (for the target value to then be looked up). Neighbouring genus is considered to be the genus of the closest phage,
unless the closest neighbour is > 0.6 away, in which case it is considered that the query is a new genus. New genus clusters are also calculated, determining clusters of
new_genus within a distance of 0.2. Phages originally under the new_genus label are also checked for size; those under 3kb have an additional label noting that they may
be fragments / incomplete genomes.

Includes also a list of representative phages, encompassing ~1500 genera. These are used for when a tree is created with such a list, for more exhaustive new_genera checking.
These are ignored for the final output, and are generally invisible. The script works in such a way that only a taxmyphage output of the query phage(s) needs to be provided, rather
than running a large taxmyphage of both the representative phages + query. If, for whatever reason, a large taxmyphage summary is provided containing both, the script should
still work.

Input: Name of directory containing output from VipTree + TaxMyPhage

Output: Update of TaxMyPhage output, with phages that have been updated in label (new_genus --> Existing) marked as changed, as well as new genus clusters described
(under the temporary name Maybevirus_suffix)

'''
import dendropy # Tree reading and manipulation
import sys
import pandas as pd # TSV reading
import argparse # Testing with command line flag implementation
from icecream import ic # Debugging
from numpy import median, std, mean # Used for median of unrecognised genus (or new_genus), as well as standard deviation for distance counts
from Bio import SeqIO # Used to parse fasta used to make tree to calculate genome lengths.
import os # File finding

def read_file_tree(file_name):
    '''
    Input: Name of newick file
    Output: Tree
    '''
    try:
        tree = dendropy.Tree.get(path = file_name, schema="newick")
    except:
        raise ValueError("Tree file not found!")
    return tree

def read_fasta_length(file_name, input_df, rep = False):
    '''
    Reads in the fasta used during ViPTree to generate sequence lengths for all of the genomes.
    Special formatting for the representative phages due to inconsistency of names between files (if used).
    '''
    merge_df = []
    if not rep:
        with open(file_name) as handle:
            for record in SeqIO.parse(handle,"fasta"):
                d = {"Genome":str(record.id),
                     "Genome size": int(len(record.seq))}
                merge_df.append(d)
        merge_df = pd.DataFrame(merge_df)
    else:
        with open(file_name) as handle:
            for record in SeqIO.parse(handle,"fasta"):
                d = {"Genome":str(record.id.split('.',1)[0]),
                     "Genome size": int(len(record.seq))}
                merge_df.append(d)
        merge_df = pd.DataFrame(merge_df)
    return pd.merge(input_df, merge_df, on = "Genome")
def read_fasta_name(file_name, input_df):
    '''
    Reads in the fasta used for representative phages to get the name used in each fasta file, as the header names are what is used in ViPTree.
    '''
    with open(file_name) as handle:
        for record in SeqIO.parse(handle,"fasta"):
            trim_name = record.id.split('.',1)
            input_df.replace(trim_name,record.id, inplace = True)
    return input_df
def check_node(name, tree):
    '''
    Input: Name of Phage to be found; string
    Output: Checks whether it exists in the tree; returns the found node
    '''
    node = tree.find_node_with_taxon_label(name)
    if node == None:
        raise ValueError(f"Query {name} not found in tree.")
    else:
        return node

def neighbour_genus(query, neighbour, df, rep_df,distance_matrix):
    '''
    Input: Neighbours of query and query
    Output: Genus of neighbours, neighbour it was calculated from
        
    Determines the genus of closest neighbour; used to determine which genus it should belong to. Selected as closest neighbour, unless
    closest neighbour is far away. Iterates through all neighbours until valid genus is found, else accepts new genus.
    '''
    
    out = "New_genus"
    for n in neighbour:
        #if n != "FO818745":
        if df.at[str(n).replace("'","").replace(" ","_"), "Genus"] != "New_genus" and distance_matrix(query.taxon, n) < 0.6:
            ic(df.at[str(n).replace("'","").replace(" ","_"), "Genus"], neighbour[0])
            return df.at[str(n).replace("'","").replace(" ","_"), "Genus"], n
    return out, None


def neighbour_search(query, tree, distance_matrix, n):
    '''
    Input: Phage query, tree, distance matrix and number of neighbours (default = 5)
    Output: Names of nearest neighbours

    Uses distance matrix to determine leaves with smallest distance. These are considered the closest neighbours.
    '''
    neighbours = []
    neighbour_distance = []
    nodes = []

    for node in tree.taxon_namespace: # Get distances of all leaves to query leaf
        neighbour_distance.append(distance_matrix(query.taxon, node))
        nodes.append(node)
    for i in range(0,n+1): # Finds an extra, as usually includes itself in the closest neighbours
        index_min = min(range(len(neighbour_distance)), key=neighbour_distance.__getitem__) # Index of smallest distance
        neighbours.append(nodes[index_min])
        t = neighbour_distance.pop(index_min)
        ic("Distances:",t)
        nodes.pop(index_min) # Remove once used (then find next smallest)
    if query.taxon in neighbours:
        neighbours.remove(query.taxon) # Remove the query from nneighbours
    for name in neighbours:
        ic("Neighbour node:", name)
    return neighbours

def new_genus_cluster(df, distance_matrix, tree,suffix, d = 0.2, threshold = 4):
    '''
    Predicts clusters of "New Genus", assinging them names Maybevirus_suffixX, allowing them to be distinguished when visualiszing the tree. Singletons
    remain under the label New_genus. Clusters are determined by low distance New_genus neighbours. Using distance of 0.2 for now, though
    a source or any better rationale would be useful. Index of suffix starts at 1. Information is best used in conjunction with labelling that does not include
    Maybevirus, as solitary labels amidst clusters of new_genus can be over-written.


    Input: The altered dataframe for the names of phages and their genus for reference, distance matrix and tree for distance.
    Default value for distance (d) is 0.2; can be changed. Suffix default is NGC (New Genus Cluster); can be changed.
    Output: Adjusted output dataframe, with a number of found clusters, and number of new clusters predicted.
    '''
    label = 1 # Used for labelling the clusters
    adjusted = [] # Ensures once a phage has been changed, it isn't further changed
    df_sub = df[df["Proposed Genus"] == "New_genus"] # Grab only new genus, as to not change any labelled phages
    for index, row in df_sub.iterrows():
        counter = 0
        temp_store = []
        if index not in adjusted:
            query = check_node(index.replace("_"," "), tree) # Obtain the leaf data of the current query
            for node in tree.taxon_namespace:
                string_node = str(node).replace("'","").replace(" ","_")
                if distance_matrix(query.taxon,node) <= d and (string_node not in adjusted) and (string_node in df_sub.index):
                    counter += 1
                    temp_store.append(string_node)
            if counter >= threshold: # If a leaf has (threshold=4) or more new_genus neighbours, mark as a cluster
                for phage in temp_store:
                    df.at[phage,"Proposed Genus"] = ("Maybevirus_"+suffix+str(label)) #Update actual df
                    adjusted.append(phage)
                label += 1
    return df, label-1
        
                
def colour_code_output(df, original_df, original = False):
    '''
    Assigns colours to outputted dataframe, merging new results with orginal dataframe, provided in a format that works with iTOL .tsv upload. Ignores the
    representative phages; too many genera to visually represent!

    Input: Dataframe output of script, original dataframe. Variable 'original' for debugging; if true it only outputs from the original df (will ignore new_genus)
    Output: Dataframe in the iTOL .tsv format, for dataset --> color strip.
    '''
    col = ["#E6194B","#3CB44B","#FFE119","#4363D8","#F58231",
           "#911EB4","#46F0F0","#F032E6","#BCF60C","#FABEBE","#E6BEFF",
           "#9A6324","#FFFAC8","#800000","#AAFFC3","#808000","#000075",
           "#808080","#FFD8B1","#CDBEBE","#A9A9A9","#00FF7F","#1E90FF",
           "#FF4500","#2F4F4F","#9400D3","#00CED1","#4B9CD3",
           "#2A52BE", "#5DADEC","#228B22", "#7FFF00", "#2E8B57","#7B68EE",
           "#8A2BE2", "#BA55D3","#CD5C5C", "#FF6EB4", "#DC143C","#D2691E",
           "#FF8C00", "#B5651D","#F0E130", "#FFD700", "#EEDC82","#708090",
           "#C0C0C0", "#2E2E2E","#FF1493", "#20B2AA", "#B22222", "#7FFFD4", "#DAA520",
            "#ADFF2F", "#6A5ACD", "#FF6347", "#4682B4", "#00BFFF",
            "#8FBC8F", "#FFB6C1", "#A52A2A", "#5F9EA0", "#DE3163",
            "#9ACD32", "#556B2F", "#483D8B", "#AFEEEE", "#F4A460",
            "#708238", "#C71585", "#6495ED", "#D8BFD8", "#FA8072",
            "#B0C4DE", "#F0FFF0", "#C19A6B", "#F5DEB3", "#BC8F8F",
            "#36454F", "#708090", "#556677", "#2C3539", "#4A646C",
            "#1B3B6F", "#0E4D92", "#2E5984", "#5DA1C4", "#7FB5B5",
            "#3A5F0B", "#4CBB17", "#2D6A4F", "#74A12E", "#8DA399",
            "#C46210", "#B87333", "#A97142", "#D4AF37", "#E1A95F",
            "#6A0DAD", "#7F00FF", "#9D4EDD", "#A98BC7", "#5D3FD3",
            "#1C1C1C", "#2F2F2F", "#4F4F4F", "#696969", "#9FA6A0",
            "#7E191B", "#A23E48", "#B04E0F", "#CC7722", "#D27D2D",
            "#355E3B", "#567E3A", "#7A9A01", "#9BBF65", "#B7CE63",
            "#0F5257", "#1C6E8C", "#3A8891", "#6CB4B8", "#9AD1D4",
            "#4B3F72", "#6C4F8C", "#8E6BAF", "#B497BD", "#D8C7DF",
            "#5A4E3C", "#7C6A4F", "#A0896F", "#C6B89E", "#E8DCC2",
            "#2B2B2B", "#474747", "#6E6E6E", "#B8B8B8", "#E5E5E5"]
    # An array of colours; should be (mostly) colour-blind friendly.
    col_dict = {}
    c = 0
    found_genus = [] # Comparisson with representative phages can introduce genera not found in the original df
    o_sub_df = original_df[original_df["Representative_phage"] == False]
    for genus in o_sub_df["Genus"].unique():
        col_dict[genus] = col[c]
        found_genus.append(genus)
        c += 1
    if not original: # A small flag to vaugely get a non-changed version of the tree. Debugging only.
        for genus in df["Proposed Genus"].unique():
            if genus.startswith("Maybevirus") or genus not in found_genus:
                col_dict[genus] = col[c]
                c += 1
    col_dict["New_genus"] = "#000000"
    output = []

    # Add back in the untouched phages
    
    for index, row in o_sub_df[o_sub_df["Genus"] != "New_genus"].iterrows():
        d = {"Phage": index,
             "Colour":col_dict[row["Genus"]],
             "Genus":row["Genus"]
            }
        output.append(d)

    # Add in the changed phage
    if not original: 
        for index, row in df.iterrows():
            if not row["Fragment"]: # Allows visual distinguishing of fragments in itol (without having to also look at genome length values)
                colour = col_dict[row["Proposed Genus"]]
                genus = row["Proposed Genus"]
            else:
                colour = "#FFE1FA"
                genus = "Potential Fragment"
            d = {"Phage": index,
                 "Colour":colour,
                 "Genus":genus
                }
            output.append(d)
    return pd.DataFrame(output)

def size_dataset_output(original_df):
    '''
    Outputs phage size information for iTOL; specific output (label) needs checking. Ignores the representative genomes.
    '''
    temp_df = []
    #o_sub_df = original_df[original_df["Representative_phage"] == False]
    o_sub_df = original_df
    for index, row in o_sub_df.iterrows():
        size = row["Genome size"] / 1000
        if size <= 50:
            label = "=<50"
        elif size > 50 and size < 100:
            label = "51-100"
        elif size > 100 and size <= 150:
            label = "101-150"
        elif size > 150 and size <= 200:
            label = "151-200"
        else:
            label = ">200"
        d = {
            "Phage": index,
            "Length (KB)": row["Genome size"]/1000,
            "Label":label
            }
        temp_df.append(d)
    return pd.DataFrame(temp_df)

def tmp_output(df, original_df):
    '''
    Returns a dataframe formatted as TaxMyPhage output. Information for those previously regarded as a new genus is chosen from a random
    phage in the genus. No information is provided about its species; kept as new_species. Changes the message at the end to read:
    "Determined to be an incmplete read within {genus name}". Also describes any new genus clusters found with optional flag -g. Message reads:
    "Phage determined to be in a new genus temporarliy labelled {Maybevirus_Suffix} due to high similiariy with other phages labelled as new genus"
    '''
    for index, row in df.iterrows():
        if row["Proposed Genus"] != "New_genus" and not row["Proposed Genus"].startswith("Maybevirus"):
            example_row = original_df.loc[original_df["Genus"] == row["Proposed Genus"]].head(1)
            original_df.at[index,"Genus"] = row["Proposed Genus"]
            original_df.at[index,"Message"] = str(f"Determined to be an incomplete read within {row['Proposed Genus']}")
            original_df.at[index,"Realm"] = example_row.iloc[0]["Realm"]
            original_df.at[index,"Kingdom"] = example_row.iloc[0]["Kingdom"]
            original_df.at[index,"Phylum"] = example_row.iloc[0]["Phylum"]
            original_df.at[index,"Class"] = example_row.iloc[0]["Class"]
            original_df.at[index,"Order"] = example_row.iloc[0]["Order"]
            original_df.at[index,"Family"] = example_row.iloc[0]["Family"]
            original_df.at[index,"Subfamily"] = example_row.iloc[0]["Subfamily"]
            full_tax = str(f"r__{example_row.iloc[0]['Realm']};k__{example_row.iloc[0]['Kingdom']};p__{example_row.iloc[0]['Phylum']};c__{example_row.iloc[0]['Class']};o__{example_row.iloc[0]['Order']};f__{example_row.iloc[0]['Kingdom']};sf__{example_row.iloc[0]['Subfamily']};g__{example_row.iloc[0]['Genus']};s__Not Defined Yet")
            original_df.at[index,"Full_taxonomy"] = full_tax
        elif row["Proposed Genus"].startswith("Maybevirus"):
            original_df.at[index,"Genus"] = row["Proposed Genus"]
            original_df.at[index,"Message"] = str(f"Phage determined to be in a new genus temporarily labelled {row['Proposed Genus']} due to high similarity with other phages labelled as new genus.")
        if row["Fragment"]:
            original_message = original_df.at[index,"Message"]
            original_df.at[index,"Message"] = (original_message + " Genome size <3kb, this phage is likely not complete and is potentially a fragment.")
    original_df = original_df[original_df.Representative_phage == False]
    return original_df.drop(["Unnamed: 0","Representative_phage","Phage","DNA_Type"], axis = 1) # Remove redundant values from metadata merge

def phylo_distance_genera(df, tree, distance_matrix, stat_df):
    '''
    Aim is to calculate the average phylogenetic distance between phages within a genus. This function iterates over all unique genera within the sample
    and calculates pairwise distances between all phages within each genus. These values are then stored as a mean, median and standard deviations.

    Can be used with Maybevirus classifications, but provides reduced information for such labels.

    Input: Finalised metadata, tree, distance matrix
    Output: .tsv with statistical information regarding the genera
    '''

    genera = df["Genus"].unique().tolist()
    genera.remove("New_genus")
    outfile = open("phage_mppt.tsv","w")
    outfile.write("Genus\tMedian\tMean\tStandardDeviation\tAv_Size\tAv_GC\n")
    
    for genus in genera:
        phages = df.index[df["Genus"] == genus].tolist()
        if len(phages) > 3:
            c = 1 # Counter to prevent checking self
            distances = []
            for n in phages: # Pairwise checking
                query = check_node(n.replace("_"," "), tree)
                for i in range(c,len(phages)):
                    distances.append(distance_matrix(query.taxon, check_node(phages[i].replace("_"," "),tree).taxon))
                    ic(genus, distance_matrix(query.taxon, check_node(phages[i].replace("_"," "),tree).taxon))
                c += 1
            try:
                size = stat_df.at[genus,"Median Genome Length (KB)"]
                gc = stat_df.at[genus,"Median molGC (%)"]
            except:
                size = None
                gc = None
            outfile.write(f"{genus}\t{median(distances)}\t{mean(distances)}\t{std(distances)}\t{size}\t{gc}\n")
    outfile.close()
                
            
            
def main():
    ## Flag Parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="Path to the folder containing newick file and its metadata; folder should contain VipTree output.")
    parser.add_argument("-v", "--verbose", help="Significantly expand outputs for debugging", action="store_true")
    parser.add_argument("-s", "--strict", help="Enable stricter criteria for determining new genus (within 5%% of median)", action="store_true")
    parser.add_argument("-l","--log", help="Outputs log of information to console, such as number of successful reads", action ="store_true")
    parser.add_argument("-i","--itol", help="Changes the output format to be compatable with iTOL dataset, color strip. Currently supports ~125 separate genera", action="store_true")
    parser.add_argument("-g","--genuscluster", help="Calculate clusters of new genus, under the label Maybevirus.", action="store_true")
    parser.add_argument("-n","--name", help="Provides the suffix added to new genera groups (Maybevirus_SUFFIX), allowing for dataset specifc names. Default is NGC")
    parser.add_argument("-t","--threshold", help="Sets the amount of phages required to form an NGC. Default value is 4")
    parser.add_argument("-q","--query", help="Input the name of a node to output a list of neighours. Bypasses usual function of the script.")
    parser.add_argument("-p","--phylodist", help="Calculate Pylogenetic MPD (Pair-wise  Distance) of genera within inputted tree; stores as phage_mppd.tsv", action="store_true")
    args = parser.parse_args()

    PATH = args.path + "/"
    print(f"Reading from {PATH}...")
    
    if not args.verbose:
        ic.disable()
    if args.strict:
        criteria = 0.95
    else:
        criteria = 0.9
    if args.name == None:
        args.name = "NGC"
    if args.threshold == None:
        args.threshold = 4
    ## File Handling
        
    dir_path = os.path.dirname(os.path.realpath(__file__)) + "/" + PATH
    try: # Read Tree
        for root, dirs, files in os.walk(dir_path):
            for file in files:
                if file == "all.bionj.asc.newick":
                    loc = root+"/"+str(file)
                    tree = read_file_tree(loc)
                    print("Found Tree!\t At: "+ PATH + loc.split(PATH)[1])
    except:
        raise ValueError(f"Expected newick file named 'all.bionj.asc.newick' within {PATH}")
    
    try: # Read in Taxonomy Data of Input Phages
        for root, dirs, files in os.walk(dir_path):
            for file in files:
                if file == "Summary_taxonomy.tsv":
                    loc = root+"/"+str(file)
                    df = pd.read_csv(loc, sep="\t")
                    print("Found Metadata!\t At: "+ PATH + loc.split(PATH)[1])
        df["Representative_phage"] = False # If file not found, this causes the error
    except:
        raise ValueError(f"Expected metadata file named Summary_taxonomy.tsv within {PATH}")
    try:
        df["Genome size"] # Check if df already has genome size calculated (field name based on previous data examples)
        print("Using Genome Size from provided metadata.")
    except:
        try:
            for root, dirs, files in os.walk(dir_path):
                for file in files:
                    if str(file) == "all.fasta":
                        loc = root+"/"+str(file)
                        df = read_fasta_length(loc,df)
            
            fasta_length_check = df["Genome size"].head(1) # Fails if file not found
            print("Found FASTA!\t At: "+ PATH + loc.split(PATH)[1])
        except:
            raise ValueError(f"Expected fasta input within ViPTree Output")
    try:
        rep_df = pd.read_csv("Tree_check_metadata/representative_phages_check_ready.tsv", sep="\t")
        rep_df = rep_df.rename(columns = {"Accession":"Genome","Genome":"DNA_Type"})
        rep_df["Representative_phage"] = True
        rep_df.set_index("Genome",inplace=True)
    except:
        raise ValueError("Expected representative phages found in metadata folder! Ensure 'representative_phages_check_ready.tsv' exists in Tree_check_metadata")
    df.set_index("Genome", inplace=True)
    df = pd.concat([df, rep_df], ignore_index = False, sort=False) # Values of representative phages added, but labelled so they can be removed from output
    try:
        stat_df = pd.read_csv("Tree_check_metadata/1612_phage_genome_stats.tsv", sep="\t")
    except:
        raise ValueError("Expected phage statistics in metadata folder! Ensure '1612_phage_genome_stats.tsv' exists in Tree_check_metadata")
    stat_df.set_index("Genus", inplace=True)
    known_genus = stat_df.index.values.tolist()

    ## Variable Setup
    
    n = 5 # Static neighbour searching; 
    fail_count = 0
    complete_count = 0
    update_count = 0
    output_df = []
    
    ## Main Script
    new_genus_subset = df[df["Genus"] == "New_genus"]
    
    if not new_genus_subset.empty and args.query == None:
        print("\nMaking Distance Matrix...")
        distance_matrix = tree.phylogenetic_distance_matrix()
        print("\nDistance Matrix Made\n")
        
        for index, row in new_genus_subset.iterrows(): # Iterate over phages labelled as new_genus
            correct_flag = False # Reset output flags
            error_flag = False
            frag_flag = False

            test_node = check_node(str(index).replace("_"," "), tree) # Query checked for existence and stored

            ic("Query:",test_node.taxon)

            neighbours = neighbour_search(test_node, tree, distance_matrix, n) # Find closest {n = 5} neighbours
            q_size = df.at[index,"Genome size"] / 1000 # df genome size in bp, not kbp
            
            if len(neighbours) > 0: # If a neighbour has been found (should always be the case)
                
                n_genus, neighbour = neighbour_genus(test_node,neighbours, df, rep_df,distance_matrix) # Determine genus of neighbours; temp labels query as that genus
                ic("Neighbour Genus:",n_genus)
                if n_genus in known_genus and n_genus != "Unclassified": # If median already calculated in statistics dataframe, use that.
                    n_median = stat_df.at[n_genus,"Median Genome Length (KB)"]
                    update_count += 1 # Assume update
                else: # Alternative: New_genus or not in database; median calculated from the phylo tree
                    n_size = []
                    for node in neighbours:
                        n_size.append(int(df.at[str(node).replace("'","").replace(" ","_"),"Genome size"]) / 1000)
                    n_median = median(n_size)
                    if n_genus == "New_genus":
                        correct_flag = True
                    else:
                        update_count += 1 # Otherwise not in database, but still needed updating into new_genus
                ic("Target Size", n_median * criteria)
                ic("Query Size:",q_size)
                # if the query is within criteria% of the median [or larger] and is not phlyogenetically close to its neighbour, return label of new_genus
                if (n_median * criteria) <= q_size and n_genus != "New_genus"and distance_matrix(test_node.taxon, neighbour) > 0.1:
                    
                    update_count -= 1 # Not being updated if being kept as New_genus
                    correct_flag = True
                    n_genus = "New_genus" # Used to override neighbouring genus for output
                    ic("Query marked as new genus")
                if q_size < 3: # If size is less than 3kb, there is a chance it is just a fragment / not complete (and thus has a new_genus tag)
                    frag_flag = True
                complete_count += 1
            else: # if no neighbours, something happened
                fail_count += 1
                error_flag = True    
            d = {"Phage":str(test_node.taxon).replace("'","").replace(" ","_"),
                    "Propose New Genus":correct_flag,
                    "Proposed Genus":n_genus,
                    "Error":error_flag,
                    "Fragment":frag_flag
                }
            output_df.append(d)
        output_df = pd.DataFrame(output_df)
        output_df.set_index("Phage", inplace=True) # Current DB is a simple output reporting the new genus phages.
        
        if args.genuscluster:
            output_df,num_new_genus = new_genus_cluster(output_df, distance_matrix,tree, suffix = args.name, threshold = args.threshold)
        else:
            num_new_genus = "NA"
        if args.itol: # ITOL input (formatting)
            output_df_itol = colour_code_output(output_df, df, False) # For use in ITol (copy/paste tsv into dataset)
            output_df_itol.to_csv("tree_checker_output_itol_genus.tsv",sep="\t", index=False)
            size_df = size_dataset_output(df)
            size_df.to_csv("tree_checker_output_itol_size.tsv",sep="\t",index=False)
        try:
            output_df = tmp_output(output_df,df) # If provided with a TMP output, return an output in the same format (prefered!)
        except:
            print("\nNote: Outputting in basic format: No TMP file used.")
        if args.phylodist:
            print("\nCalculating Phylogenetic MPD")
            phylo_distance_genera(output_df, tree, distance_matrix, stat_df)
        output_df.to_csv("tree_checker_output.tsv",sep="\t")
        if args.log:
            print(f"\nSuccessfully parsed {complete_count} labelled New_genus\nFailed to parse {fail_count}\nUpdated {update_count} to neighbouring genus.\nPredict {num_new_genus} new genera clusters.")
        print("\n>>>\tResults exported to file 'tree_checker_output.tsv'")
    else:
        print("\nNo Phages labelled 'New_genus' detected; nothing to compute.")
        if args.query != None:
            distance_matrix = tree.phylogenetic_distance_matrix()
            print("\nDistance Matrix Made\n")
            test = check_node(str(args.query).replace("_"," "),tree)
            neigh = neighbour_search(test, tree, distance_matrix, n=9)
            print(neigh)
        #size_df = size_dataset_output(df)
        #size_df.to_csv("tree_checker_output_itol_size.tsv",sep="\t",index=False)
    sys.exit(0)
        

if __name__ == "__main__":
    main()
