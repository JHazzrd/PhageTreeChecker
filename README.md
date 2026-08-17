A repo containing my Bioinformatics Masters project: A script for checking phage phylogenetic trees. 

_TreePhix_ is a script that is designed to read in results from taxMyPhage (Millard et al., 2025) and VipTreeGen (Nishimura et al., 2017) in order to update incorrect labels of new genera caused by the mislabelling of incomplete bacteriophage genomes. _TreePhix_ relies on phylogenetic distance matricies to determine leaf pair-wise distance, as well as bacteriophage genome size, in order to correct mislabelled 'new genera' to the closest neighbouring genus.

Use:

In terminal: python treephix.py _directory_containing_input_files_

This directory must contain three files:
* Summary_taxonomy.tsv - From taxMyPhage output
* all.bionj.asc.newick - From VipTreeGen output
* all.fasta            - From VipTreeGen output

The names of these files must be left unchanged, though _TreePhix_ will search for these files within the directory provided.
For ease of use, it is recommended to run _TreePhix_ in the same directory as where the input directory is stored.

_TreePhix_ also has various input flags to adjust outputs which can be utilised. This should be added after the directory is provided:
* -h, --help            show this help message and exit
*  -v, --verbose         Significantly expand outputs for debugging
*  -s, --strict          Enable stricter criteria for determining new genus (within 5% of median)
*  -l, --log             Outputs log of information to console, such as number of successful reads
*  -i, --itol            Changes the output format to be compatible with iTOL dataset, colour
                        strip. Currently supports ~125 separate genera
*  -g, --genuscluster    Calculate clusters of new genus, under the label Maybevirus.
*  -n NAME, --name NAME  Provides the suffix added to new genera groups (Maybevirus_SUFFIX),
                        allowing for dataset specifc names. Default is NGC
*  -t THRESHOLD, --threshold THRESHOLD
                        Sets the amount of phages required to form an NGC. Default value is 4
*  -q QUERY, --query QUERY
                        Input the name of a node to output a list of neighours. Bypasses usual
                        function of the script.
*  -p, --phylodist       Calculate Pylogenetic MPD (Pair-wise Distance) of genera within inputted
                        tree; stores as phage_mppd.tsv

Though _TreePhix_ does not require the list of representative genomes to work, the process of how this was generated is shown via the three extra scripts:
* representative_phage_get.py generates a list of representative phage genomes from the VMR dataset
* representative_phage_genome_get retrieves the genomes associated
* rep_phage_outliers_get.py retrieves the genomes that could not be found in the above process

Please note that each script has numerous comments to state how errors should be dealt with, as some bacteriophages could not be found. These should be followed
as to produce the same set of genomes used. 

A full list of phage genomes used to create the representative phage list can be found in Tree_check_metadata/representative_phages_check_ready.tsv. A full list of all accession numbers can be found,
as well as basic taxonomic information.

treephix_testdata includes files for testing _TreePhix_; this can be done by parsing this directory as input for _TreePhix_. 

References:
Millard, A. et al. (2025) 'taxMyPhage: Automated taxonomy of dsDNA phage genomes at the genus and species level', PHAGE (New Rochelle, N.Y.), 6(1), pp. 5–11. Available at: https://doi.org/10.1089/phage.2024.0050 .

Nishimura, Y. et al. (2017) 'ViPTree: The viral proteomic tree server', Bioinformatics (Oxford, England), 33(15), pp. 2379–2380. Available at: https://doi.org/10.1093/bioinformatics/btx157 . 


