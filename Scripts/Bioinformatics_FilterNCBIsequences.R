# Date: March 26th 2026
# Purpose: Clean up NCBI labyrinthula fasta sequences
# Seagrass Disease in Palau

##### Load Packages ############################################################
library(tidyverse)
library(seqinr)

##### Load data ################################################################
fasta_file <- read.fasta(file = "G:/My Drive/_PhD/Research/Seagrass in Palau/Genetics/Extracting Labyrinthula 18S from NCBI/Search_1_combo.fasta.txt")

##### Filter out duplicates ####################################################
# check for duplicated names
names <- names(fasta_file)
which(duplicated(names)) # no duplicated names

# check for duplicated sequences
sequences <- getSequence(fasta_file)
index_for_duplications <- which(duplicated(sequences)) # many! 

# Remove duplications 
names_dup <- names[index_for_duplications] # find the names of those duplicated 
unique_fasta <- fasta_file[-index_for_duplications]
length(unique_fasta)
length(index_for_duplications)

##### Filter length of sequences ###############################################
unique_fasta_long <- unique_fasta[getLength(unique_fasta) > 800]
length(unique_fasta_long) # 171 unique sequences greater than 800 base pairs 

##### Remove indivudal accession numbers #######################################
unique_names <- names(unique_fasta_long)

# Find and remove the accession numbers from Table S1 of Agnew-Camiener et al. 2025 
# https://doi.org/10.1111/jeu.13073

accession_names <- substr(unique_names, 1, nchar(unique_names) - 2)
n_to_remove <- which(accession_names == "AF348522" | # "rRNA-like"
        accession_names == "EU431330" | #"Non-axenic culture"
        accession_names == "EU431329") #"Non-axenic culture"

fasta_final <- unique_fasta_long[-n_to_remove]
length(fasta_final)

##### Save the filtered sequences ##############################################
write.fasta(sequences = fasta_final, 
            names = names(fasta_final),
            file.out = "G:/My Drive/_PhD/Research/Seagrass in Palau/Genetics/Extracting Labyrinthula 18S from NCBI/Search_1_filtered.fasta")

# Next step is to align and then calculate how similar they are to remove those with >= 99% similarity 
