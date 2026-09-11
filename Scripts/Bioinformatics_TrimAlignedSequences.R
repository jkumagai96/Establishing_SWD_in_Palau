# Date: March 18th 2026
# Purpose: Open aligned fasta and trim to 100 coverage
# Seagrass Disease in Palau

##### Load Packages ############################################################
library(tidyverse)
library(seqinr)

##### Run Palau only aligned trim ###
fasta_file_raw <- read.fasta(file = "Data/Genetics/clustalo_2026_03_18.fa")
file_to_delete <- which(names(fasta_file_raw) == "JK_902937-1027_39A_JK_003_C04") # Based on alignment, does not cover important 2 SNP region
fasta_files <- fasta_file_raw[-file_to_delete]

# need to remove that one
for (i in 1:length(fasta_files)) {
  sequence <- fasta_files[[i]]
  new_sequence <- sequence[117:655]
  
  fasta_files[[i]] <- new_sequence
  
}

write.fasta(sequences = fasta_files, names = names(fasta_files), file.out ="Data/Genetics/clustalo_2026_03_18_trimmed.fa")

##### Run NCBI Samples plus Palau representative samples trim ####
##### Run Palau only aligned trim ###
fasta_file_raw <- read.fasta(file = "Data/Genetics/aligned_fasta_20260608.fa")
file_to_delete <- which(names(fasta_file_raw) == "FJ536742.1")
fasta_files <- fasta_file_raw[-file_to_delete]


for (i in 1:length(fasta_files)) {
  sequence <- fasta_files[[i]]
  new_sequence <- sequence[90:990]

  fasta_files[[i]] <- new_sequence
  
}

write.fasta(sequences = fasta_files, names = names(fasta_files), file.out ="Data/Genetics/aligned_fasta_20260608_trimmed.fa")

# These files were then sent to iqtree3 to create the phylogenetic tree based on maximum likelihood
# Returns as .treefile in Processed_data
D
