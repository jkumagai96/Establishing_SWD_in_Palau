# Date: March 17th 2026
# Purpose: Truncate FASTA files
# Seagrass Disease in Palau

##### Load Packages ############################################################
library(tidyverse)
library(stringr)
library(seqinr)

##### Create Function ##########################################################
# Load Data 
file_fasta <- read.fasta(file = "G:/My Drive/_PhD/Research/Seagrass in Palau/Genetics/2026_Sequences_18S_SSU/FASTA files/Martin_et_all_plus_2025_cultures_and_2026_all_no45B.txt", seqtype = "DNA")
fasta_final <- file_fasta

# Start for large for loop that goes through each line
for (j in 1:length(fasta_final)) {

  # Select file
  sequence <- fasta_final[[j]]
  
  # Decide whether it needs to be trimmed, only trim the files that start with jk
  letters <- substr(attributes(sequence)$name, 1, 2) # get the first two letters
  if(letters != "JK") {next}
  
  # Find the index to trim the beginning of the seq 
  indexes_n_beg <- which(sequence[1:100] == "n")
  cut_off_beg <- max(indexes_n_beg)
  
  # Find the index to trim the end of the sequence that has two n's consecutively 
  indexes_n_end <- which(sequence[101:length(sequence)] == "n") + 100
  
  for (i in 1:length(indexes_n_end)) {
    v <- indexes_n_end[i]
    v2 <- v + 1
    v2_real <- indexes_n_end[i+1]
    
    if (v2_real == v2) {
      cut_off_end <- v
      break }
  }
  
  new_sequence <- sequence[(cut_off_beg + 1):(cut_off_end - 1)]
  attributes(new_sequence) <- attributes(sequence)
  
  fasta_final[[j]] <- new_sequence
}

write.fasta(sequences = fasta_final, names = names(fasta_final), file.out ="G:/My Drive/_PhD/Research/Seagrass in Palau/Genetics/2026_Sequences_18S_SSU/FASTA files/2026_samples_18S_and_Martin_et_al_truncated.txt")

# Erick didn't like this approach, so we are going to go with sangarseqR package to do the same thing but with the quality scores from the ab1 file 















