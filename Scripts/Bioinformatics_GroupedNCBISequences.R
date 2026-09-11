# Date: March 27th 2026
# Purpose: Filter for similarity of aligned NCBI labyrinthula fasta sequences
# Seagrass Disease in Palau

##### Load Packages ############################################################
library(tidyverse)
library(ape)
library(igraph)
library(seqinr)

##### Load data ################################################################
aligned_fasta <- read.fasta(file = "G:/My Drive/_PhD/Research/Seagrass in Palau/Genetics/Extracting Labyrinthula 18S from NCBI/Search_1_aligned.fa")

dna_bin <- read.dna(file = "G:/My Drive/_PhD/Research/Seagrass in Palau/Genetics/Extracting Labyrinthula 18S from NCBI/Search_1_aligned.fa", 
                    format = "fasta")

filename <- "G:/My Drive/_PhD/Research/Seagrass in Palau/Genetics/Extracting Labyrinthula 18S from NCBI/muscle_pim.txt"
t <- read.table(filename, header = FALSE)
muscle_names <- t[,2]
t <- t[,-c(1:2)]
muscle_matrix <- as.matrix(t)
rownames(muscle_matrix ) <- muscle_names
colnames(muscle_matrix ) <- muscle_names

##### Calculate pairwise distances #############################################
# Calculate pairwise distances (raw = proportion of differences)
dist_matrix <- dist.dna(dna_bin, model = "raw", pairwise.deletion = TRUE)

# Convert to matrix for easier manipulation
dist_mat <- as.matrix(dist_matrix)

##### Identify sequences to remove (≥99% similar) ##############################
# 99% similar = 1% different = distance of 0.01
# We want distance ≤ 0.01 (excluding self-comparisons where distance = 0)

# Find pairs of sequences that are ≥99% similar
similar_pairs <- which(dist_mat <= 0.01 & dist_mat > 0, arr.ind = TRUE)

dist_mat <- muscle_matrix
similar_pairs <- which(dist_mat >= 99.0 & dist_mat < 100, arr.ind = TRUE)


# View the similar pairs
if(nrow(similar_pairs) > 0) {
  cat("Found", nrow(similar_pairs), "pairwise comparisons ≥99% similar\n")
  
  # Create a summary of similar sequences
  similar_df <- data.frame(
    seq1 = rownames(dist_mat)[similar_pairs[,1]],
    seq2 = rownames(dist_mat)[similar_pairs[,2]],
    distance = dist_mat[similar_pairs],
    percent_identity = (1 - dist_mat[similar_pairs]) * 100
  )
  
  # View some examples
  print(head(similar_df))
  
  ##### Keep longest sequence from each similar group ##########################
  # Create a graph where similar sequences are connected
  g <- graph_from_edgelist(similar_pairs, directed = FALSE)
  
  # Find connected components (groups of similar sequences)
  clusters <- components(g)
  
  cat("\nFound", clusters$no, "groups of similar sequences\n")
  
  # For each cluster, keep the longest sequence
  to_remove <- c()
  removal_info <- list()
  
  for(cluster_id in 1:clusters$no) {
    # Get indices in this cluster
    cluster_members <- which(clusters$membership == cluster_id)
    
    # Get actual lengths of sequences in cluster (excluding gaps)
    lengths <- sapply(cluster_members, function(idx) {
      seq <- aligned_fasta[[idx]]
      sum(seq != "-")  # Count non-gap characters
    })
    
    # Get names
    cluster_names <- names(aligned_fasta)[cluster_members]
    
    # Keep the longest, remove others
    longest_idx <- cluster_members[which.max(lengths)]
    to_remove_from_cluster <- cluster_members[cluster_members != longest_idx]
    
    # Store info about this cluster
    removal_info[[cluster_id]] <- data.frame(
      cluster = cluster_id,
      sequence_name = cluster_names,
      length = lengths,
      kept = cluster_members == longest_idx
    )
    
    to_remove <- c(to_remove, to_remove_from_cluster)
  }
  
  # Combine removal info for reporting
  removal_report <- do.call(rbind, removal_info)
  
  cat("Removing", length(to_remove), "sequences that are ≥99% similar to others\n")
  cat("Keeping", length(aligned_fasta) - length(to_remove), "unique sequences\n\n")
  
  # Show what's being kept vs removed
  print(removal_report)
  
  # Create final filtered set
  aligned_fasta_unique <- aligned_fasta[-to_remove]
  
} else {
  cat("No sequences found with ≥99% similarity\n")
  aligned_fasta_unique <- aligned_fasta
  similar_df <- NULL
  removal_report <- NULL
}

##### Save the final filtered aligned sequences ################################
write.fasta(sequences = aligned_fasta_unique, 
            names = names(aligned_fasta_unique),
            file.out = "G:/My Drive/_PhD/Research/Seagrass in Palau/Genetics/Extracting Labyrinthula 18S from NCBI/Search_1_aligned_unique_muscle.fasta")

# Save reports
if(!is.null(similar_df)) {
  write.csv(similar_df, 
            file = "G:/My Drive/_PhD/Research/Seagrass in Palau/Genetics/Extracting Labyrinthula 18S from NCBI/similar_sequences_pairwise_muscle.csv",
            row.names = FALSE)
  
  write.csv(removal_report,
            file = "G:/My Drive/_PhD/Research/Seagrass in Palau/Genetics/Extracting Labyrinthula 18S from NCBI/removal_report.csv",
            row.names = FALSE)
}

##### Summary statistics #######################################################
cat("After removing ≥99% similar:", length(aligned_fasta_unique), "\n")

# Next step is to make lists of each accession numbers of the paper and of mine and 
# compare where I may be missing five of their numbers and where that happened.

##### Compare to Agnew-Camiener et al. 2025 Accession numbers ##################
library(readxl)
agnew_df <- read_xlsx("G:/My Drive/_PhD/Research/Seagrass in Palau/Genetics/Agnew_Camiener et al. 2025 Accession numbers.xlsx") %>% 
  mutate(agnew = 1)

original_accession <- str_sub(names(aligned_fasta), end = -3)
accesion_n_kept <- str_sub(names(aligned_fasta_unique), end = -3)
df_n_kept <- data.frame(accesion_n_kept, NCBI_search = 1)
colnames(df_n_kept) <- c("Accession_number", "NCBI")

df <- full_join(agnew_df, df_n_kept)

write.csv(df,
          file = "G:/My Drive/_PhD/Research/Seagrass in Palau/Genetics/Extracting Labyrinthula 18S from NCBI/comparision_Agnew_mine_muscle.csv",
          row.names = FALSE)
