# Date: March 18th 2026
# Purpose: Truncate ab1 files based on quality scores from ab1 files
# Seagrass Disease in Palau

##### Load Packages ############################################################
library(tidyverse)
library(seqinr)
library(sangerseqR)

##### Create Function ##########################################################
trim_sanger <- function(ab1_file, quality_scores, min_quality = 20, window_size = 10) {
  # Get the primary sequence and quality scores
  sequence <- primarySeq(ab1_file)
  
  # If quality vector is shorter, stop!
  if (length(quality_scores) < length(sequence)) {
    print("Problem! Quality scores are less than sequence length")
    break
  }
  
  # Find trimming positions
  # Trim from start
  start_pos <- 1
  for (i in 1:(length(quality_scores) - window_size)) {
    window_mean <- mean(quality_scores[i:(i + window_size - 1)])
    if (window_mean >= min_quality) {
      start_pos <- i
      break
    }
  }
  
  # Trim from end
  end_pos <- length(quality_scores)
  for (i in length(quality_scores):(window_size + 1)) {
    window_mean <- mean(quality_scores[(i - window_size + 1):i])
    if (window_mean >= min_quality) {
      end_pos <- i
      break
    }
  }
  
  # Return trimmed sequence
  list(
    trimmed_seq = subseq(sequence, start_pos, end_pos),
    start = start_pos,
    end = end_pos,
    original_length = length(sequence),
    trimmed_length = end_pos - start_pos + 1
  )
}

##### Example of the function ##################################################
# Get sequence data
filename <- "G:/My Drive/_PhD/Research/Seagrass in Palau/Genetics/2026_Sequences_18S_SSU/Sangar_samples_Palau_and_M2B2/JK_902937-1001_40_JK_003_A01.ab1"
ab1_file <- readsangerseq(filename)

# Obtain quality scores
abif_data <- sangerseqR::read.abif(filename)
quality_scores <- abif_data@data$PCON.1

# Run the trimming function
result <- trim_sanger(ab1_file, quality_scores, 40, 20)

# View results
result

# Visualize the trimming 
library(ggplot2)

# Prepare data frame
quality_df <- data.frame(
  position = 1:length(quality_scores),
  quality = quality_scores,
  region = ifelse(1:length(quality_scores) < result$start, "Trimmed (5')",
                  ifelse(1:length(quality_scores) > result$end, "Trimmed (3')", "Kept"))
)

# Create plot
ggplot(quality_df, aes(x = position, y = quality)) +
  geom_ribbon(aes(ymin = 0, ymax = quality, fill = region), alpha = 0.3) +
  geom_line(size = 0.8) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = result$start, linetype = "solid", color = "darkgreen", size = 1) +
  geom_vline(xintercept = result$end, linetype = "solid", color = "darkgreen", size = 1) +
  scale_fill_manual(values = c("Kept" = "green", 
                               "Trimmed (5')" = "red", 
                               "Trimmed (3')" = "red")) +
  labs(title = "Sanger Sequence Quality and Trimming",
       subtitle = sprintf("Original: %d bp | Kept: %d bp (%.1f%%) | Trimmed: %d bp",
                          result$original_length, 
                          result$trimmed_length,
                          result$trimmed_length/result$original_length*100,
                          result$original_length - result$trimmed_length),
       x = "Base Position",
       y = "Phred Quality Score",
       fill = "Region") +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11))

###### Now run loop over the file directory and export as a fasta file #########
##### Batch Processing Function #####
process_all_ab1_files <- function(directory, 
                                  output_fasta = "trimmed_sequences.fasta",
                                  min_quality = 20, 
                                  window_size = 10,
                                  generate_report = TRUE) {
  
  # Get all AB1 files in the directory
  ab1_files <- list.files(directory, pattern = "\\.ab1$", 
                          full.names = TRUE, ignore.case = TRUE)
  
  if (length(ab1_files) == 0) {
    stop("No AB1 files found in the specified directory!")
  }
  
  cat(sprintf("Found %d AB1 files to process\n", length(ab1_files)))
  
  # Initialize storage
  fasta_sequences <- c()
  summary_stats <- data.frame(
    filename = character(),
    original_length = numeric(),
    trimmed_length = numeric(),
    percent_kept = numeric(),
    trim_5prime = numeric(),
    trim_3prime = numeric(),
    mean_quality = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process each file
  for (i in seq_along(ab1_files)) {
    file_path <- ab1_files[i]
    file_name <- basename(file_path)
    
    cat(sprintf("[%d/%d] Processing: %s\n", i, length(ab1_files), file_name))
    
    tryCatch({
      # Read AB1 file
      ab1_file <- readsangerseq(file_path)
      abif_data <- sangerseqR::read.abif(file_path)
      
      # Get quality scores
      quality_scores <- abif_data@data$PCON.1
      
      # Trim the sequence
      result <- trim_sanger(ab1_file, quality_scores, min_quality, window_size)
      
      if (result$trimmed_length > 900) {
        print("Length trimmed is too high on file -- skipped", file_name)
        next
      }
      
      # Create FASTA header (remove .ab1 extension and any spaces)
      seq_name <- gsub("\\.ab1$", "", file_name, ignore.case = TRUE)
      seq_name <- gsub(" ", "_", seq_name)
      
      # Add to FASTA output
      fasta_sequences <- c(fasta_sequences,
                           paste0(">", seq_name),
                           as.character(result$trimmed_seq))
      
      # Store summary statistics
      summary_stats <- rbind(summary_stats, data.frame(
        filename = file_name,
        original_length = result$original_length,
        trimmed_length = result$trimmed_length,
        percent_kept = round(result$trimmed_length / result$original_length * 100, 2),
        trim_5prime = result$start - 1,
        trim_3prime = result$original_length - result$end,
        mean_quality = round(mean(quality_scores, na.rm = TRUE), 2),
        stringsAsFactors = FALSE
      ))
      
      cat(sprintf("  ✓ Original: %d bp | Trimmed: %d bp (%.1f%% kept)\n",
                  result$original_length, result$trimmed_length,
                  result$trimmed_length / result$original_length * 100))
      
    }, error = function(e) {
      cat(sprintf("  ✗ ERROR processing %s: %s\n", file_name, e$message))
      
      # Add failed file to summary with NA values
      summary_stats <<- rbind(summary_stats, data.frame(
        filename = file_name,
        original_length = NA,
        trimmed_length = NA,
        percent_kept = NA,
        trim_5prime = NA,
        trim_3prime = NA,
        mean_quality = NA,
        stringsAsFactors = FALSE
      ))
    })
  }
  
  # Write FASTA file
  writeLines(fasta_sequences, output_fasta)
  cat(sprintf("\n✓ FASTA file written to: %s\n", output_fasta))
  
  # Write summary report
  if (generate_report) {
    report_file <- gsub("\\.fasta$", "_summary.csv", output_fasta)
    write.csv(summary_stats, report_file, row.names = FALSE)
    cat(sprintf("✓ Summary report written to: %s\n", report_file))
  }
  
  # Print overall summary
  cat("\n=== SUMMARY ===\n")
  cat(sprintf("Total files processed: %d\n", nrow(summary_stats)))
  cat(sprintf("Successfully trimmed: %d\n", sum(!is.na(summary_stats$trimmed_length))))
  cat(sprintf("Failed: %d\n", sum(is.na(summary_stats$trimmed_length))))
  cat(sprintf("Mean original length: %.1f bp\n", 
              mean(summary_stats$original_length, na.rm = TRUE)))
  cat(sprintf("Mean trimmed length: %.1f bp\n", 
              mean(summary_stats$trimmed_length, na.rm = TRUE)))
  cat(sprintf("Mean percent kept: %.1f%%\n", 
              mean(summary_stats$percent_kept, na.rm = TRUE)))
  
  return(summary_stats)
}

##### RUN THE BATCH PROCESSING #####

# Set your directory
ab1_directory <- "G:/My Drive/_PhD/Research/Seagrass in Palau/Genetics/2026_Sequences_18S_SSU/Sangar_samples_Palau_and_M2B2"

# Process all files
results <- process_all_ab1_files(
  directory = ab1_directory,
  output_fasta = "G:/My Drive/_PhD/Research/Seagrass in Palau/Genetics/2026_Sequences_18S_SSU/Trimming_2026/trimmed_sequences.fasta",
  min_quality = 40,
  window_size = 20,
  generate_report = TRUE
)

# View the summary
View(results)
