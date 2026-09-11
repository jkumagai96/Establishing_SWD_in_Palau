# A reproducible phylogenetic tree visualization
# p -> base_plot -> plot_with_points -> final_plot

# Import libraries
library(tidyverse)
library(ggtree)
library(ape)
library(treeio)
library(readxl)

################################################################################
# Import tree and metadata
################################################################################



# Updates May 21st
# Change color to pathogenic_clade with the label not N, T, and P, but rather Nonpathogenic (Blue), Terrestrial (Green), and Pathogenic (Red)
# Add a star (*) or bold the labels of those that have a confirmed status
# Change the name to include Accession number followed by place, such as "AF265335 — United States: Washington, San Juan Islands"
# Remove the T, P, and N labels


#Tree
# tree_file <- "aligned_fasta_20260403_trimmed.fa.contree"
# tree_file <- "tree (1).nwk"
# tree_file <- "aligned_fasta_20260608_trimmed.fa.contree"
tree_file <- "Processed_data/aligned_fasta_20260608_trimmed.fa.treefile"

tree_df <- read.tree(tree_file)

#Metafile
meta_file <- "Data/Metadata_for_Laby_phylogenetic_tree.xlsx"
meta_df <- read_excel(meta_file)

# Build a tip from tree data
tip_labels <- tree_df$tip.label

# Strip accession trailing digits to match meta_df
tip_clean <- sub("\\..*$", "", tip_labels)
tip_map <- data.frame(label = tip_labels, Accession_number = tip_clean, stringsAsFactors = FALSE)

# Join tip labels to meta data
tip_meta <- tip_map %>%
  left_join(meta_df, by = "Accession_number")

# Attach metadata to the ggtree object (so labels and columns are available in the plot data)
p <- ggtree(tree_df, branch.length = "branch.length") %<+% tip_meta

################################################################################
#Add tip point markers
################################################################################

# Set levels
# pathogenic_levels <- c("Pathogenic", "Non pathogenic")

# Add tip labels to base tree (use ID or fallback to label)
# p$data$tip_label_to_plot <- ifelse(is.na(p$data$ID), p$data$label, p$data$ID)
p$data$tip_label_to_plot <- paste0(sub("\\..*$", "", p$data$label), " \u2014 ", case_when(is.na(p$data$Place_country) ~ "Unknown", .default = as.character(p$data$Place_country)))

p$data <- p$data %>%
  mutate(
    Pathogenic_confirmed = replace_na(Pathogenic_confirmed, "Unknown"),
    Pathogenic_confirmed = trimws(Pathogenic_confirmed),
    tip_fontface = case_when(
      isTip & Pathogenic_confirmed == "Pathogenic" ~ "bold",
      isTip & Pathogenic_confirmed == "Non pathogenic" ~ "bold", #"Non pathogenic" ~ "italic"
      isTip & Place_country == "Palau" ~ "bold",
      TRUE ~ "plain"
    ),
    tip_color = case_when(
      isTip & Place_country == "Palau" ~ "goldenrod",
      .default = "black"
    )
  )

tip_color_map <- c("goldenrod" = "goldenrod","black" = "black")


base_plot <- p + geom_tiplab(aes(label = tip_label_to_plot, fontface = tip_fontface, color = tip_color), align = TRUE, offset = 0.001, size = 3, show.legend = FALSE) +
  scale_color_manual(values = tip_color_map)
base_plot
# Prepare tip points data
tip_points <- p$data %>%
  # If you want just a few points (similar to the first graph I sent)
  # just filter down the remaining sites
  # filter(!is.na(Pathogenic_confirmed)) %>% 
  filter(isTip == TRUE) %>% #remove internal sites
  mutate(Pathogenic_confirmed = replace_na(Pathogenic_confirmed, "Unknown")) %>%
  
  # Ensure Pathogenic_confirmed values are consistent strings
  mutate(Pathogenic_confirmed = trimws(Pathogenic_confirmed)) %>% 
  
  # Also change other labels to better descriptors
  mutate(Pathogenic_clade = Pathogenic_clade |>
           replace_values(
             "P" ~ "Pathogenic",
             "N" ~ "Nonpathogenic",
             "T" ~ "Terrestrial"
           ))

  

# Map shapes/colors
shape_map <- c("Pathogenic" = 22, "Non pathogenic" = 24, "Unknown" = 21)
# shape_map <- c("Pathogenic" = 22, "Nonpathogenic" = 24, "Terrestrial" =23,"Unknown" = 21)

# color_map <- c("Pathogenic" = "red", "Non pathogenic" = "blue")
color_map <- c("Pathogenic" = "firebrick", "Non pathogenic" = "blue","Unknown" = "grey") #, "Terrestrial" ="forestgreen"
plot_with_points <- base_plot +
  geom_point(data = tip_points,
             aes(x = x, y = y,  shape = Pathogenic_confirmed, fill = Pathogenic_confirmed), #
             size = 3, 
             color = "black", 
             inherit.aes = FALSE) + #,
  # shape = 21
  scale_shape_manual(name = "Status", values = shape_map, na.value = NA) +
  scale_fill_manual(name = "Status", values = color_map, na.value = NA)
  # scale_fill_discrete(name = "Source") +
  # guides(
  #   shape = guide_legend(override.aes = list(fill = "grey80", color = "black")),
  #   fill = guide_legend(override.aes = list(shape = 21, color = "black", stroke = 0.5))
  # )
plot_with_points

################################################################################
# Add the clade level description
################################################################################
# Helper function to get all tip descendants of a node
# get_tip_descendants <- function(tree, node) {
#   ntips <- length(tree$tip.label)
#   if (node <= ntips) {
#     return(node) 
#   }
#   # Get all descendants recursively
#   edge_matrix <- tree$edge
#   children <- edge_matrix[edge_matrix[,1] == node, 2]
#   tips <- c()
#   for (child in children) {
#     tips <- c(tips, get_tip_descendants(tree, child))
#   }
#   return(tips)
# }
# 
# Ntip <- length(tree_df$tip.label)
# internal_nodes <- (Ntip + 1):(Ntip + tree_df$Nnode)
gdat <- p$data
# 
# clade_list <- lapply(internal_nodes, function(node) {
#   # Get all tip descendants
#   desc_tip_indices <- get_tip_descendants(tree_df, node)
#   
#   if (length(desc_tip_indices) < 3) return(NULL)  # Need at least 3 tips for a clade
#   
#   desc_labels <- tree_df$tip.label[desc_tip_indices]
#   
#   # Lookup metadata for descendant tips
#   desc_meta <- tip_meta %>% filter(label %in% desc_labels)
#   
#   # Get unique Pathogenic_clade values (excluding NAs and "Unknown")
#   cl_vals <- desc_meta %>% 
#     filter(!is.na(Pathogenic_clade) & 
#              Pathogenic_clade != "" & 
#              Pathogenic_clade != "Unknown") %>%
#     pull(Pathogenic_clade) %>%
#     unique()
#   
#   # Require exactly ONE clade letter among all non-NA descendants
#   # (NAs/Unknowns are allowed and don't break the clade)
#   if (length(cl_vals) != 1) return(NULL)
#   
#   # Get y positions from ggtree data
#   y_tips <- gdat$y[gdat$label %in% desc_labels & gdat$isTip]
#   if (length(y_tips) == 0) return(NULL)
#   
#   xnode <- gdat$x[gdat$node == node]
#   if (length(xnode) != 1) return(NULL)
#   
#   # Change this line if you want to relablel your clades
#   clade_label <- cl_vals[1]
#   # clade_label <- switch(cl_vals[1],
#   #                       "P" = "Pathogenic",
#   #                       "N" = "Non-pathogenic", 
#   #                       "T" = "Terrestrial",
#   #                       cl_vals[1])  # fallback to original
#   
#   data.frame(node = node,
#              x = xnode,
#              ymin = min(y_tips),
#              ymax = max(y_tips),
#              clade = clade_label,
#              clade_code = cl_vals[1],
#              n_tips = length(desc_tip_indices),
#              stringsAsFactors = FALSE)
# })
# 
# clade_df <- bind_rows(Filter(Negate(is.null), clade_list))
# 
# # Keep only the most ancestral clade for each clade type
# # Removing nested clades by keep the one with the smallest x
# 
# if (nrow(clade_df) > 0) {
#   clade_df <- clade_df %>%
#     group_by(clade_code) %>%
#     arrange(x) %>%  # Smaller x = more ancestral (closer to root)
#     slice(1) %>%    # Keep only the most ancestral
#     ungroup()
# }
# 
# 
# # Add vertical lines and clade text labels
# x_offset <- max(gdat$x, na.rm = TRUE) * 0.02 # for bars and labels
# y_offset <- 1.5
# 
# 
# if (nrow(clade_df) > 0) {
#   # We'll draw the vertical line at the node x coordinate, and put the clade label slightly to the right
#   final_plot <- plot_with_points +
#     geom_segment(data = clade_df,
#                  aes(x = x - x_offset,
#                      xend = x - x_offset,
#                      y = ymin,
#                      yend = ymax),
#                  color = "black", size = 1.8) +
#     geom_text(data = clade_df,
#               aes(x = x - (x_offset * 2), #  
#                   y = ymax + y_offset, # (ymin + ymax) / 2,
#                   label = clade),
#               hjust = 0, size = 5.5,
#               fontface = "bold")
# } else {
#   final_plot <- plot_with_points
# }
final_plot <- plot_with_points


# Final theme and design decisions
# Change scalar to add space for label text
scalar_x_plot <- 1.3

final_plot <- final_plot +
  # ggtitle("Labyrinthula Phylogenetic Tree") +
  theme(legend.position = "right") +
  xlim(NA, max(gdat$x, na.rm = TRUE) * scalar_x_plot)


print(final_plot)

