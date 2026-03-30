##
install.packages("ggdist")
library(ggdist)
install.packages("ggrastr")
library(ggrastr)
install.packages("Cairo", type = "source")


setwd("/your_env/")  # need to be changed if you have another user 
getwd()

# combined
input_base <- "mash_distances_data.tsv"
species <- readLines("list_of_species.txt")
output_base <- "your_output_dir/" # need to be changed if you have another user

mash_distance_files <- paste0(input_base, species, "/mash_distances_data.tsv")
output_plot_files <- paste0(output_base, species, ".pdf")
output_plot_files_heatmaps <- paste0(output_base, species, ".combined_heatmap.pdf")

extract_family_id <- function(sample_name) {
  sub("^sub\\.(\\d+)\\..*$", "\\1", sample_name)
}

extract_sample_id <- function(sample_name) {
  sub("^sub\\.(\\d+\\.\\d+)_.*$", "\\1", sample_name)
}

extract_bin <- function(sample_name) {
  sub("^.*\\.(\\d+)$", "\\1", sample_name)
}


###########################################################
# Main execution code
###########################################################

# Create an empty dataframe to collect all results
results_df <- data.frame(
  species = character(),
  n_ap_samples = numeric(),
  n_noap_samples = numeric(),
  n_related = numeric(),
  n_unrelated = numeric(),
  raw_p_value_AP_vs_NoAP = numeric(),
  adjusted_p_value_AP_vs_NoAP = numeric(),
  raw_p_value_related_vs_unrelated = numeric(),
  adjusted_p_value_related_vs_unrelated = numeric(),
  stringsAsFactors = FALSE
)

# Create a list to store all plots
all_plots <- list()

# Process each species
for (i in seq_along(mash_distance_files)) {
  # Run the improved function and get results
  species_result <- box_plot_ap_control_relation(
    mash_distance_file = mash_distance_files[[i]],
    output_plot_file = output_plot_files[[i]],
    title = species[i],
    species_meta_data = species_dataframes_combined[[species[i]]]
  )
  
  # Add the results to our collection dataframe
  results_df <- rbind(results_df, species_result)
  
  # Store the plot for later use
  all_plots[[species[i]]] <- attr(species_result, "plot")
  
  # Report progress 
  
  cat("\nProcessed", species[i], "\n")
}

print(results_df)

write_csv(results_df, file = "your_output_file.csv")

## density plots:

density_ap_control <- function(mash_distance_file, output_plot_file, title, species_meta_data) {
  print(title)
  
  # Read and process mash distance data
  mash_dist <- read.delim(mash_distance_file)
  
  mash_dist$genome1 <- mash_dist[[1]]  
  mash_dist$genome1 <- gsub(".fa", "", mash_dist$genome1)
  mash_dist$genome2 <- mash_dist[[2]]  
  mash_dist$genome2 <- gsub(".fa", "", mash_dist$genome2)
  
  data_mash <- mash_dist[, c("genome1", "genome2", "X0")]
  data_mash <- data_mash[data_mash$genome1 != data_mash$genome2, ]
  
  # Get unique pairs (upper triangle)
  data_mash_unique <- data_mash %>%
    mutate(pair_id = pmin(genome1, genome2), pair_id2 = pmax(genome1, genome2)) %>%
    distinct(pair_id, pair_id2, .keep_all = TRUE) %>%
    select(-pair_id, -pair_id2)
  
  # Add relation information
  data_mash_relation <- data_mash_unique %>%
    mutate(row_family_id = extract_family_id(as.character(genome1)),
           col_family_id = extract_family_id(as.character(genome2)),
           relation = if_else(row_family_id == col_family_id, "related", "unrelated"))
  
  # Filter and join with metadata for AP analysis
  filtered_df <- subset(data_mash_relation, relation == "related")
  metadata <- inner_join(species_meta_data, filtered_df, by = c("BinName" = "genome1"))
  filtered_df_and_meta <- metadata[, c("X0", "case_id")]
  
  # Extract data for statistical tests
  no_ap_values <- filtered_df_and_meta %>% filter(case_id == "No AP") %>% pull(X0)
  ap_values <- filtered_df_and_meta %>% filter(case_id == "AP Case") %>% pull(X0)
  related <- data_mash_relation %>% filter(relation == "related") %>% pull(X0)
  unrelated <- data_mash_relation %>% filter(relation == "unrelated") %>% pull(X0)
  
  print("no ap, ap:")
  print(length(no_ap_values))
  print(length(ap_values))
  print("related, unrelated:")
  print(length(related))
  print(length(unrelated))
  
  # Collect p-values for  Bonferroni correction
  p_values <- c()
  test_names <- c()
  
  # AP vs No AP test
  if (length(ap_values) > 1 && length(no_ap_values) > 1) {
    wilcox_test <- wilcox.test(ap_values, no_ap_values, exact = FALSE)
    p_values <- c(p_values, wilcox_test$p.value)
    test_names <- c(test_names, "AP_vs_NoAP")
  } else {
    p_values <- c(p_values, NA)
    test_names <- c(test_names, "AP_vs_NoAP")
  }
  
  # Related vs Unrelated test
  if (length(unrelated) > 1 && length(related) > 1) {
    wilcox_test2 <- wilcox.test(unrelated, related, exact = FALSE)
    p_values <- c(p_values, wilcox_test2$p.value)
    test_names <- c(test_names, "Related_vs_Unrelated")
  } else {
    p_values <- c(p_values, NA)
    test_names <- c(test_names, "Related_vs_Unrelated")
  }
  
  # Apply Bonferroni correction across all tests
  adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
  names(adjusted_p_values) <- test_names
  
  # Print adjusted p-values
  for (i in 1:length(test_names)) {
    if (!is.na(adjusted_p_values[i])) {
      print(paste(test_names[i], "adjusted p-value:", signif(adjusted_p_values[i], 5)))
    } else {
      print(paste(test_names[i], "adjusted p-value: NA (insufficient data)"))
    }
  }
  
  # Set factor levels for plot ordering
  data_mash_relation$relation <- factor(data_mash_relation$relation, levels = c("unrelated", "related"))
  filtered_df_and_meta$case_id <- factor(filtered_df_and_meta$case_id, levels = c("No AP", "AP Case"))
  
  
  # Plot density
  g <- ggplot(filtered_df_and_meta, aes(x = X0, fill = case_id)) +
    geom_density(alpha = 0.5) +
    labs(
      title = paste(title),
      fill = "Relation between compared samples"
    ) +
    scale_fill_manual(values = c(
      "No AP" = rgb(128/255, 128/255, 128/255),
      "AP Case" = rgb(155/255, 79/255, 54/255)
    )) +
    theme_bw() +
    coord_flip()
  
  # 
  # Save the plot
  output_plot_file_with_relation <- sub("(\\.[a-z]+)$", "_density_case_id\\1", output_plot_file)
  print(g)
  ggsave(output_plot_file_with_relation, g, width = 12, height = 12)
}


for (i in seq_along(mash_distance_files)) {
  density_ap_control(mash_distance_files[[i]], output_plot_files[[i]], species[i], species_dataframes_combined[[species[i]]])
}