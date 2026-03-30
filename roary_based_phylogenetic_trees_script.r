#########
tree_object <- ape::read.tree("your_roary_output_tree/accessory_binary_genes.fa.newick")

extract_subtrees_df <- function(tree) {
  
  # Get the number of tips and internal nodes
  num_tips <- length(tree$tip.label)
  total_nodes <- num_tips + tree$Nnode
  
  # Initialize an empty data frame
  subtrees_df <- data.frame(Main_Node = integer(),
                            Tip_Labels = character(),
                            Num_Tips = integer(),
                            stringsAsFactors = FALSE)
  
  # Iterate over each internal node to extract subtrees
  for (node in (num_tips + 1):total_nodes) {
    
    # Extract the subtree rooted at the current node
    subtree <- extract.clade(tree, node)
    
    # Get the tip labels of the extracted subtree
    tip_labels <- subtree$tip.label
    
    # Create a new row for the data frame
    new_row <- data.frame(Main_Node = node,
                          Tip_Labels = paste(tip_labels, collapse = ", "),
                          Num_Tips = length(tip_labels),
                          stringsAsFactors = FALSE)
    
    # Append the new row to the data frame
    subtrees_df <- rbind(subtrees_df, new_row)
  }
  
  return(subtrees_df)
}

result <- extract_subtrees_df(tree_object)

transform_to_long_df <- function(subtrees_df) {
  
  # Initialize an empty data frame for the long format
  long_df <- data.frame(Main_Node = integer(),
                        Tip_Label = character(),
                        num_tips = character(), ## added now
                        stringsAsFactors = FALSE)
  
  # Iterate over each row of the input data frame
  for (i in seq_len(nrow(subtrees_df))) {
    # Extract the main node and tip labels for the current subtree
    main_node <- subtrees_df$Main_Node[i]
    tip_labels <- unlist(strsplit(subtrees_df$Tip_Labels[i], ", "))
    subtrees_df$Num_Tips[i] <- as.character(subtrees_df$Num_Tips[i]) ## added now
    num_tips <- unlist(strsplit(subtrees_df$Num_Tips[i], ", ")) ## added now
    
    # Create a data frame for the current subtree's tip labels
    temp_df <- data.frame(Main_Node = rep(main_node, length(tip_labels)),
                          Tip_Label = tip_labels,
                          num_tips = num_tips, ## added now
                          stringsAsFactors = FALSE)
    
    # Append the temporary data frame to the long format data frame
    long_df <- rbind(long_df, temp_df)
  }
  
  # merge with meta-data: 
  sample_names <- "your_sample_ids"  # Replace with actual sample IDs
  sample_names <- str_to_lower(sample_names)
  case_id_df <- data.frame(sample_name = sample_names, case_id = your_meta_data$case_id, sample_id = your_meta_data$family_id)
  extracted_samples <- sub("^(.*?)_.*$", "\\1", long_df$Tip_Label)
  long_df <- long_df %>%
    mutate(sample_extracted = extracted_samples) %>%
    left_join(case_id_df, by = c("sample_extracted" = "sample_name"))
  
  return(long_df)
}

long_df <- transform_to_long_df(result)

# calculate statist in entire tree:
df_tree <- long_df %>%
  distinct(Tip_Label, .keep_all = TRUE)

# function to perform the hypothesis test for each subtree 
main_counts <- table(df_tree$case_id)
total_main <- sum(main_counts)

prop_AP <- 0.5
prop_No_AP <- 0.5

# a function to perform the Chi-square goodness-of-fit test for each subtree with FDR correction
chi_square_goodness_of_fit2 <- function(df, main_prop_AP, main_prop_No_AP) {
  results <- data.frame(Main_Node = integer(), AP = integer(), No_AP = integer(),
                        Chi_Square = numeric(), P_Value = numeric(), 
                        FDR_Adjusted_P_Value = numeric(), infant_based = logical(), stringsAsFactors = FALSE)
  
  unique_nodes <- unique(df$Main_Node)
  
  p_values <- numeric(length(unique_nodes))
  
  unique_families <- long_df %>%
    distinct(sample_id, .keep_all = TRUE)
  
  for (i in seq_along(unique_nodes)) {
    node <- unique_nodes[i]
    # Subset data for each subtree
    subtree_data <- subset(df, Main_Node == node)
    
    # Only include subtrees with at least 4 samples
    if (nrow(subtree_data) < 2) {
      next  # Skip this subtree if it has fewer than 4 samples
    }
    
    # Observed counts in the subtree
    subtree_counts <- table(subtree_data$case_id)
    observed_AP <- ifelse("AP Case" %in% names(subtree_counts), subtree_counts["AP Case"], 0)
    observed_No_AP <- ifelse("No AP" %in% names(subtree_counts), subtree_counts["No AP"], 0)
    total_subtree <- observed_AP + observed_No_AP
    
    # Expected counts based on the main tree distribution
    expected_AP <- total_subtree * main_prop_AP
    expected_No_AP <- total_subtree * main_prop_No_AP
    
    # Create observed and expected vectors
    observed <- c(observed_AP, observed_No_AP)
    expected <- c(expected_AP, expected_No_AP)
    
    # Calculate Chi-square statistic manually
    chi_square_stat <- sum((observed - expected)^2 / expected)
    
    # Calculate p-value - CDF 
    p_value <- pchisq(chi_square_stat, df = 1, lower.tail = FALSE)
    
    # Store the p-value for FDR correction
    p_values[i] <- p_value
    
    # Get family IDs using merge instead of match
    sample_id <- subtree_data %>%
      left_join(unique_families, by = "sample_id") %>%
      pull(sample_id)
    
    # Check infant_based condition
    infant_based <- if(length(sample_id) == 0) {
      NA
    } else {
      length(unique(sample_id)) == 1 &&  # All same family
        !any(is.na(sample_id))             # No missing values
    }
    
    # Append results to the data frame without FDR adjustment (to be added later)
    results <- rbind(results, data.frame(Main_Node = node, 
                                         AP = observed_AP, 
                                         No_AP = observed_No_AP, 
                                         Chi_Square = chi_square_stat, 
                                         P_Value = p_value,
                                         infant_based = infant_based))
  }
  results$FDR_Adjusted_P_Value <- p.adjust(results$P_Value, method = "fdr")
  return(results)
}

# Apply the function
test_results <- chi_square_goodness_of_fit2(long_df, prop_AP, prop_No_AP)

## add pie charts to nodes
p <- ggtree(ape_tree_object)

test_results$node <- test_results$Main_Node

my_colors2 <- c("AP" = rgb(155/255, 79/255, 54/255), "No_AP" =  rgb(128/255, 128/255, 128/255))
is_feature <- test_results$infant_based
### circler
pies <- nodepie(test_results, cols = c("AP", "No_AP"), color = my_colors2, alpha=.8)
pies_with_shapes <- lapply(seq_along(pies), function(i) {
  pie <- pies[[i]]
  if (is_feature[i]) {
    pie <- pie + 
      geom_point(aes(x = 0, y = 0), shape = 21, size = 3, fill = "white", color = "black")
  }
  return(pie)
})
df <- tibble::tibble(node=as.numeric(test_results$node), pies=pies_with_shapes)
p1 <- ggtree(ape_tree_object, layout = "circular") %>% as.treedata() %>% as_tibble() 

meta_and_tree <-  p1 %>% inner_join(your_meta_data, by=c("label" = "BinName"))
tree_meta_x <- ggtree(ape_tree_object, layout = "circular") %<+% meta_and_tree
tree_meta_x$data <- tree_meta_x$data %>%
  mutate(label = gsub("sub.", "", label)) %>% 
  mutate(label = sub(".gff.*", "", label))
p2 <- tree_meta_x %<+% df

p_circle <- p2 + geom_plot(data=td_filter(!isTip), 
                           mapping=aes(x=x, y=y, label=pies), 
                           vp.width=0.05, vp.height=0.05, 
                           hjust=0.5, vjust=0.5) +
  geom_tiplab(size = 2) +
  geom_tippoint(size = 1, aes(color = case_id)) +
  scale_color_manual(values = c("AP Case" = rgb(155/255, 79/255, 54/255), "No AP" =  rgb(128/255, 128/255, 128/255))) +
  xlim(-1, NA)+
  labs(title = title)

p_circle
output_plot_file_p2 <- sub("(\\.[a-z]+)$", "_pangenome_roary_tree_new\\1", output_plot_file)
ggsave(output_plot_file_p2, p_circle, width=12, height=12)

# only significant subtrees:
sig_results <- test_results[test_results$FDR_Adjusted_P_Value <= 0.05,]
is_feature <- sig_results$infant_based
pies <- nodepie(sig_results, cols = c("AP", "No_AP"), color = my_colors2, alpha=.8)
pies_with_shapes <- lapply(seq_along(pies), function(i) {
  pie <- pies[[i]]
  if (is_feature[i]) {
    pie <- pie + 
      geom_point(aes(x = 0, y = 0), shape = 21, size = 3, fill = "white", color = "black")
  }
  return(pie)
})
df <- tibble::tibble(node=as.numeric(sig_results$node), pies=pies_with_shapes)
p1 <- ggtree(ape_tree_object, layout = "circular") %>% as.treedata() %>% as_tibble() 

meta_and_tree <-  p1 %>% inner_join(your_meta_data, by=c("label" = "BinName"))
tree_meta_x <- ggtree(ape_tree_object, layout = "circular") %<+% meta_and_tree
tree_meta_x$data <- tree_meta_x$data %>%
  mutate(label = gsub("sub.", "", label)) %>% 
  mutate(label = sub(".gff.*", "", label))
p2 <- tree_meta_x %<+% df
p_circle2 <- p2 + geom_plot(data=td_filter(!isTip), 
                            mapping=aes(x=x, y=y, label=pies), 
                            vp.width=0.03, vp.height=0.03, 
                            hjust=0.5, vjust=0.5) +
  geom_tiplab(size = 2) +
  geom_tippoint(size = 1, aes(color = case_id)) +
  scale_color_manual(values = c("AP Case" = rgb(155/255, 79/255, 54/255), "No AP" =  rgb(128/255, 128/255, 128/255))) +
  xlim(-1, NA)+
  labs(title = "E. coli")


p_circle2
output_plot_file_p3 <- sub("(\\.[a-z]+)$", "_pangenome_sig_roary_tree_new\\1", output_plot_file)

##### test on e coli: ####
roary_e_coli <- "your_roary_results/gene_presence_absence.csv"
e_coli_roary <- read.csv(roary_e_coli, sep = ",", header = TRUE)
e_coli_genes <- e_coli_roary[, c(15:ncol(e_coli_roary))]
rownames(e_coli_genes) <- e_coli_roary$Gene
e_coli_genes <- e_coli_genes %>%
  mutate_all(~ ifelse(. == "", 0, 1))

df_filtered <- e_coli_genes[rowSums(e_coli_genes) > 1, ]
df_filtered <- df_filtered[rowSums(df_filtered) != ncol(df_filtered), ]


colnames(df_filtered) <- sub("bin.*", "bin", colnames(df_filtered))
e_coli_df_combined$BinName <- sub("bin.*", "bin", e_coli_df_combined$BinName)

annotation_cols <- data.frame(
  symptoms = as.factor(e_coli_df_combined$symptoms),
  AP = e_coli_df_combined$case_id,
  row.names = e_coli_df_combined$BinName
)

print(all(colnames(df_filtered) %in% rownames(annotation_cols)))

annotation_colors <- list(
  AP = c("AP Case" = rgb(155/255, 79/255, 54/255), "No AP" = rgb(128/255, 128/255, 128/255)),
  symptoms = c(
    "Control" = rgb(128/255, 128/255, 128/255),
    "Pre-symptoms" = rgb(249/255, 222/255, 157/255),
    "Symptomatic" = rgb(221/255, 113/255, 77/255),
    "Resolved" = rgb(233/255, 164/255, 100/255)
  )
)

annotation_cols$AP <- factor(annotation_cols$AP, levels = c("AP Case", "No AP"))
colors <- colorRampPalette(c("#fff7bc", "red"))(n = 100)
p <- pheatmap(df_filtered,
              color=colors,  width = 40, height = 20,
              annotation_col = annotation_cols,
              annotation_colors = annotation_colors,
              fontsize = 10, border_color = NA,
              show_colnames = F, show_rownames = F,
              main = "Escherichia_coli")

p
row_dend <- p$tree_row
col_dend <- p$tree_col


main_counts <- table(your_meta_data$case_id)
total_main <- sum(main_counts)
prop_AP <- main_counts["AP Case"] / total_main
prop_No_AP <- main_counts["No AP"] / total_main


# Extract the row dendrogram (or column dendrogram, depending on your needs)
col_dend <- p$tree_col
dend <- as.dendrogram(p$tree_col)

# Function to extract subtrees and their leaves
extract_subtrees <- function(dend) {
  subtrees <- list()
  
  traverse <- function(node, id) {
    if (is.leaf(node)) {
      return(list(id = id, leaves = labels(node)))
    } else {
      left <- traverse(node[[1]], paste0(id, "1"))
      right <- traverse(node[[2]], paste0(id, "2"))
      subtrees[[id]] <<- c(left$leaves, right$leaves)
      return(list(id = id, leaves = c(left$leaves, right$leaves)))
    }
  }
  
  traverse(dend, "1")
  return(subtrees)
}

subtrees <- extract_subtrees(dend)

# Modify the chi_square_goodness_of_fit2 function
chi_square_goodness_of_fit2 <- function(df, subtrees, main_prop_AP, main_prop_No_AP) {
  results <- data.frame(Main_Node = character(), AP = integer(), No_AP = integer(),
                        Chi_Square = numeric(), P_Value = numeric(), 
                        FDR_Adjusted_P_Value = numeric(), stringsAsFactors = FALSE)
  
  p_values <- numeric(length(subtrees))
  
  for (i in seq_along(subtrees)) {
    node <- names(subtrees)[i]
    subtree_leaves <- subtrees[[i]]
    
    # Subset data for each subtree
    subtree_data <- df[df$BinName %in% subtree_leaves, ]
    
    if (nrow(subtree_data) < 5) {
      next
    }
    
    # Rest of the function remains the same
    subtree_counts <- table(subtree_data$case_id)
    observed_AP <- ifelse("AP Case" %in% names(subtree_counts), subtree_counts["AP Case"], 0)
    observed_No_AP <- ifelse("No AP" %in% names(subtree_counts), subtree_counts["No AP"], 0)
    total_subtree <- observed_AP + observed_No_AP
    
    expected_AP <- total_subtree * main_prop_AP
    expected_No_AP <- total_subtree * main_prop_No_AP
    
    observed <- c(observed_AP, observed_No_AP)
    expected <- c(expected_AP, expected_No_AP)
    
    chi_square_stat <- sum((observed - expected)^2 / expected)
    p_value <- pchisq(chi_square_stat, df = 1, lower.tail = FALSE)
    
    p_values[i] <- p_value
    
    results <- rbind(results, data.frame(Main_Node = node, 
                                         AP = observed_AP, 
                                         No_AP = observed_No_AP, 
                                         Chi_Square = chi_square_stat, 
                                         P_Value = p_value))
  }
  
  # Add FDR adjustment
  results$FDR_Adjusted_P_Value <- p.adjust(results$P_Value, method = "fdr")
  
  return(results)
}

# Run the analysis
results <- chi_square_goodness_of_fit2(e_coli_df_combined, subtrees, prop_AP, prop_No_AP)

###### circular trees #####
library(ggtree)
library(ggplot2)
library(dendextend)

# Define the function for node color determination
get_node_color <- function(results) {
  node_colors <- setNames(rep("black", nrow(results)), results$Main_Node)
  
  for (i in 1:nrow(results)) {
    if (results$AP[i] > results$No_AP[i]) {
      node_colors[results$Main_Node[i]] <- "red"  # More AP
    } else {
      node_colors[results$Main_Node[i]] <- "blue"  # More No AP
    }
  }
  
  return(node_colors)
}

# Get node colors
node_colors <- get_node_color(results)

# Function to color dendrogram nodes
color_dendrogram <- function(dend, node_colors) {
  traverse <- function(node, id) {
    if (!is.leaf(node) && id %in% names(node_colors)) {
      attr(node, "edgePar") <- list(col = node_colors[[id]], lwd = 2)
    }
    
    if (!is.leaf(node)) {
      node[[1]] <- traverse(node[[1]], paste0(id, "1"))
      node[[2]] <- traverse(node[[2]], paste0(id, "2"))
    }
    
    return(node)
  }
  
  return(traverse(dend, "1"))
}

# Apply colors to the dendrogram
colored_dend <- color_dendrogram(dend, node_colors)

# Convert the colored dendrogram into hclust and then phylo format
hc <- as.hclust(colored_dend)
dend_tree <- as.phylo(hc)

# Create a circular dendrogram using ggtree
p <- ggtree(dend_tree, layout = "circular", size = 1) +
  theme_void() +  # Remove axes and text labels
  theme(plot.background = element_rect(fill = "white"),  # Set background to white
        panel.grid.major = element_line(color = "gray", size = 0.5)) +  # Add grid lines
  ggtitle("E.coli Circular Dendrogram") +  # Add title
  scale_color_identity()  # Use the node colors directly

# Display the plot
print(p)    