#### Get started ####
rm(list=ls())

setwd("/your_env/")  # need to be changed if you have another user 
getwd()

if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("treeio")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtreeExtra")

library(ggtree) 
library(treeio)
library(ggtreeExtra)

merged_df <- read.delim("your_data")
merged_df$probiotic_new <- as.character(merged_df$probiotic_new)
merged_df$symptoms <- as.character(merged_df$symptoms)

##### script for paper figure 5.a and 5.b: probiotic_case id: #####
probiotic_case_id <- merged_df[, c("case_id", "probiotic_new")]
merged_df$new_name <- paste(merged_df$family_id, merged_df$visit_age_mo.x, sep = "_")
rownames(probiotic_case_id) <- merged_df$sample_id
rownames(probiotic_case_id) <- gsub("sub.", "",rownames(probiotic_case_id))

# e coli:
tree <- ape::read.tree("your/raxml/tree/e_coli/RAxML_bestTree.tre")
p <- ggtree(tree) + 
  geom_tiplab(size=2, align=TRUE, linesize=.5) + 
  theme_tree2()

p$data$label <- gsub("sub.", "",p$data$label)

p1 <- gheatmap(p, probiotic_case_id, offset=0, width=0.2, font.size=3, colnames=TRUE,
         colnames_angle=-45, hjust=0) +
  scale_fill_manual(breaks=c("0", "1", "AP Case", "No AP"), 
                    values=c("blueviolet", "darkgoldenrod1", rgb(155/255, 79/255, 54/255), rgb(128/255, 128/255, 128/255)), name="probiotic, case id")
print(p1)
output_path <- "your/output/path/ecoli_tree.pdf"
ggsave(output_path, p1, width=12, height=15)

# lacticasebacillus:
tree <- ape::read.tree("your/raxml/tree/lacticasebacillus/RAxML_bestTree.tre")
p <- ggtree(tree) + 
  geom_tiplab(size=2, align=TRUE, linesize=.5) + 
  theme_tree2()

p$data$label <- gsub("sub.", "",p$data$label)

p2 <- gheatmap(p, probiotic_case_id, offset=0, width=0.2, font.size=3, colnames=TRUE,
         colnames_angle=-45, hjust=0) +
  scale_fill_manual(breaks=c("0", "1", "AP Case", "No AP"), 
                    values=c("blueviolet", "darkgoldenrod1", rgb(155/255, 79/255, 54/255), rgb(128/255, 128/255, 128/255)), name="probiotic, case id")
p2
output_path <- "your/output/path/lacticasebacillus_tree.pdf"
ggsave(output_path, p2, width=12, height=15)

## lacticasebacillus tree with symptomatic stage:
probiotic_case_id <- merged_df[, c("case_id", "probiotic_new", "symptoms")]
merged_df$new_name <- paste(merged_df$family_id, merged_df$visit_age_mo.x, sep = "_")
rownames(probiotic_case_id) <- merged_df$sample_id
rownames(probiotic_case_id) <- gsub("sub.", "",rownames(probiotic_case_id))

tree <- ape::read.tree("your/raxml/tree/lacticasebacillus/RAxML_bestTree.tre")
p <- ggtree(tree) + 
  geom_tiplab(size=2, align=TRUE, linesize=.5) + 
  theme_tree2()

p$data$label <- gsub("sub.", "",p$data$label)

p2 <- gheatmap(
  p,
  probiotic_case_id[, c("probiotic_new", "case_id", "symptoms")],
  offset = 0,
  width = 0.35,
  font.size = 3,
  colnames = TRUE,
  colnames_angle = -45,
  hjust = 0
)

p2 <- p2 +
  scale_fill_manual(
    name = "Metadata",
    values = c(
      "0" = "blueviolet",
      "1" = "darkgoldenrod1",
      "AP Case" = rgb(155/255, 79/255, 54/255),
      "No AP" = rgb(128/255, 128/255, 128/255),
      "Control" = rgb(128/255, 128/255, 128/255),
      "Pre-symptoms" = "pink",
      "Symptomatic" = "#D73027",
      "Resolved" = "#D73027"
    )
  )

p2
output_path <- "your/output/path/lacticasebacillus_tree_symptoms.pdf"
ggsave(output_path, p2, width=12, height=15)