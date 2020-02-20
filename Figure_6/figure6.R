# Figure 6

#########
# Library
#########
library("checkpoint")
#checkpoint("2020-01-01")

library(tidyverse)
library(pheatmap)


#############
# data import
#############


## Step 1: filtering scaled counts using target genes

# scaled counts data preparation
df <- read_tsv("Figure_6/abundance_tidy.tsv", col_names = TRUE)
# create a locus/gene column to prepare future filtering using target_genes 
colnames(df)[1] = "transcript" 
df_parsed = df %>% mutate(gene = substr(transcript, start = 1, stop = 14)) 

# target genes 
target_genes <- read.delim("Figure_6/targets.tsv", header = T, stringsAsFactors = F)

# filter the scaled counts using the target genes
df_filtered <- inner_join(df_parsed, target_genes, by = "gene")

# transform into wide format 
df_filtered_wide <- pivot_wider(df_filtered, id_cols = "gene", names_from = "sample", values_from = "est_counts") 

## Step two: importing the sample information
samples <- read_tsv("Figure_6/samples.tsv", col_names = TRUE)[c("sample", "condition")]

#############
# scaling
#############

# convert to matrix
mat = as.data.frame(df_filtered_wide[,-1]) 
row.names(mat) = df_filtered_wide$gene

# center data
# matrix has to be transposed for scaling ...
# .. and transposed again for heatmap visualisation
# Check success with rowMeans(mat_scaled) which should return values close to 0
# scale = F to retain the real gene variation
mat_scaled <- t(scale(x = t(mat), center = T, scale = F)) 

# heatmap

annotation_rows = as.data.frame(target_genes[,-1])
row.names(annotation_rows) <- target_genes$gene

annotation_cols = as.data.frame(samples)
row.names(annotation_cols) <- samples$sample
annotation_cols$sample <- NULL

pheatmap(mat = mat_scaled, 
         scale = "none", 
         cluster_rows = F, 
         cluster_cols = F,
         annotation_col = annotation_cols, 
         annotation_row = annotation_rows)

