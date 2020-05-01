# Figure 7

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
df <- read_tsv("figure_7/abundance_tidy.tsv", col_names = TRUE)
# create a locus/gene column to prepare future filtering using target_genes 
colnames(df)[1] = "transcript" 
df_parsed = df %>% mutate(gene = substr(transcript, start = 1, stop = 14)) 

# target genes 
target_genes <- read.delim("figure_7/targets.tsv", header = T, stringsAsFactors = F)

# filter the scaled counts using the target genes
df_filtered <- inner_join(target_genes,df_parsed, by = "gene")
df_filtered = df_filtered[order(df_filtered$pathway),]

# transform into wide format 
df_filtered_wide <- pivot_wider(df_filtered, id_cols = "gene", names_from = "sample", values_from = "est_counts") 

## Step two: importing the sample information and re-ordering it
samples <- read_tsv("Figure_7/samples.tsv", col_names = TRUE)[c("sample", "condition")]
samples$condition <- with(samples, factor(condition, levels = c("elite","F1","wild","F2")))
samples_ordered <- with(samples, samples[order(condition),])

#############
# scaling
#############

# convert to matrix
mat = as.data.frame(df_filtered_wide[,-1]) 
row.names(mat) = df_filtered_wide$gene

mat_log2_scaled <- log2(mat + 1)
mat_log2_scaled <- mat_log2_scaled[,samples_ordered$sample]

# heatmap

annotation_rows = as.data.frame(target_genes[,-1])
row.names(annotation_rows) <- target_genes$gene
annotation_rows = annotation_rows[order(annotation_rows$pathway),]

annotation_cols = as.data.frame(samples)
row.names(annotation_cols) <- samples$sample
annotation_cols$sample <- NULL

hm = pheatmap(mat = mat_log2_scaled, 
         scale = "none", 
         cluster_rows = F, 
         cluster_cols = F,
         annotation_col = annotation_cols, 
         annotation_row = annotation_rows)

