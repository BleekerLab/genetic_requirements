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
df = df %>% filter(!sample %in% c("Elite_02", "F1-hab", "LA1777_F1", "LA1777"))
# create a locus/gene column to prepare future filtering using precursor_genes 
colnames(df)[1] = "transcript" 
df_parsed = df %>% mutate(gene = substr(transcript, start = 1, stop = 14)) 

################
# target genes #
################

###########################
# Layout for all heatmaps #
###########################

annotation_cols = as.data.frame(samples)
row.names(annotation_cols) <- samples$sample
annotation_cols$sample <- NULL

my_colour = list(
  condition = c(elite = "#5977ff", F1= "#f74747", wild = "gray31", F2 = "orangered3"),
  pathway= c(MEP = "#82ed82", MVA = "#9e82ed")
)

# MEP/MVA precursors

precursor_genes <- read.delim("figure_7/precursor_genes.tsv", header = T, stringsAsFactors = F)
# Remove (putative) Nudix and IPK genes as this is still unclear in tomato
# precursor_genes = precursor_genes %>% filter(!name %in% c("Nudix", "IPK"))
# filter the scaled counts using the target genes
df_filtered <- inner_join(precursor_genes,df_parsed, by = "gene")
df_filtered = df_filtered[order(df_filtered$pathway),]
# transform into wide format 
df_filtered_wide <- pivot_wider(df_filtered, id_cols = "gene", names_from = "sample", values_from = "est_counts") 

## Step two: importing the sample information and re-ordering it
samples <- read_tsv("Figure_7/samples.tsv", col_names = TRUE)[c("sample", "condition")]
samples$condition <- with(samples, factor(condition, levels = c("elite","F1","wild","F2"), ordered = TRUE))
samples_ordered <- with(samples, samples[order(condition),])

# Scaling
# convert to matrix
mat = as.data.frame(df_filtered_wide[,-1]) 
row.names(mat) = df_filtered_wide$gene

mat_log2_scaled <- log2(mat + 1)
mat_log2_scaled <- mat_log2_scaled[,samples_ordered$sample]

# heatmap

annotation_rows = as.data.frame(precursor_genes[,-1])
row.names(annotation_rows) <- precursor_genes$gene
annotation_rows = annotation_rows[order(annotation_rows$pathway),]


        pheatmap(mat = mat_log2_scaled, 
         scale = "none", 
         cluster_rows = F, 
         cluster_cols = F,
         fontsize = 10,
         cellwidth = 15,
         cellheight = 15,
         annotation_col = annotation_cols, 
         annotation_row = annotation_rows,
         annotation_colors = my_colour,
         gaps_row = 10,
         filename = "Figure_7/precursors.pdf")
         
         #filename = "Figure_7/precursors.pdf")
        

        
#####################
# Terpene Synthases #
#####################

TPS <- read.delim("figure_7/terpene_synthases_zhou2020.tsv", header = T, stringsAsFactors = F)
TPS = TPS %>% filter(structure == "functional") # Remove pseudogenes 
        
# filter the scaled counts using the target genes
df_filtered_TPS <- inner_join(TPS,df_parsed, by = "gene")
df_filtered_TPS = df_filtered_TPS[order(df_filtered_TPS$annotation),] # order from low to high TPS number
# transform into wide format 
df_filtered_TPS_wide <- pivot_wider(df_filtered_TPS, id_cols = "gene", names_from = "sample", values_from = "est_counts") 

## Step two: importing the sample information and re-ordering it
samples <- read_tsv("Figure_7/samples.tsv", col_names = TRUE)[c("sample", "condition")]
samples$condition <- with(samples, factor(condition, levels = c("elite","F1","wild","F2"), ordered = TRUE))
samples_ordered <- with(samples, samples[order(condition),])

# Scaling
# convert to matrix
mat_TPS = as.data.frame(df_filtered_TPS_wide[,-1]) 
row.names(mat_TPS) = df_filtered_TPS_wide$gene

mat_TPS_log2_scaled <- log2(mat_TPS + 1)
# mat_log2_scaled <- mat_log2_scaled[,samples_ordered$sample]

# heatmap

TPS_expressed = as.data.frame(left_join(df_filtered_TPS_wide, TPS, by = "gene"))
row.names(TPS_expressed) <- TPS_expressed$gene
TPS_expressed = TPS_expressed %>% select(order("annotation"))
annotation_rows = TPS_expressed 


pheatmap(mat = mat_TPS_log2_scaled, 
         scale = "none", 
         cluster_rows = F, 
         cluster_cols = F,
         fontsize = 10,
         cellwidth = 15,
         cellheight = 15,
         annotation_row = annotation_rows,
         annotation_col = annotation_cols,
         annotation_colors = my_colour,
         filename = "Figure_7/TPS_heatmap.pdf")

###########################
# Trans-prenyltranferases #
###########################

TPT <- read.delim("figure_7/trans_prenyltransferases_zhou2020.tsv", header = T, stringsAsFactors = F)

# filter the scaled counts using the target genes
df_filtered_TPT <- inner_join(TPT,df_parsed, by = "gene")
df_filtered_TPT = df_filtered_TPT[order(df_filtered_TPT$annotation),] # order from low to high TPT number
# transform into wide format 
df_filtered_TPT_wide <- pivot_wider(df_filtered_TPT, id_cols = "gene", names_from = "sample", values_from = "est_counts") 

## Step two: importing the sample information and re-ordering it
samples <- read_tsv("Figure_7/samples.tsv", col_names = TRUE)[c("sample", "condition")]
samples$condition <- with(samples, factor(condition, levels = c("elite","F1","wild","F2"), ordered = TRUE))
samples_ordered <- with(samples, samples[order(condition),])

# Scaling
# convert to matrix
mat_TPT = as.data.frame(df_filtered_TPT_wide[,-1]) 
row.names(mat_TPT) = df_filtered_TPT_wide$gene

mat_TPT_log2_scaled <- log2(mat_TPT + 1)
# mat_log2_scaled <- mat_log2_scaled[,samples_ordered$sample]

# heatmap

TPT_expressed = as.data.frame(left_join(df_filtered_TPT_wide, TPT, by = "gene"))
row.names(TPT_expressed) <- TPT_expressed$gene
TPT_expressed = TPT_expressed %>% select("annotation")
annotation_rows = TPT_expressed 


pheatmap(mat = mat_TPT_log2_scaled, 
         scale = "none", 
         cluster_rows = F, 
         cluster_cols = F,
         fontsize = 10,
         cellwidth = 15,
         cellheight = 15,
         annotation_row = annotation_rows,
         annotation_col = annotation_cols,
         annotation_colors = my_colour,
         filename = "Figure_7/TPT_heatmap.pdf")

