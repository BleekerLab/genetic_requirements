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
df_parsed = df %>% mutate(target_id = substr(target_id, start = 1, stop = 14))


## Step two: importing the sample information and re-ordering it
samples <- read_tsv("Figure_7/samples.tsv", col_names = TRUE )[c("sample", "condition")]
samples$condition <- with(samples, factor(condition, levels = c("elite","F1","wild","F2"))) 
samples = samples %>% filter(!sample %in% c("Elite_02", "F1-hab", "LA1777_F1", "LA1777"))

################
# target genes #
################

###########################
# Layout for all heatmaps #
###########################

annotation_cols = as.data.frame(samples)
col_order = c("F2-151", "F2-411", "F2-445", "PI127826_F1","Elite_01", "PI127826", "F2-28", "F2-73", "F2-127")
row.names(annotation_cols) <- samples$sample
annotation_cols$sample <- NULL

my_colour = list(
  condition = c(elite = "gray99", F1= "gray99", wild = "gray0", F2 = "gray75"),
  pathway= c(MEP = "forestgreen", MVA = "midnightblue")
)

# MEP/MVA precursors

precursor_genes <- read.delim("figure_7/precursor_genes.tsv", header = T, stringsAsFactors = F)
# Remove (putative) Nudix and IPK genes as this is still unclear in tomato
# precursor_genes = precursor_genes %>% filter(!name %in% c("Nudix", "IPK"))
# filter the scaled counts using the target genes
df_filtered <- inner_join(precursor_genes,df_parsed, by = "target_id")
df_filtered <- df_filtered[order(df_filtered$pathway),]
# transform into wide format 
df_filtered_wide <- pivot_wider(df_filtered, id_cols = "target_id", names_from = "sample", values_from = "est_counts") 


# Prepare df for heatmap
mat = as.data.frame(df_filtered_wide[,-1]) 
row.names(mat) = df_filtered_wide$target_id

mat_log2_scaled <- log2(mat + 1)
mat_log2_scaled = mat_log2_scaled[,col_order] 
# heatmap

annotation_rows = as.data.frame(precursor_genes[,-1])
annotation_rows = annotation_rows[order(annotation_rows$pathway),] # this makes the rows separate in MEP/MVA
row.names(annotation_rows) <- precursor_genes$target_id # Connect genes to row annotation


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
         gaps_col = 5 ,filename = "Figure_7/precursors.pdf")
        

        
#####################
# Terpene Synthases #
#####################

TPS <- read.delim("figure_7/terpene_synthases_zhou2020.tsv", header = T, stringsAsFactors = F)
TPS = TPS %>% filter(structure == "functional") # Remove pseudogenes 
        
# filter the scaled counts using the target genes
df_filtered_TPS <- inner_join(TPS,df_parsed, by = "target_id")
df_filtered_TPS = df_filtered_TPS[order(df_filtered_TPS$annotation),] # order from low to high TPS number
# transform into wide format 
df_filtered_TPS_wide <- pivot_wider(df_filtered_TPS, id_cols = "target_id", names_from = "sample", values_from = "est_counts") 

# convert to matrix
mat_TPS = as.data.frame(df_filtered_TPS_wide[,-1]) 
row.names(mat_TPS) = df_filtered_TPS_wide$target_id

mat_TPS_log2_scaled <- log2(mat_TPS + 1)
mat_TPS_log2_scaled <- mat_TPS_log2_scaled[,col_order]

# heatmap
# Order on localisation
annotation_rows = left_join(df_filtered_TPS_wide, TPS, by = "target_id") %>% select(target_id, annotation, localisation) %>% column_to_rownames(var = "target_id")
annotation_rows = annotation_rows[order(annotation_rows$localisation),]
mat_TPS_log2_scaled = mat_TPS_log2_scaled[row.names(annotation_rows),]

my_colour2 = list(
  condition = c(elite = "gray99", F1= "gray99", wild = "gray0", F2 = "gray75"),
  pathway= c(MEP = "forestgreen", MVA = "midnightblue"),
  localisation = c(cytosol = "midnightblue", mitochondria = "firebrick", plastid = "forestgreen")
)



pheatmap(mat_TPS_log2_scaled, 
         cluster_rows = F, 
         cluster_cols = F,
         fontsize = 10,
         cellwidth = 15,
         cellheight = 15,
         gaps_col = 5,
         gaps_row = c(11,13),
         annotation_row = annotation_rows,
         annotation_col = annotation_cols,
         annotation_colors = my_colour2,
         height = 9,
         filename = "Figure_7/TPS_heatmap.pdf")



###########################
# Trans-prenyltranferases #
###########################

TPT <- read.delim("figure_7/trans_prenyltransferases_zhou2020.tsv", header = T, stringsAsFactors = F)

# filter the scaled counts using the target genes
df_filtered_TPT <- inner_join(TPT,df_parsed, by = "target_id")
df_filtered_TPT = df_filtered_TPT[order(df_filtered_TPT$annotation),] # order from low to high TPT number
# transform into wide format 
df_filtered_TPT_wide <- pivot_wider(df_filtered_TPT, id_cols = "target_id", names_from = "sample", values_from = "est_counts") 

# Scaling
# convert to matrix
mat_TPT = as.data.frame(df_filtered_TPT_wide[,-1]) 
row.names(mat_TPT) = df_filtered_TPT_wide$target_id

mat_TPT_log2_scaled <- log2(mat_TPT + 1)
mat_TPT_log2_scaled <- mat_TPT_log2_scaled[,col_order]

# heatmap

TPT_expressed = as.data.frame(left_join(df_filtered_TPT_wide, TPT, by = "target_id"))
row.names(TPT_expressed) <- TPT_expressed$target_id
TPT_expressed = TPT_expressed %>% select("annotation")
annotation_rows = TPT_expressed

annotation_rows = left_join(df_filtered_TPT_wide, TPT, by = "target_id") %>% select(target_id, annotation, localisation) %>% column_to_rownames(var = "target_id")
annotation_rows = annotation_rows[order(annotation_rows$localisation),]
mat_TPT_log2_scaled = mat_TPT_log2_scaled[row.names(annotation_rows),]


pheatmap(mat = mat_TPT_log2_scaled, 
         scale = "none", 
         cluster_rows = F, 
         cluster_cols = F,
         fontsize = 10,
         cellwidth = 15,
         cellheight = 15,
         gaps_col = 5,
         gaps_row = c(2,4),
         annotation_row = annotation_rows,
         annotation_col = annotation_cols,
         annotation_colors = my_colour2,
         height = 6,
         filename = "Figure_7/TPT_heatmap.pdf")


