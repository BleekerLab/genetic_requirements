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
df <- read.delim(file = "Figure_8/normalised_counts_by_DEseq2.tsv", header = T)


df = df %>% filter(!genotype %in% c("Elite_02", "F1_hab", "LA1777_F1", "LA1777", "F2_28"))

## Step two: importing the sample information and re-ordering it
samples <- read_tsv("Figure_7/samples.tsv", col_names = TRUE )[c("genotype", "condition")]
samples$condition <- with(samples, factor(condition, levels = c("elite","F1","wild","F2"))) 
samples = samples %>% filter(!sample %in% c("Elite_02", "F1-hab", "LA1777_F1", "LA1777","F2-28"))

################
# target genes #
################

###########################
# Layout for all heatmaps #
###########################

annotation_cols = as.data.frame(samples)
col_order = c("F2_151", "F2_411", "F2_445", "PI127826_F1","Elite_01", "PI127826", "F2_73", "F2_127")
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
df_filtered <- inner_join(precursor_genes,df, by = "target_id")
df_filtered <- df_filtered[order(df_filtered$pathway),]
# transform into wide format 
df_filtered_wide <- pivot_wider(df_filtered, id_cols = "target_id", names_from = "genotype", values_from = "count") 


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
         #annotation_col = annotation_cols, 
         annotation_row = annotation_rows,
         nnotation_colors = my_colour,
         gaps_row = 12,
         gaps_col = 5 ,filename = "Figure_7/heatmaps/MEP_MVA_heatmap_normalised_counts.pdf")
        

        
#####################
# Terpene Synthases #
#####################

TPS <- read.delim("figure_7/terpene_synthases_zhou2020.tsv", header = T, stringsAsFactors = F)
TPS = TPS %>% filter(structure == "functional") # Remove pseudogenes 
        
# filter the scaled counts using the target genes
df_filtered_TPS <- inner_join(TPS,df_parsed, by = "target_id")
df_filtered_TPS = df_filtered_TPS[order(df_filtered_TPS$annotation),] # order from low to high TPS number
# transform into wide format 
df_filtered_TPS_wide <- pivot_wider(df_filtered_TPS, names_from = "sample", values_from = "est_counts") 
  

# convert to matrix
mat_TPS = as.data.frame(df_filtered_TPS_wide[,-(1:4)]) 
row.names(mat_TPS) = df_filtered_TPS_wide$target_id

mat_TPS_log2_scaled <- log2(mat_TPS + 1)
mat_TPS_log2_scaled <- mat_TPS_log2_scaled[,col_order]

# heatmap
# Order on localisation
annotation_rows <- df_filtered_TPS_wide %>% select(target_id, annotation, localisation) %>% column_to_rownames(var = "target_id")
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
         filename = "Figure_7/heatmaps/TPS_heatmap.pdf")



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
         filename = "Figure_7/heatmaps/TPT_heatmap.pdf")

#####################
# Citrate shuttle  #
####################
CIT<- read.delim("figure_7/citrate_shuttle.tsv", header = T, stringsAsFactors = F)

# filter the scaled counts using the target genes
df_filtered_CIT <- inner_join(CIT,df_parsed, by = "target_id")
# transform into wide format 
df_filtered_CIT_wide <- pivot_wider(df_filtered_CIT, id_cols = "target_id", names_from = "sample", values_from = "est_counts") 

# Scaling
# convert to matrix
mat_CIT <- as.data.frame(df_filtered_CIT_wide[,-1]) 
row.names(mat_CIT) <- df_filtered_CIT_wide$target_id

mat_CIT_log2_scaled <- log2(mat_CIT + 1)
mat_CIT_log2_scaled <- mat_CIT_log2_scaled[,col_order]

# heatmap
# Order on localisation
annotation_rows = left_join(df_filtered_CIT_wide, CIT, by = "target_id") %>% select(target_id, annotation, localisation) %>% column_to_rownames(var = "target_id")
annotation_rows = annotation_rows[order(annotation_rows$localisation),]
#mat_TPS_log2_scaled = mat_TPS_log2_scaled[row.names(annotation_rows),]

my_colour2 = list(
  condition = c(elite = "gray99", F1= "gray99", wild = "gray0", F2 = "gray75"),
  pathway= c(MEP = "forestgreen", MVA = "midnightblue"),
  localisation = c(cytosol = "midnightblue", mitochondria = "firebrick", plastid = "forestgreen")
)

pheatmap(mat = mat_CIT_log2_scaled, 
         scale = "none", 
         cluster_rows = F, 
         cluster_cols = F,
         fontsize = 10,
         cellwidth = 15,
         cellheight = 15,
         gaps_col = 5,
         gaps_row =11,
         annotation_row = annotation_rows,
         annotation_col = annotation_cols,
         annotation_colors = my_colour2,
         height = 8,
         filename = "Figure_7/heatmaps/Cit_shuttle_heatmap.pdf")

##############
# Glycolysis #
##############

GLY<- read.delim("figure_7/glycolysis.tsv", header = T, stringsAsFactors = F)

# filter the scaled counts using the target genes
df_filtered_GLY <- inner_join(GLY,df_parsed, by = "target_id")
# transform into wide format 
df_filtered_GLY_wide <- pivot_wider(df_filtered_GLY, id_cols = "target_id", names_from = "sample", values_from = "est_counts") 

# Scaling
# convert to matrix
mat_GLY <- as.data.frame(df_filtered_GLY_wide[,-1]) 
row.names(mat_GLY) <- df_filtered_GLY_wide$target_id

mat_GLY_log2_scaled <- log2(mat_GLY + 1)
mat_GLY_log2_scaled <- mat_GLY_log2_scaled[,col_order]

# heatmap
# Order on localisation
annotation_rows = left_join(df_filtered_GLY_wide, GLY, by = "target_id") %>% select(target_id, annotation, localisation) %>% column_to_rownames(var = "target_id")
annotation_rows = annotation_rows[order(annotation_rows$localisation),]
#mat_TPS_log2_scaled = mat_TPS_log2_scaled[row.names(annotation_rows),]

my_colour2 = list(
  condition = c(elite = "gray99", F1= "gray99", wild = "gray0", F2 = "gray75"),
  pathway= c(MEP = "forestgreen", MVA = "midnightblue"),
  localisation = c(cytosol = "midnightblue", mitochondria = "firebrick", plastid = "forestgreen")
)

pheatmap(mat = mat_GLY_log2_scaled, 
         scale = "none", 
         cluster_rows = F, 
         cluster_cols = F,
         fontsize = 10,
         cellwidth = 15,
         cellheight = 15,
         gaps_col = 5,
         gaps_row =17,
         annotation_row = annotation_rows,
         annotation_col = annotation_cols,
         annotation_colors = my_colour2,
         height = 11,
         filename = "Figure_7/heatmaps/Glycolysis_heatmap.pdf")




