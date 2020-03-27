library("checkpoint")

checkpoint("2020-01-01")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))

source("Figure_8/R/normalise_counts.R")

###############
# Import counts
###############

# Counts based on the Heinz1706 genome 
heinz <- read.delim("Figure_8/counts_on_heinz.tsv", stringsAsFactors = F, check.names = F)

# Counts based on the PI127826 genome
habro <- read.delim("Figure_8/counts_on_habrochaites.tsv", stringsAsFactors = F, check.names = F)

# remove columns not useful (Chr, Start, End, etc.)
heinz <- heinz[,-c(2,3,4,5,6)] 
habro <- habro[,-c(2,3,4,5,6)] 

#######################
# Normalise with DESeq2
#######################

row.names(heinz) <- heinz$target
row.names(habro) <- habro$target

cts_heinz <- heinz[,-1]
cts_habro <- habro[,-1]

col_data = data.frame(
  sample = colnames(heinz)[2:ncol(heinz)], 
  condition = c(rep("A",4), rep("B","5"))) # arbitrary conditions (not used)

dds_heinz <- DESeqDataSetFromMatrix(countData = cts_heinz, colData = col_data, design = ~ condition)
dds_habro <- DESeqDataSetFromMatrix(countData = cts_habro, colData = col_data, design = ~ condition)

dds_heinz <- estimateSizeFactors(dds_heinz)
dds_habro <- estimateSizeFactors(dds_habro)

norm_heinz <- DESeq2::counts(dds_heinz, normalized = TRUE)
norm_habro <- DESeq2::counts(dds_habro, normalized = TRUE)

###################
# Keep common genes
###################

# keep only common target genes
norm_heinz <- norm_heinz[row.names(norm_heinz) %in% habro$target,]

# reorder genes in the same order
norm_habro <- norm_habro[order(row.names(norm_habro)),]

# sanity check
if (! all_equal(row.names(norm_habro), row.names(norm_heinz)) ) {
  stop("problems with matching gene names")

###################################
# Conversion to matrix and addition
###################################

total_counts <- norm_heinz + norm_habro

# % of counts coming from ELite (Heinz1706 genome)
perc_heinz <- norm_heinz / total_counts * 100 

# % of counts coming from habro (PI127826 genome)
perc_habro <- norm_habro / total_counts * 100 


# add genome ref
perc_habro <- as.data.frame(perc_habro)
perc_heinz <- as.data.frame(perc_heinz)
perc_habro$genome = "habrochaites"
perc_heinz$genome = "lycopersicum"

# add genes names
perc_habro$target <- targets
perc_heinz$target <- targets

# add genome ref
perc_habro$genome = "habrochaites"
perc_heinz$genome = "lycopersicum"

#############################################################
# Creation of a combined dataframe with a "genome_ref" column
# Goal = plot comparative barplots
#############################################################

final_df = bind_rows(perc_habro, perc_heinz) %>% 
  pivot_longer(., 
               cols = - c(target, genome),
               names_to = "genotype", 
               values_to = "percentage_of_counts") 

p_barplot <- final_df %>% 
  ggplot(., aes(x = genotype, y = percentage_of_counts, fill = genome)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ target) + 
    coord_flip() +
    labs(x = "Genotype", y = "Percentage of read counts from parental genome") +
    theme(axis.text.x = element_text(angle = 0))
  
ggsave(filename = "Figure_8/figure8.png", plot = p_barplot, width = 10, height = 7)
