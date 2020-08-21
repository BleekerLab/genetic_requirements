library(tidyverse)
library(pheatmap)
# Load data
df <- read.delim(file = "Figure_9/DE_genes_detected_alleles.txt",
                 header = TRUE,
                 check.names = FALSE)
# Make df tidy
df <- pivot_longer(df, 
                   cols = -target_id,
                   names_to = "genotype",
                   values_to = "allele")

# Remove transcripts of which allele version could no be determined ("NA") 
df <- na.omit(df)

df.wide <- df %>% pivot_wider(names_from = genotype, values_from = allele)
df.wide = df.wide %>% column_to_rownames(var = "target_id")

pheatmap(df.wide,
         cluster_rows = FALSE,
         cluster_cols = FALSE)

