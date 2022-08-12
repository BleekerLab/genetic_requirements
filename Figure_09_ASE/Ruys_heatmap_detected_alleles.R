library(tidyverse)
library(pheatmap)
# Load data
df <- read.delim(file = "Figure_9/DE_genes_detected_alleles.txt",
                 header = TRUE,
                 check.names = FALSE)
# Make df tidy
df.long <- pivot_longer(df, 
                   cols = -c(target_id, chromosome),
                   names_to = "genotype",
                   values_to = "allele")

# Remove transcripts of which allele version could no be determined ("NA") 
df.long.chromosome <- na.omit(df.long) %>% 
  filter(chromosome == "Chr08") 



df.wide <- df.long.chromosome %>% select(-chromosome) %>%
  pivot_wider(names_from = genotype, values_from = allele)

df.wide = df.wide %>% column_to_rownames(var = "target_id")
pheatmap(df.wide,
         cluster_rows = FALSE,
         cluster_cols = FALSE)

