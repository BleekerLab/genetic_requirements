library(ape)
library(tidyverse)
library(chromoMap)

# load gff data with gene locations
gene.models <- read.gff(file = "/Volumes/Samsung_T5/F2_RNA_seq/ITAG4.0_gene_models.gff")

# Shape datafile to cotain usefull info
gene.models <- gene.models %>% filter(type == "gene") %>%
  select(seqid, start, end, attributes) %>%
  mutate(attributes = substr(attributes, start = 9, stop = 22))

# load DE genes

# Filter gene.models dataset for DE genes

# Save aas annotation file




#############
# ChromoMap #
#############

  
