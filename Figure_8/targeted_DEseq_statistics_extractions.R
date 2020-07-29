library(tidyverse)

# Load statistical results from DESeq2 #
DE <- read.table(file = "Figure_8/F2_RNAseq_DEseq_resuts.tsv", header = T) %>%
  select(target_id, log2FoldChange, padj) #select only their genes and


##############
# Precursors #
##############

precursors <- read.delim(file = "Figure_7/precursor_genes.tsv", header = TRUE) %>% 
  mutate(target_id = substr(target_id, start = 1, stop = 14)) 

# fish out expression levels from datafram

precursors.stats<- left_join(precursors, DE, by = "target_id")

write.table(precursors.stats, file = "Figure_8/targeted_DE_stats/precursor_stats.tsv", sep = "\t", row.names = FALSE)

######################
# Prenyltransferases #
######################
prenyl <- read.delim(file = "Figure_7/trans_prenyltransferases_zhou2020.tsv", header = TRUE) %>% 
  mutate(target_id = substr(target_id, start = 1, stop = 14)) 

# fish out expression levels from datafram

prenyl.stats <- left_join(prenyl, DE, by = "target_id")

write.table(prenyl.stats, file = "Figure_8/targeted_DE_stats/prenyltransferases_stats.tsv", sep = "\t", row.names = FALSE)

#####################
# Terpene synthases #
#####################

TPS <- read.delim(file = "Figure_7/terpene_synthases_zhou2020.tsv", header = TRUE) %>% 
  mutate(target_id = substr(target_id, start = 1, stop = 14)) 

# fish out expression levels from datafram

TPS.stats <- left_join(TPS, DE, by = "target_id")

write.table(TPS.stats, file = "Figure_8/targeted_DE_stats/terpene_synthases_stats.tsv", sep = "\t", row.names = FALSE)

##############
# Glycolysis #
##############

GLY<- read.delim("figure_7/glycolysis.tsv", header = T, stringsAsFactors = F)

# filter the scaled counts using the target genes
GLY.stats<- inner_join(GLY,DE, by = "target_id")

write.table(GLY.stats, file = "Figure_8/targeted_DE_stats/glycolysis_stats.tsv", sep = "\t", row.names = FALSE)


#####################
# Citrate shuttle  #
####################
CIT<- read.delim("figure_7/citrate_shuttle.tsv", header = T, stringsAsFactors = F)

# filter the scaled counts using the target genes
CIT.stats <- inner_join(CIT,DE, by = "target_id")

write.table(CIT.stats, file = "Figure_8/targeted_DE_stats/citrate_shuttle_stats.tsv", sep = "\t", row.names = FALSE)

