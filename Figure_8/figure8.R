library("checkpoint")

checkpoint("2020-01-01")

library(tidyverse)

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

###################################
# Keep common genes
###################################

# keep only common target genes
heinz <- heinz[heinz$target %in% habro$target,]

# reorder genes in the same order
habro <- habro %>% arrange(target) 

# sanity check
if (all_equal(habro$target, heinz$target) == FALSE){
  stop()
}

###################################
# Conversion to matrix and addition
###################################
# keep genes names
targets <-heinz$target

heinz <- heinz %>% select(- target) 
habro <- habro %>% select(- target) 

total_counts <- heinz + habro

# % of counts coming from ELite (Heinz1706 genome)
perc_heinz <- heinz / total_counts * 100 

# % of counts coming from habro (PI127826 genome)
perc_habro <- habro / total_counts * 100 

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
               values_to = "counts") 

  #ggplot(., aes(x = counts)) + 
  #geom_density() 
