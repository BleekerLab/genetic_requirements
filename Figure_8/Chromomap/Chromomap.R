library(ape)
library(tidyverse)
library(chromoMap)

# load gff data with gene locations
# Download link:ftp://ftp.solgenomics.net/tomato_genome/Heinz1706/annotation/ITAG4.0_release/

gene.models <- read.gff(file = "/Volumes/Samsung_T5/F2_RNA_seq/ITAG4.0_gene_models.gff")

# Shape datafile to cotain usefull info
gene.models <- gene.models %>% filter(type == "gene") %>%
  select(seqid, start, end, attributes) %>%
  mutate(attributes = substr(attributes, start = 9, stop = 22))


# load DE genes
DE <- read.delim(file = "Figure_8/DESeq_analysis_without_parents/F2_RNAseq_DEseq_resuts_significant_genes.tsv",
                 header = TRUE)

names(DE)[1] <- "attributes"

# Filter gene.models dataset for DE genes

annotation_file<- left_join(DE, gene.models, by = "attributes")
annotation_file <- annotation_file %>% select(attributes, seqid, start, end)
# Save aas annotation file

write.table(annotation_file, "Figure_8/Chromomap/annotation_file.txt",
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

chroms <- read.delim(file = "Figure_8/Chromomap/S_lycopersicum_chromosomes.4.00.txt",
                     header = FALSE, sep = "\t")

chroms = as.data.frame(chroms)

write.table(chroms, file = "Figure_8/Chromomap/chroms.txt",
            sep = "\t", col.names = FALSE, row.names = FALSE,
            quote = FALSE)


#############
# ChromoMap #
#############

ch <- chromoMap("Figure_8/Chromomap/chroms.txt",
          "Figure_8/Chromomap/annotation_file.txt")


  
