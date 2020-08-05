library(ape)
library(tidyverse)
library(chromoMap)
library(RIdeogram)

# load gff data with gene locations
# Download link:ftp://ftp.solgenomics.net/tomato_genome/Heinz1706/annotation/ITAG4.0_release/

gene.models <- read.gff(file = "/Volumes/Samsung_T5/F2_RNA_seq/ITAG4.0_gene_models.gff")

# Shape datafile to cotain usefull info
gene.models <- gene.models %>% filter(type == "gene") %>%
  select(seqid, start, end, attributes) %>%
  mutate(attributes = substr(attributes, start = 9, stop = 22))


# load DE genes
DE <- read.delim(file = "Figure_8/F2_RNAseq_DEseq_resuts_significant_genes.tsv",
                 header = TRUE)

names(DE)[1] <- "attributes"

# Filter gene.models dataset for DE genes

annotation_file<- left_join(DE, gene.models, by = "attributes")
annotation_file <- annotation_file %>% select(attributes, seqid, start, end, log2FoldChange)
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
png(file = "Figure_8/Chromomap/Chr_map_diff_expression.png", width = 5, height = 10)
 chromoMap("Figure_8/Chromomap/chroms.txt",
          "Figure_8/Chromomap/annotation_file.txt",
          data_based_color_map = TRUE,
          data_type = "numeric",
          canvas_height = 1400,
          canvas_width = 1200,
          chr_length = 5,
          chr_width = 15,
          ch_gap = 20,
          labels = F,
          legend= T,
          lg_x = 20,
          lg_y = 20,
          v_align = T,
          chr_text = F,
          chr_color = c("darkslategray"))
  dev.off()

#############
# RIdeogram #
#############

######################
# Create Chromosomes #
######################

chromosomes <- read.delim(file = "Figure_8/Chromomap/S_lycopersicum_chromosomes.4.00.txt",
                     header = FALSE, sep = "\t")
chromosomes$V5 = chromosomes$V4 + 100 # creares an "alternative centromer"
colnames(chromosomes) <- c("Chr", "Start", "End", "CE_start", "CE_end")
chromosomes$Chr = as.character(chromosomes$Chr)

####################################################
# Get gene coordinates and their expression values #
####################################################

gene.models <- read.gff(file = "/Volumes/Samsung_T5/F2_RNA_seq/ITAG4.0_gene_models.gff")

# Shape datafile to cotain usefull info
gene.models <- gene.models %>% filter(type == "gene") %>%
  select(seqid, start, end, attributes) %>%
  mutate(attributes = substr(attributes, start = 9, stop = 22))


# load all results from DE analysis
gene.expression <- read.delim(file = "Figure_8/F2_RNAseq_DEseq_resuts_significant_genes.tsv",
                 header = TRUE)

names(gene.expression)[1] <- "attributes"

# Filter gene.models dataset for DE genes
genes.for.ideogram <- left_join(gene.expression, gene.models, by = "attributes")
genes.for.ideogram  <- genes.for.ideogram  %>% select(seqid, start, end, log2FoldChange)
colnames(genes.for.ideogram) <- c("Chr", "Start", "End", "Value")

ideogram(karyotype = chromosomes, overlaid = genes.for.ideogram, label_type = "marker")
convertSVG("chromosome.svg", device = "png")



  
