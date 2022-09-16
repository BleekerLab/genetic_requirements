library(tidyverse)
library(DESeq2)

# Open non-normalised counts 
raw.counts <- read.delim("data/RNA_seq_2022/raw_counts.parsed.tsv")

raw.counts$Geneid <- str_sub(raw.counts$Geneid, start = 6, 19)

raw.counts <- raw.counts %>% column_to_rownames("Geneid")

raw.counts.precursors <- raw.counts[rownames(raw.counts) %in% precursors$gene_id, ]

# Open sample information
samples <- read.delim("data/RNA_seq_2022/samples_to_conditions.csv", sep = ",") %>% 
  select(sample_id, genotype) %>%
  column_to_rownames("sample_id") 


########################
# normalise the counts #
########################

# create a DESeq2 dataset 
dds.data <- DESeqDataSetFromMatrix(countData = raw.counts,
                              colData = samples,
                              design = ~genotype)

##################
# DESeq analysis #
##################

dds <- DESeq(dds.data)

res <- results(dds, contrast = c("genotype", "P28", "Elite"))

res <- lfcShrink(dds,
          contrast = c("genotype", "P28", "Elite"), res=res, type = 'ashr')

# extract a table with the DESeq results
# To be used showing the MEP/MVA pathway genes
fold.changes <- bind_cols(gene_id = res@rownames, 
                          bind_cols(res@listData))

m <- fold.changes.precursors %>% select(gene_id, log2FoldChange) %>% 
  column_to_rownames("gene_id")

fold.changes.precursors <- fold.changes %>% 
  filter(gene_id %in% precursors$gene_id) %>% 
  left_join(.,precursors, by = "gene_id")

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')

############################
# Obtain normalised counts #
############################

dds.size.factors <- estimateSizeFactors(dds.data)

df.normalised <- counts(dds.size.factors, normalized = TRUE)

#############################################
# plot normalised counts of precursor genes #
#############################################


df <- as.data.frame(df.normalised) %>% 
  rownames_to_column(var = "gene_id") 
  

# Open genes of interest and description
precursors <- read.delim("figures/Figure_6_RNA_seq/precursor_genes.tsv") %>% 
  dplyr::rename(gene_id = target_id)

samples <- samples %>% rownames_to_column(var = "sample")

df.precursors <- df %>% 
  filter(gene_id %in% precursors$gene_id) %>% 
  pivot_longer(-gene_id, names_to = "sample", values_to = "value") %>% 
  left_join(., samples, by = "sample") %>% 
  left_join(., precursors, by = "gene_id") %>% 
  mutate(genotype = if_else(genotype == "Elite", "Cultivar", genotype))
  



df.precursors %>% 
  filter(pathway == "MEP") %>% 
  ggplot(aes(x = genotype, y = value)) + 
  geom_boxplot(aes(fill = genotype))+
  geom_point()+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~gene_id)+
  theme_bw()
  

df.precursors %>% 
  filter(pathway == "MVA") %>% 
  ggplot(aes(x = genotype, y = value)) + 
  geom_boxplot(aes(fill = genotype))+
  geom_point()+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~gene_id)

# Selected precursors
df.precursors %>% 
  filter(gene_id %in% c("Solyc01g109300",
                        "Solyc08g005680",
                        "Solyc11g010850",
                        "Solyc11g069380",
                        "Solyc02g038740",
                        "Solyc03g032020",
                        "Solyc07g045350",
                        "Solyc08g080170",
                        "Solyc11g007020")) %>% 
  ggplot(aes(x = genotype, y = value)) + 
  geom_boxplot(aes(fill = genotype))+
  geom_point()+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(pathway~gene_id, ncol = 5)+
  theme_bw()

