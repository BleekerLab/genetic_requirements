suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(patchwork))

raw <- read.delim("Supplemental_data_RNA-seq/raw_counts.tsv", header = T) 
scaled <- read.delim("Supplemental_data_RNA-seq/scaled_counts.tsv", header = T) 


tidy_raw <- raw %>% pivot_longer(cols = - Geneid, names_to = "genotype", values_to = "raw_counts") %>% 
  rename(gene = Geneid)

tidy_scaled <- scaled %>% pivot_longer(cols = - gene, names_to = "genotype", values_to = "scaled_counts") 
  
df <- inner_join(tidy_raw, tidy_scaled, by = c("gene","genotype"))

ggplot(data = df, aes(x = raw_counts + 1, y = scaled_counts + 1, color = genotype)) +
  geom_point() +
  facet_wrap(~ genotype) + scale_x_log10() + scale_y_log10()


ggsave(filename = "Supplemental_data_RNA-seq/scatterplots_raw_vs_scaled.png", width = 8, height = 8)

p_boxplot_raw <- ggplot(data = df, aes(x = genotype, y = raw_counts + 1)) +
  geom_boxplot() +
  scale_y_log10() + 
  ggtitle("Raw counts") +
  theme(axis.text.x = element_text(angle = 90))

p_boxplot_scaled <- ggplot(data = df, aes(x = genotype, y = scaled_counts + 1)) +
  geom_boxplot() +
  scale_y_log10() + 
  ggtitle("Scaled counts") +
  theme(axis.text.x = element_text(angle = 90))

p_boxplot_raw + p_boxplot_scaled

ggsave(filename = "Supplemental_data_RNA-seq/raw_vs_scaled_boxplot_counts.png", width = 10, height = 7)


