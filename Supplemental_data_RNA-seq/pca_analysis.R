################
# PCA score plot
################

suppressPackageStartupMessages(library("tidyverse"))
source("scripts/mypca.R")

df <- read.delim("Supplemental_data_RNA-seq/scaled_counts.tsv", stringsAsFactors = F) %>% 
  mutate(gene = gsub(pattern = "mRNA:", replacement = "", x = gene)) %>%
  column_to_rownames("gene") %>% 
  t(.)

df[1:5,1:5]

pca_res <- mypca(df, center = TRUE, scale = TRUE)


############
# Scree plot
############
dfev <- data.frame(PC = seq(1,7), exp_var  = pca_res$explained_var)
head(dfev)

# make the plot
scree_plot <- ggplot(dfev, aes(x = PC, y = exp_var)) +
  ylab('explained variance (%)') + 
  ggtitle('explained variance per component') + 
  geom_bar(stat = "identity")

# display it
scree_plot

############
# Score plot
############
scores <- pca_res$scores
scores$genotype <- rownames(df)

p <- ggplot(scores, aes(x = PC1, y = PC2, col = genotype, label = genotype)) + 
  geom_point(size = 3) + 
  xlab(paste0('PC1(',pca_res$explained_var[1,1],'%)')) + 
  ylab(paste0('PC2(',pca_res$explained_var[2,1],'%)')) + 
  ggtitle('PCA score plot') +
  ggrepel::geom_text_repel() +
  scale_color_brewer(type = "qual", palette = "Dark2")
p

ggsave("Supplemental_data_RNA-seq/pca_score_plot.pdf", width = 10, height = 8)

