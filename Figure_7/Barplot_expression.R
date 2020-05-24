library(tidyverse)
library(ggbiplot)
library(ggfortify)

df = read.delim(file = "Figure_7/abundance_tidy.tsv", header = T)

# Shorten target_id names. Then Remove LA1777  related samples from df
df.parsed = df %>% mutate(target_id = substr(target_id, start = 1, stop = 14)) %>%
  filter(sample %in% c("Elite_01", "PI127826_F1", "F2-151", "F2-411", "F2-445",
                       "PI127826", "F2-28", "F2-73", "F2-127"))



# Create df with active/lazy condition per samples
con = data.frame(sample = c("Elite_01", "PI127826_F1", "F2-151", "F2-411", "F2-445",
                            "PI127826", "F2-28", "F2-73", "F2-127"),
                 condition = c("lazy","lazy","lazy","lazy","lazy",
                               "active","active","active","active"))

# Fuse the active/lazy condition with the main df
df.parsed = left_join(df.parsed, con, by = "sample")

# Custom order of the samples
df.parsed$sample = factor(df.parsed$sample, levels = c("Elite_01", "PI127826_F1", "F2-151", "F2-411", "F2-445",
                                                       "PI127826", "F2-28", "F2-73", "F2-127"),
                          ordered = TRUE)

##################################################
# Load Differentially expressed genes from DEse2 #
##################################################

diff <- read.delim(file = "Figure_7/DEseq_analysis/F2_RNAseq_DEseq_resuts_significant_genes.tsv", header = TRUE) %>% 
  mutate(target_id = substr(target_id, start = 1, stop = 14)) %>% 
  arrange(diff$log2FoldChange)

diff.top10 <- diff[1:10,1]

# Filter by solycnumber and plot
p.top10 = df.parsed %>% filter(target_id %in% diff.top10) %>%
  ggplot(aes(x = reorder(sample, condition), y = est_counts, fill = condition))+
  geom_bar(stat = "identity")+
  facet_wrap(~target_id, scale = "free")

ggsave(file = "Figure_7/plots/top10_diff_expressed.pdf", plot = p.top10)


#############################
# Selection based on counts #
#############################

df.wide = df.parsed  %>% select(-condition) %>% pivot_wider(names_from = sample, values_from = est_counts)
df.wide = distinct(df.wide, target_id, .keep_all = TRUE)
write.table(df.wide, file = "Figure_7/abundance_wide.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)


#######
# PCA #
#######
df.for.pca <- df.wide %>% column_to_rownames(var = "target_id")

#Select for top 1000 highest expressed genes
df.for.pca$sum <- rowSums(df.for.pca)
df.for.pca = df.for.pca[order(-df.for.pca$sum),]
df.for.pca.top1000 = df.for.pca[1:10000,] %>% select(-sum)

df.for.pca <- t(df.for.pca.top1000)
df.for.pca <- log(df.for.pca +1)

pca <- prcomp(df.for.pca, scale = T, center = T)

conditions = c("lazy", "active", "lazy", "active", "lazy", "lazy", "active", "active", "lazy")
# ggbiplot(pca, groups = conditions, obs.scale = 1, var.scale = 1)
autoplot(pca, label = TRUE)
