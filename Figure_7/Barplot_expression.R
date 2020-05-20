library(tidyverse)
library(ggbiplot)

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

# load differentially expressed genes

diff <- read.delim(file = "Figure_7/DEseq_analysis/F2_RNAseq_DEseq_resuts_significant_genes.tsv", header = TRUE) %>% 
  mutate(target_id = substr(target_id, start = 1, stop = 14)) %>% 
  arrange(log2FoldChange)

diff.top10 <- diff[1:10,1]

# Filter by solycnumber and plot
df.parsed %>% filter(target_id %in% diff.top10) %>%
  ggplot(aes(x = reorder(sample, condition), y = est_counts, fill = condition))+
  geom_bar(stat = "identity")+
  facet_wrap(~target_id, scale = "free")


#############################
# Selection based on counts #
#############################

df.wide counts = df.parsed  %>% pivot_wider(names_from = sample, values_from = est_counts)
write.table(df.wide, file = "Figure_7/abundance_wide.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)

df.filtered = df.wide %>% filter(
                                   
                                   "F2-28" > "F2-151" &
                                   "F2-28" > "F2-411") &
                                   "F2-28" > "F2-455" 
)

#######
# PCA #
#######
df.for.pca <- df.wide %>% column_to_rownames(var = "target_id")
df.for.pca <- t(df.for.pca)
df.for.pca <- log(df.for.pca +1)
pca <- prcomp(df.for.pca)

ggbiplot(pca, groups = conditions)
