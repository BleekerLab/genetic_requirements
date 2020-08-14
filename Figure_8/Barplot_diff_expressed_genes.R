library(tidyverse)

######################
# Theme for plotting #
######################

my_theme = theme_bw()+
  theme(text = element_text(),
        axis.text.x = element_text(size = 8, colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8, colour = "black"),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.text.x = element_text(size=8, colour = "black"))

#######################
# Load and parse data #
#######################

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

diff <- read.delim(file = "Figure_8/F2_RNAseq_DEseq_resuts_significant_genes.tsv", header = TRUE) %>% 
  arrange(diff$padj)

diff.top10 <- diff[1:10,1]

#####################################################
# Filter dataframe by diff expressed genes and plot #
#####################################################
p.top10 =
  df.parsed %>% filter(target_id %in% diff.top10) %>%
  ggplot(aes(x = reorder(sample, condition), y = est_counts, fill = condition))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("lazy" = "grey", "active" = "black"))+
  facet_wrap(~target_id, scale = "free")+
  labs(x = "Sample" , y = "Gene expression (counts)")+
  my_theme

ggsave(file = "Figure_7/plots/top10_diff_expressed.pdf", plot = p.top10)

df.parsed %>% filter(target_id %in% diff.top10) %>%
  ggplot(aes(x = condition, y = est_counts, fill = condition))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~target_id, scale = "free")+
  labs(x = "Group" , y = "Gene expression (counts)")+
  my_theme
