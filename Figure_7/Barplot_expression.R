library(tidyverse)
library(ggbiplot)
library(ggfortify)

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


# Save in wide format (for PCA)
# df.wide = df.parsed  %>% select(-condition) %>% pivot_wider(names_from = sample, values_from = est_counts)
# df.wide = distinct(df.wide, target_id, .keep_all = TRUE)
# write.table(df.wide, file = "Figure_7/abundance_wide.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)


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

#####################################################
# Filter dataframe by diff expressed genes and plot #
#####################################################
  p.top10 =
 df.parsed %>% filter(target_id %in% diff.top10) %>%
  ggplot(aes(x = reorder(sample, condition), y = est_counts, fill = condition))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("lazy" = "steelblue", "active" = "darkred"))+
  facet_wrap(~target_id, scale = "free")+
  labs(x = "Sample" , y = "Gene expression (counts)")+
  my_theme

ggsave(file = "Figure_7/plots/top10_diff_expressed.pdf", plot = p.top10)

###########################################
# Load Target genes list for barplotting  #
###########################################


##############
# Precursors #
##############

precursors <- read.delim(file = "Figure_7/precursor_genes.tsv", header = TRUE) %>% 
  mutate(target_id = substr(target_id, start = 1, stop = 14)) 

# fish out expression levels from datafram

precursors.expression <- left_join(precursors, df.parsed, by = "target_id")

# Make sure order of samples is OK
precursors.expression$sample = factor(precursors.expression$sample, levels = c("Elite_01", "PI127826_F1", "F2-151", "F2-411", "F2-445",
                                                       "PI127826", "F2-28", "F2-73", "F2-127"),
                          ordered = TRUE)

# plot MEP genes

p.mep =
precursors.expression %>% filter(pathway == "MEP")%>%
  ggplot(aes(x = sample, y = est_counts, fill = condition))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("lazy" = "steelblue", "active" = "darkred"))+
  facet_wrap(~name, scale = "free")+
  labs(x = "Sample" , y = "Gene expression (counts)")+
  my_theme

ggsave(file = "Figure_7/plots/MEP_genes_barplot.pdf", plot = p.mep)

# plot MVA genes  

p.mva =
  precursors.expression %>% filter(pathway == "MVA")%>%
  ggplot(aes(x = sample, y = est_counts, fill = condition))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("lazy" = "steelblue", "active" = "darkred"))+
  facet_wrap(~name, scale = "free")+
  labs(x = "Sample" , y = "Gene expression (counts)")+
  my_theme

ggsave(file = "Figure_7/plots/MVA_genes_barplot.pdf", plot = p.mva)

######################
# Prenyltransferases #
######################
prenyl <- read.delim(file = "Figure_7/trans_prenyltransferases_zhou2020.tsv", header = TRUE) %>% 
  mutate(target_id = substr(target_id, start = 1, stop = 14)) 

# fish out expression levels from datafram

prenyl.expression <- left_join(prenyl, df.parsed, by = "target_id")

prenyl.expression$sample = factor(prenyl.expression$sample, levels = c("Elite_01", "PI127826_F1", "F2-151", "F2-411", "F2-445",
                                                                               "PI127826", "F2-28", "F2-73", "F2-127"),
                                      ordered = TRUE)

p.prenyl =
  
  prenyl.expression %>% 
  ggplot(aes(x = sample, y = est_counts, fill = condition))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("lazy" = "steelblue", "active" = "darkred"))+
  facet_wrap(~annotation, scale = "free")+
  labs(x = "Sample" , y = "Gene expression (counts)")+
  my_theme

ggsave(file = "Figure_7/plots/prenyltransferases_genes_barplot.pdf", plot = p.prenyl)



#####################
# Terpene synthases #
#####################

TPS <- read.delim(file = "Figure_7/terpene_synthases_zhou2020.tsv", header = TRUE) %>% 
  mutate(target_id = substr(target_id, start = 1, stop = 14)) 

# fish out expression levels from datafram

TPS.expression <- left_join(TPS, df.parsed, by = "target_id")

TPS.expression$sample = factor(TPS.expression$sample, levels = c("Elite_01", "PI127826_F1", "F2-151", "F2-411", "F2-445",
                                                                       "PI127826", "F2-28", "F2-73", "F2-127"),
                                  ordered = TRUE)

p.tps =
  
  TPS.expression %>% filter(structure == "functional")%>% filter(est_counts > 0) %>%
  ggplot(aes(x = sample, y = est_counts, fill = condition))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("lazy" = "steelblue", "active" = "darkred"))+
  facet_wrap(~annotation, scale = "free")+
  labs(x = "Sample" , y = "Gene expression (counts)")+
  my_theme

ggsave(file = "Figure_7/plots/TPS_genes_free_scale_barplot.pdf", width = 12, height = 12, plot = p.tps)


