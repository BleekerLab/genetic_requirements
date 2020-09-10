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

# Load normalised_counts dataset from DEseq2 analysis
df = read.delim(file = "Figure_8/normalised_counts_by_DEseq2.tsv", header = T)


# Create df with active/lazy condition per samples
con = data.frame(genotype = c("Elite_01", "PI127826_F1", "F2_151", "F2_411", "F2_445",
                            "PI127826", "F2_73", "F2_127"),
                 condition = c("lazy","lazy","lazy","lazy","lazy",
                               "active","active","active"))

# Fuse the active/lazy condition with the main df
df = left_join(df, con, by = "genotype")

# Custom order of the samples
df$genotype = factor(df$genotype, levels = c("Elite_01", "PI127826_F1", "F2_151", "F2_411", "F2_445",
                                             "PI127826", "F2_73", "F2_127"),
                          ordered = TRUE)

##################################################
# Load Differentially expressed genes from DEse2 #
##################################################

diff <- read.delim(file = "Figure_8/F2_RNAseq_DEseq_resuts_significant_genes.tsv", header = TRUE) %>% 
  dplyr::arrange(pvalue)

diff.top10 <- diff[1:10,1] # top 10 diff expressed genes
diff.top10 <- c(diff.top10, "Solyc08g079840")

#####################################################
# Filter dataframe by diff expressed genes and plot #
#####################################################
  p.top10 =
 df %>% filter(target_id %in% diff.top10) %>%
  ggplot(aes(x = reorder(genotype, condition), y = count, fill = condition))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("lazy" = "steelblue", "active" = "darkred"))+
  facet_wrap(~target_id, scale = "free")+
  labs(x = "Sample" , y = "Gene expression (counts)")+
  my_theme

ggsave(file = "Figure_7/plots/top10_diff_expressed.pdf", plot = p.top10)

###########
# Heatmap #
###########

#create the df for pheatmap

diff.top25 <- diff[1:25,1] # top 25 diff expressed genes

# Create df for include annotation of the genes in hetmap
annotation.top.25 <- diff[1:25,c(1,8)] 
rownames(annotation.top.25) <- annotation.top.25$target_id # annotations for the heatmap
annotation.top.25 = annotation.top.25 %>% select(-target_id)

df.heatmap <- df %>% filter(target_id %in% diff.top25) %>%
                    select(-condition) %>%
                    pivot_wider(names_from = genotype,
                                values_from = count) %>%
                    column_to_rownames(var = "target_id")

df.heatmap <- df.heatmap[, c("F2_151", "F2_411", "F2_445", "PI127826_F1","Elite_01", "PI127826", "F2_73", "F2_127")]
log.df.heatmap <- log(df.heatmap+1)

###################
# Layout heatmap  #
###################




pheatmap(log.df.heatmap,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize = 8,
         annotation_row = annotation.top.25,
         fontsize_row = 6
         )
###########################################
# Load Target genes list for barplotting  #
###########################################


##############
# Precursors #
##############

precursors <- read.delim(file = "Figure_7/precursor_genes.tsv", header = TRUE) %>% 
  mutate(target_id = substr(target_id, start = 1, stop = 14)) 

# fish out expression levels from datafram

precursors.expression <- left_join(precursors, df, by = "target_id")

# Make sure order of samples is OK
precursors.expression$genotype= factor(precursors.expression$genotype, levels = c("Elite_01", "PI127826_F1", "F2_151", "F2_411", "F2_445",
                                                       "PI127826", "F2_73", "F2_127"),
                          ordered = TRUE)

# plot MEP genes

p.mep =
precursors.expression %>% filter(pathway == "MEP")%>%
  ggplot(aes(x = genotype, y = log(count), fill = condition))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("lazy" = "steelblue", "active" = "darkred"))+
  facet_wrap(~name, scale = "free")+
  labs(x = "Sample" , y = "Gene expression (counts)")+
  my_theme

ggsave(file = "Figure_7/plots/MEP_genes_barplot.pdf", plot = p.mep)

# plot MVA genes  

p.mva =
  precursors.expression %>% filter(pathway == "MVA")%>%
  ggplot(aes(x = genotype, y = count, fill = condition))+
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

prenyl.expression <- left_join(prenyl, df, by = "target_id")

prenyl.expression$genotype = factor(prenyl.expression$genotype, levels = c("Elite_01", "PI127826_F1", "F2_151", "F2_411", "F2_445",
                                                                               "PI127826", "F2_73", "F2_127"),
                                      ordered = TRUE)

p.prenyl =
  
  prenyl.expression %>% 
  ggplot(aes(x = genotype, y = count, fill = condition))+
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

TPS.expression <- left_join(TPS, df, by = "target_id")

TPS.expression$genotype = factor(TPS.expression$genotype, levels = c("Elite_01", "PI127826_F1", "F2_151", "F2_411", "F2_445",
                                                                       "PI127826", "F2_73", "F2_127"),
                                  ordered = TRUE)

p.tps =
  
  TPS.expression %>% filter(structure == "functional")%>% filter(count > 0) %>%
  ggplot(aes(x = genotype, y = count, fill = condition))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("lazy" = "steelblue", "active" = "darkred"))+
  facet_wrap(~annotation)+
  labs(x = "Sample" , y = "Gene expression (counts)")+
  my_theme

ggsave(file = "Figure_7/plots/TPS_genes_free_scale_barplot.pdf", width = 12, height = 12, plot = p.tps)


