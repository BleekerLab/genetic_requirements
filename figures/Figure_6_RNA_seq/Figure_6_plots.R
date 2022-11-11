library(pheatmap)
library(tidyverse)

# Load raw data
df.raw <- read.delim("data/Deseq2_normalised_counts_P28_CV.txt") %>% 
  rename(gene = X) %>%
  pivot_longer(-gene, names_to = "sample", values_to = "normalised_counts")

# Get sample information
sample.description <- read.csv("data/RNA_Seq_2022/samples_to_conditions.csv", sep = ",")

# fuse scaled-counts data with sample info
df <- 
  left_join(df.raw, sample.description, by = c("sample" = "sample_id")) %>% 
  select(gene, sample, genotype, normalised_counts) %>% 
  mutate(genotype = if_else(genotype == "Elite", "Cultivar", genotype)) %>%
  mutate(genotype_sample= paste0(genotype, "_", sample))



# Open genes of interest and description
targets <- read.delim("figures/Figure_6_RNA_seq/precursor_genes.tsv")
targets2 <- targets %>% filter(target_id %in% df$gene) %>% distinct(target_id, .keep_all =T)

# Fuse gene info of the targets and filter them from the dataset
df.targets <- 
  df %>% 
  filter(gene %in% targets$target_id) %>% 
  left_join(., targets, by = c("gene" = "target_id")) %>% 
  mutate(normalised_counts = log(normalised_counts+1, base = 10)) # log scale transformation


#############
# MEP / MVA #
#############

df.mep.mva <- 
  df.targets %>% 
  filter(pathway %in% c("MEP", "MVA")) %>% 
  arrange(factor(gene, levels = as.vector(targets %>% 
                                            filter(pathway %in% c("MEP", "MVA")) %>% 
                                            filter(target_id %in% df.targets$gene) %>% 
                                            .$target_id))) %>%
  select(genotype_sample, gene, normalised_counts) %>% 
  pivot_wider(names_from = genotype_sample, values_from = normalised_counts) %>% 
  column_to_rownames(var = "gene")

# Colors for the plot
my_colour = list(
  pathway= c(MEP = "forestgreen", MVA = "midnightblue")
)

# Annotation of the rows in the heatmap
annotation_rows <- targets %>% filter(pathway %in% c("MEP", "MVA")) %>% column_to_rownames(var = "target_id") %>% select(pathway) # this makes the rows separate in MEP/MVA

pheatmap(mat = df.mep.mva, 
         scale = "none", 
         cluster_rows = F, 
         cluster_cols = F,
         fontsize = 10,
         cellwidth = 15,
         cellheight = 15,
         #annotation_col = annotation_cols, 
         annotation_row = annotation_rows,
         annotation_colors = my_colour,
         gaps_row = 12,
         gaps_col = 5)
        #, filename = "MEP_MVA_heatmap_normalised_counts2.pdf")

###########
# Boxplot #
###########

df.mean <- df %>% 
  group_by(genotype, gene) %>% 
  summarise(mean_counts = mean(normalised_counts))

df %>% 
#  filter(gene %in% c("Solyc11g010850","Solyc11g069380", "Solyc01g109300","Solyc08g005680")) %>%
  ggplot(aes(x = genotype, y = normalised_counts, fill = genotype)) +
  geom_boxplot()+
  geom_point(color = "darkgrey")+
  facet_grid(~gene)+
  scale_fill_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none")
  
ggsave("MEP_genes.svg")

#####################
# Multiple heatmaps #
#####################

for (i in unique(df.targets$pathway)){
  
  df.i<- 
    df.targets %>% 
    filter(pathway %in% i) %>% 
    arrange(factor(gene, levels = as.vector(targets2 %>% 
                                              filter(pathway %in% i) %>% 
                                              filter(target_id %in% df.targets$gene) %>% 
                                              .$target_id))) %>%
    select(genotype_sample, gene, scaled_counts) %>% 
    distinct(genotype_sample, gene, .keep_all = TRUE) %>% 
    pivot_wider(names_from = genotype_sample, values_from = scaled_counts) %>%
    column_to_rownames(var = "gene")
  
  pheatmap(mat = df.i, 
           scale = "none", 
           cluster_rows = F, 
           cluster_cols = F,
           fontsize = 10,
           cellwidth = 15,
           cellheight = 15,
           gaps_col = 5,
           filename = paste0("RNA_Seq_2022/Heatmaps/", i, "_heatmap.pdf"))
  
}

df %>% 
  filter(gene %in% c("Solyc04g039670", "Solyc12g099260", "Solyc05g006520")) %>%
  ggplot(aes(x = genotype, y = normalised_counts, fill = genotype)) +
  geom_boxplot()+
  geom_point(color = "darkgrey")+
  facet_grid(~gene)+
  scale_fill_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none")

ggsave("RNA_Seq_2022/citrate_shuttle_genes.png", width = 5, height = 3)
