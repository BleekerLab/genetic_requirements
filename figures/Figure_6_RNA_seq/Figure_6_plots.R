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
  mutate(genotype = if_else(genotype == "Elite", "Cultivar", "F2-28")) %>%
  mutate(genotype_sample= paste0(genotype, "_", sample))


# Open genes of interest and description
targets <- read.delim("figures/Figure_6_RNA_seq/precursor_genes.tsv") %>% 
  filter(name != "TPS20")
targets2 <- targets %>% 
  filter(target_id %in% df$gene) %>% 
distinct(target_id, .keep_all =T)

# Fuse gene info of the targets and filter them from the dataset
df.targets <- 
  df %>% 
  filter(gene %in% targets$target_id) %>% 
  left_join(., targets, by = c("gene" = "target_id")) %>% 
  mutate(normalised_counts = log(normalised_counts+1, base = 2)) # log scale transformation


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
  select(genotype_sample, gene, name, normalised_counts) %>% 
  pivot_wider(names_from = genotype_sample, values_from = normalised_counts) %>% 
  column_to_rownames(var = "gene")

# Colors for the plot
my_colour = list(
  pathway= c(MEP = "black", MVA = "grey")
)

# Annotation of the rows in the heatmap
annotation_rows <- targets %>% filter(pathway %in% c("MEP", "MVA")) %>% column_to_rownames(var = "target_id") %>% select(pathway) # this makes the rows separate in MEP/MVA
annotation_rows$name <- targets %>% filter(pathway %in% c("MEP", "MVA")) %>% column_to_rownames(var = "target_id") %>% .$name 

pheatmap(mat = df.mep.mva %>% select(-name), 
         color = colorRampPalette(c("white","yellow" ,"red"))(20),
         scale = "none", 
         cluster_rows = F, 
         cluster_cols = F,
         fontsize = 10,
         cellwidth = 15,
         cellheight = 15,
         #annotation_col = annotation_cols, 
         annotation_row = annotation_rows,
         annotation_colors = my_colour,
         gaps_row = 11,
         gaps_col = 5,
        filename = "figures/Figure_6_RNA_seq/MEP_MVA_heatmap_normalised_counts2.pdf")
##########
# Ratios #
##########

df.mean <- df %>% 
  group_by(genotype, gene) %>% 
  summarise(mean_counts = mean(normalised_counts)) %>% 
  filter(gene %in% targets$target_id)

df.ratio <- df.mean %>% 
  pivot_wider(names_from = genotype, values_from = mean_counts) %>% 
  mutate(Log2_ratio = log(P28/Cultivar,2)) %>% 
  left_join(., targets, by = c("gene" = "target_id")) %>% 
  arrange(-P28, -Log2_ratio)

df.ratio %>% 
  pivot_longer(cols = c(Cultivar, P28), names_to = "plant", values_to = "normalised_counts") %>% 
  ggplot(aes(x = log(normalised_counts,2), fill = plant)) + 
  geom_density(alpha = 0.5)+
  facet_wrap(~pathway, ncol = 1)+
  theme_bw()

mat.ratio <-
df.ratio %>% 
  select(gene, Log2_ratio) %>%
  arrange(factor(gene, levels = as.vector(targets %>% 
                                            filter(pathway %in% c("MEP", "MVA")) %>% 
                                            filter(target_id %in% df.targets$gene) %>% 
                                            .$target_id))) %>%
  column_to_rownames("gene")

pheatmap(mat = mat.ratio, 
         color = colorRampPalette(c("blue","white" ,"red"), bias = 1.75)(20),
         scale = "none", 
         cluster_rows = F, 
         cluster_cols = F,
         fontsize = 10,
         cellwidth = 15,
         cellheight = 15,
         #annotation_col = annotation_cols, 
         annotation_row = annotation_rows,
         annotation_colors = my_colour,
         gaps_row = 11,
         filename = "figures/Figure_6_RNA_seq/ratio_heatmap.pdf")

###########
# Boxplot #
###########

df %>% 
  filter(gene %in% targets$target_id) %>%
  left_join(., targets, by = c("gene" = "target_id")) %>%
  mutate(gene = paste(name, gene)) %>%
  ggplot(aes(x = genotype, y = normalised_counts, fill = genotype)) +
  geom_boxplot()+
  geom_point(color = "black")+
  facet_wrap(pathway~gene)+
  scale_fill_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 6)
  )
  
ggsave("figures/Figure_6_RNA_seq/boxplot_all_MVA_MEP.pdf", height = 11, width = 8)

top_genes <- c("Solyc11g010850", 
               "Solyc01g109300", 
               "Solyc11g069380",
               "Solyc04g056390", 
               "Solyc08g005680") 

df %>% 
  filter(gene %in% top_genes) %>%
  left_join(., targets, by = c("gene" = "target_id")) %>%
  mutate(gene = paste(name, gene)) %>%
  ggplot(aes(x = genotype, y = normalised_counts, fill = genotype)) +
  geom_boxplot()+
  geom_point(color = "black")+
  xlab(NULL)+
  ylab("Normalised counts")+
  facet_grid(~gene)+
  scale_fill_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(size = 7),
        axis.text  = element_text(colour = "black", size = 8),
        axis.title.y = element_text(size = 8, colour = "black")
  )

ggsave("figures/Figure_6_RNA_seq/boxplot_top_candidates.pdf", height = 2.25, width = 6.5)


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
  filter(gene %in% c("Solyc04g039670", "Solyc12g099260", "Solyc05g006520", "Solyc12g015860")) %>%
  ggplot(aes(x = genotype, y = normalised_counts, fill = genotype)) +
  geom_boxplot()+
  geom_point(color = "darkgrey")+
  facet_grid(~gene)+
  scale_fill_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none")

ggsave("RNA_Seq_2022/citrate_shuttle_genes.png", width = 5, height = 3)
