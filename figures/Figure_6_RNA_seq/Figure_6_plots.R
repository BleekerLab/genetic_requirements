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
  filter(name != "TPS20") #remove TPS20 from precursor list

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
annotation_rows$name <- targets %>% filter(pathway %in% c("MEP", "MVA")) %>% column_to_rownames(var = "target_id") %>% .$name # To annotate the heatmap with gene-annotations

# Plot the heatmap 
pheatmap(mat = df.mep.mva %>% select(-name), 
         color = colorRampPalette(c("white","yellow" ,"red"))(15), 
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
         gaps_col = 5)#,
      #  filename = "figures/Figure_6_RNA_seq/MEP_MVA_heatmap_normalised_counts2.pdf")

####################################################
# Differential expression (ratio F2-28 / Cultivar) #
####################################################

# Calculate the mean expression over the replicated
df.mean <- df %>% 
  group_by(genotype, gene) %>% 
  summarise(mean_counts = mean(normalised_counts)) %>% 
  filter(gene %in% targets$target_id) %>% 
  left_join(., targets, by = c("gene" = "target_id")) 

# Calculate the ratio of the mean expression levels 
df.ratio <- df.mean %>% 
  pivot_wider(names_from = genotype, values_from = mean_counts) %>% 
  mutate(Log2_ratio = log(`F2-28`/Cultivar,2)) %>% 
  arrange( -`F2-28`)

# Create a dataframe suitable for to make a heatmap of the ratios
mat.ratio <-
df.ratio %>% 
  select(gene, Log2_ratio) %>%
  arrange(factor(gene, levels = as.vector(targets %>% 
                                            filter(pathway %in% c("MEP", "MVA")) %>% 
                                            filter(target_id %in% df.targets$gene) %>% 
                                            .$target_id))) %>%
  column_to_rownames("gene")

# Plot the heatmap 
pheatmap(mat = mat.ratio, 
         color = colorRampPalette(c("blue","white" ,"red"), bias = 1.75)(20),
         scale = "none", 
         cluster_rows = F, 
         cluster_cols = F,
         fontsize = 10,
         cellwidth = 15,
         cellheight = 15,
         annotation_row = annotation_rows,
         annotation_colors = my_colour,
         gaps_row = 11)
        # filename = "figures/Figure_6_RNA_seq/ratio_heatmap.pdf")

##########################################
# Statistics: T- Test F2-28 vs. Cultivar #
##########################################
df.t.test <-
df %>% 
  filter(gene %in% targets$target_id) %>%
  left_join(., targets, by = c("gene" = "target_id")) %>%
  mutate(gene = paste(name, gene)) %>% 
  select(gene, genotype, normalised_counts) %>%
  group_by(genotype, gene) %>% 
  summarise(normalised_counts = list(normalised_counts)) %>% 
  pivot_wider(names_from = genotype, values_from = normalised_counts) %>% 
  group_by(gene) %>% 
  mutate(shapiro_cultivar = shapiro.test(unlist(Cultivar))$p.value,
         shapiro_F2_28 = shapiro.test(unlist(`F2-28`))$p.value,
         p_value = t.test(unlist(Cultivar), unlist(`F2-28`))$p.value,
         t_value = t.test(unlist(Cultivar), unlist(`F2-28`))$statistic)

 # openxlsx::write.xlsx(df.t.test, "figures/Figure_6_RNA_seq/T_text.results.xlsx")

###############################
# Boxplots of precursor genes #
###############################

# Plot all precursor genes for Supplemental Figure S8
df %>% 
  filter(gene %in% targets$target_id) %>%
  left_join(., targets, by = c("gene" = "target_id")) %>%
  mutate(gene = paste(name, gene)) %>%
  ggplot(aes(x = genotype, y = normalised_counts, fill = genotype)) +
  geom_boxplot()+
  geom_point(color = "black")+
  scale_y_continuous(labels = scales::scientific)+
  facet_wrap(pathway~gene, scale = "free", ncol = 5)+
  scale_fill_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white"),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        strip.text = element_text(size = 6)
  )
  
# ggsave("figures/Figure_6_RNA_seq/boxplot_all_MVA_MEP.pdf", height = 11, width = 8)

# Plot top candidates - Figure 6B
top_genes <- c("Solyc11g010850", 
               "Solyc01g109300", 
               "Solyc11g069380",
               "Solyc04g056390", 
               "Solyc08g005680",
               "Solyc12g099900") 

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

# ggsave("figures/Figure_6_RNA_seq/boxplot_top_candidates.pdf", height = 2.25, width = 6.5)

#################################
# Boxplots of known regulators  #
#################################

# SlEOT1 (Solyc02g062400) 
# SlMYC1 (Solyc08g005050)
# SlSCL3 (Solyc12g099900) 

regulator_genes <- data.frame(target_id = 
                               c("Solyc02g062400", 
                                 "Solyc08g005050", 
                                 "Solyc12g099900"),
                             name = c("EOT1",
                                      "MYC1",
                                      "SCL3")
)

df %>% 
  filter(gene %in% regulator_genes$target_id) %>%
  left_join(., regulator_genes, by = c("gene" = "target_id")) %>%
  mutate(gene = paste(name, gene)) %>%
  ggplot(aes(x = genotype, y = normalised_counts, fill = genotype)) +
  geom_boxplot()+
  geom_point(color = "black")+
  xlab(NULL)+
  ylab("Normalised counts")+
  facet_wrap(~gene, scale = "free")+
  scale_fill_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(size = 7),
        axis.text  = element_text(colour = "black", size = 8),
        axis.title.y = element_text(size = 8, colour = "black")
  )

ggsave("figures/Figure_6_RNA_seq/boxplot_regulators.pdf", height = 2.25, width = 5)

# T-Test
  df %>% 
  filter(gene %in% regulator_genes$target_id) %>%
  left_join(., regulator_genes, by = c("gene" = "target_id")) %>%
  mutate(gene = paste(name, gene)) %>% 
  select(gene, genotype, normalised_counts) %>%
  group_by(genotype, gene) %>% 
  summarise(normalised_counts = list(normalised_counts)) %>% 
  pivot_wider(names_from = genotype, values_from = normalised_counts) %>% 
  group_by(gene) %>% 
  mutate(shapiro_cultivar = shapiro.test(unlist(Cultivar))$p.value,
         shapiro_F2_28 = shapiro.test(unlist(`F2-28`))$p.value,
         p_value = t.test(unlist(Cultivar), unlist(`F2-28`))$p.value,
         t_value = t.test(unlist(Cultivar), unlist(`F2-28`))$statistic)
  
  df %>% 
    filter(gene %in% regulator_genes$target_id) %>%
    left_join(., regulator_genes, by = c("gene" = "target_id")) %>% 
    group_by(genotype, name) %>% 
    summarise(avg = mean(normalised_counts)) %>% 
    arrange(name)
  