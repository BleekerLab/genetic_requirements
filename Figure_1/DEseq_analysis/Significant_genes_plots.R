library(tidyverse)
library(ggrepel)

#############################
# Custom theme for plotting #
#############################

my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 6, colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6, colour = "black"), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.text.x = element_text(size=6, colour = "black")
  )+
  theme_bw()

###############
# import data #
###############

df.counts = read.delim("Figure_1/F2_RNAseq_mapped_genes_abundance_tidy.tsv", header = T, stringsAsFactors = T)
df.counts$tpm = as.numeric(df.counts$tpm)
names(df.counts)[3] = "Sequence_ID" #change the same of one column
df.counts$genotype = factor(df.counts$genotype, levels =c("PI127826","F1_hab","F2_73","F2_28","F2_127","CV","F1","F2_151","F2_411","F2_445"), ordered = T)

df.fold.change = na.omit(read.delim("figure_1/DEseq_analysis/DEseq_resuts_Active_vs_Lazy_annotated.tsv", header = T, stringsAsFactors = T, check.names = F))
df.fold.change$Sequence_ID = as.character(df.fold.change$Sequence_ID)
str(df.fold.change)

threshold = df.fold.change$padj < 0.01
df.fold.change$threshold = as.factor(threshold)

p.volcano =
ggplot(df.fold.change, aes(x = log2FoldChange, y = -log(padj)))+
  geom_point(aes(x = log2FoldChange, y = -log(padj), colour = threshold))+
  scale_color_manual(values = c("grey50", "red"))+
  geom_text_repel(aes(label= ifelse(-log(df.fold.change$padj) >8, Sequence_ID, "")), size = 1.5, hjust = 1.2, vjust = 0.5)+
  my.theme

ggsave(file = "Figure_1/Deseq_analysis/plots/Volacano_p_001.svg", plot = p.volcano, height = 4, width = 4)
ggsave(file = "Figure_1/Deseq_analysis/plots/Volacano_p_001.png", plot = p.volcano, height = 4, width = 4)

##########################
# Plot significant genes #
##########################

# Filter on significant genes
df.significant = df.fold.change %>% filter(padj < 0.01) 
df.significant = df.significant[order(df.significant$padj),]

#take top-20 with lowest p-value
df.top.significant = df.significant[1:20,]

# Merge with df containing the acutal counts
df.top.significant.counts = left_join(df.top.significant, df.counts, by = "Sequence_ID") %>% filter(!genotype == "F1_hab" )
df.top.significant.counts$Functional.annotation = as.character(df.top.significant.counts$Functional.annotation)

# Plot

p.bar = 
ggplot(df.top.significant.counts, aes(x = genotype, y = est_counts))+
  geom_bar(aes(x = genotype, y = est_counts, fill = group), stat = "identity")+
  facet_wrap(~Sequence_ID, scale = "free", ncol = 4)+
  my.theme +
  annotate(geom = "text", x = 1, y = 50, hjust = "left", vjust = "middle", # Gene annotation is plottet
            label = df.significant[1:20,]$Functional.annotation, size = 2)+
  scale_fill_manual(values = c("black", "grey"))

ggsave(file = "Figure_1/Deseq_analysis/plots/Top20_most_significant_barplot.svg", plot = p.bar, height = 11, width = 11)

