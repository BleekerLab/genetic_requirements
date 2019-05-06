df = read.csv(file = "F2_RNAseq_mapped_genes_abundance_tidy.tsv", sep = '\t', header = TRUE)
str(df)
df$tpm = as.factor(df$tpm)
df$genotype = factor(df$genotype, levels = c("CV", "PI127826", "F1", "F2_151", "F2_411", "F2_445", "F2_28", "F2_73", "F2_127"), 
                        ordered = TRUE)
df %>% filter(target_id == "Solyc10g075090.2.1") %>%
ggplot(., aes(x = genotype,
                                     y = est_counts,
                                     fill = group )) +
  geom_bar(stat = "identity")
  ggsave("Heatmap_ggplot.jpg", device = "jpg", scale = 1, width = 28, height = 28, units = "cm", dpi = 100) + options(device=NULL)

df2 = df[,1;2;3;7]
reshape(dat1, idvar = "name", timevar = "numbers", direction = "wide")

df.parsed <- dplyr::select(df, genotype, group, target_id, tpm)
df.parsed$tpm <- as.numeric(df.parsed$tpm)
df.wide2 = spread(df.parsed, 
                 key = c("genotype","group"),
                 value = "tpm",
                 convert = F,
                 fill = NA)
df.parsed = na.omit(df.parsed)
str(df.parsed)
head(df.parsed)
