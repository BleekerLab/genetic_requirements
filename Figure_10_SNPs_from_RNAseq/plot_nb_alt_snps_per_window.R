library("tidyverse")
library("patchwork")

################################
# Import and separate end/counts
################################
df1 <- read.delim("Figure_10_SNPs_from_RNAseq/counts_from_snps/PI127826.counts.tsv", 
                 stringsAsFactors = F, 
                 col.names = c("chrom","start","end", "count")) 

df2 <- read.delim("Figure_10_SNPs_from_RNAseq/counts_from_snps/F2-28_intersect_PI127826.counts.tsv", 
                  stringsAsFactors = F, 
                  col.names = c("chrom","start","end", "count")) 

df3 <- read.delim("Figure_10_SNPs_from_RNAseq/counts_from_snps/Elite_2020.counts.tsv", 
                  stringsAsFactors = F, 
                  col.names = c("chrom","start","end", "count")) 

###########################################
# Plot number of SNPs per window start site
###########################################

max_of_all_dfs <- map_int(list(df1,df2,df3), function(x) max(x$count)) %>% max()

p_PI127826 <- df1 %>% 
  filter(chrom != "SL2.50ch00") %>% 
  ggplot(., aes(x = start, y = count)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ chrom, scales = "free_x") +
  ggtitle("PI127826 number of SNPs per 1Mb") +
  ylim(0, max_of_all_dfs)

p_F2 <- df2 %>% 
  filter(chrom != "SL2.50ch00") %>% 
  ggplot(., aes(x = start, y = count)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ chrom, scales = "free_x") +
  ggtitle("F2-28 number of SNPs per 1Mb") +
  ylim(0, max_of_all_dfs)

p_Elite <- df3  %>% 
  filter(chrom != "SL2.50ch00") %>% 
  ggplot(., aes(x = start, y = count)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ chrom, scales = "free_x") +
  ggtitle("Elite line number of SNPs per 1Mb") +
  ylim(0, max_of_all_dfs)

p_PI127826 + p_F2 + p_Elite

ggsave(filename = "Figure_10_SNPs_from_RNAseq/plots/nb_snps_per_1Mb_window.pdf", 
       width = 20, 
       height = 8)
