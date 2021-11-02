library("tidyverse")

################################
# Import and separate end/counts
################################
df <- read.delim("Figure_10_SNPs_from_RNAseq/F2-28_snp_counts.tsv", 
                 stringsAsFactors = F, 
                 col.names = c("chrom","start","end", "count")) 

###########################################
# Plot number of SNPs per window start site
###########################################

df %>% 
  filter(chrom != "SL2.50ch00") %>% 
  ggplot(., aes(x = start, y = count)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ chrom, scales = "free_x")

ggsave(filename = "Figure_10_SNPs_from_RNAseq/F2-28_number_of_snps_per_1Mb_window.pdf", width = 12, height = 7)
