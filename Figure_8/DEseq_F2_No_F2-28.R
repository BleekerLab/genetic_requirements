library(tidyverse)
library(DESeq2)
library(pheatmap)
library(ggrepel)
library(EnhancedVolcano)

##################
# Import dataset #
##################

# Load reads
counts <- read.table("Figure_8/raw_counts.txt", header = TRUE, sep = "\t", check.names = FALSE) %>%
  mutate(target_id = substr(Geneid, start = 6, stop = 19)) %>% 
   select(target_id, PI127826, Elite_01,PI127826_F1, F2_73, F2_127, F2_151, F2_411, F2_445)

# Lazy / Active information on the samples
sampleinfo <- data.frame(cbind(colnames(counts[,2:9]),
                               c("active","lazy","lazy","active","active","lazy","lazy","lazy")), stringsAsFactors = F)

colnames(sampleinfo) <- c("SampleName","Condition")

#Create a df for DESeq2 and add sample names
counts4DE <- counts
colnames(counts4DE)[2:9] <- sampleinfo[,1] # Add samples names 

####################
# Clean-up dataset #
####################

# create function to calculate how many samples have >200 counts per gene
numOverTen <- function(x) {sum(x > 200)} 
ExpressionNum <- apply(counts4DE[2:9], 1, numOverTen)
counts4DE <- counts4DE[which(ExpressionNum > 1),] # Keep genes of which at least 2 samples have > 10 counts

# Remove rows with duplicated gene names
counts4DE <- distinct(counts4DE, target_id, .keep_all = TRUE)  


###################
# DEseq2 Analysis #
###################

DES <- DESeqDataSetFromMatrix(counts4DE, sampleinfo, ~ Condition, tidy = TRUE)
head(DES)
DES <- DESeq(DES, parallel = T) #creates a normalised dataset
# plotDispEsts(DES)
# plotMA(DES, main = "Lazy and Active Differences in Gene Expression") # Overview of differences in expression

# Calculate the results of DEseq2 analysis
res = results(DES, contrast = c("Condition","lazy","active"))
res_for_volcano = lfcShrink(DES, contrast = c("Condition","lazy","active"), res=res, type = 'normal') # This is for creating the volcanoplot

# Import gene annotatoins
annotations <- read.csv(file= "Figure_8/ITAG4.1_descriptions.txt", header = F, sep = "\t") 
annotations <- separate(annotations, V1, sep = " ", c("target_id", "annotation"),extra = "merge") %>%
  mutate(target_id = substr(target_id, start = 1, stop = 14))

# filter significant DE genes from results and add the annotations
res_df <- as.data.frame(res) %>% rownames_to_column(., var = "target_id")
res.significant <- res_df %>% filter(padj < 0.01) %>%
left_join(annotations, by = "target_id")

# Export result to txt files
 write.table(res_df, file = "Figure_8/F2_RNAseq_DEseq_resuts.tsv", sep = "\t", row.names = FALSE)
 write.table(res.significant, file = "Figure_8/F2_RNAseq_DEseq_resuts_significant_genes.tsv", sep = "\t", row.names = FALSE)

 ###############
 # Volcanoplot #
 ###############
 
 p.volcano =
 EnhancedVolcano(res_for_volcano,
                 x = 'log2FoldChange',
                 y = 'pvalue',
                 pCutoff = max(res.significant$pvalue),
                 lab = rownames(res),
                 xlab =  bquote(~Log[2]~ "fold change"),
                 ylab = bquote(~-Log[10]~italic(Pvalue))
 )
 
 ggsave(file = "Figure_8/plots/volcanoplot.pdf", plot = p.volcano, width = 6, height = 6)
 
 
 ##########################################
 # Plotting significantly expressed genes #
 ##########################################
 
 # Extracting normalised counts from DEseq2 dataset
 normalised.counts <- assay(DES) 
 normalised.counts.tidy = 
   normalised.counts %>% 
   data.frame() %>%
   rownames_to_column(var="target_id") %>% 
   pivot_longer(-target_id, names_to = "genotype", values_to = "count") 

 # set conditions per genotype
 con = data.frame(genotype = c("Elite_01", "PI127826_F1", "F2_151", "F2_411", "F2_445",
                             "PI127826", "F2_73", "F2_127"),
                  condition = c("lazy","lazy","lazy","lazy","lazy",
                                "active","active","active"))
 
 # Fuse the active/lazy condition with the main df
normalised.counts.tidy = left_join(normalised.counts.tidy, con, by = "genotype")

# Fuse with the annotations
normalised.counts.tidy = left_join(normalised.counts.tidy, annotations, by = "target_id")

# Set genotypes in proper order
normalised.counts.tidy$genotype = factor(normalised.counts.tidy$genotype, 
                                         levels = c("Elite_01", "PI127826_F1", "F2_151", "F2_411", "F2_445",
                                                    "PI127826", "F2_73", "F2_127"),
                                         ordered = TRUE)
 
 
 # Custom theme for plotting 
my_theme = theme_bw()+
   theme(text = element_text(),
         axis.text.x = element_text(size = 8, colour = "black", angle = 45, hjust = 1),
         axis.text.y = element_text(size = 8, colour = "black"),
         strip.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_rect(),
         panel.background = element_rect(fill = NA, color = "black"),
         strip.text.x = element_text(size=8, colour = "black"))
 
 
# Determine target genes to plot 
 res.significant =  res.significant %>% arrange(padj,-(baseMean))
 diff.top <- res.significant[1:5,1]


 # Barplot per genotype
 p.barplot.top5 = 
   normalised.counts.tidy %>% filter(target_id %in% diff.top) %>%
   ggplot(aes(x = genotype, y = count, fill = condition))+
   geom_bar(stat = "identity")+
   scale_fill_manual(values = c("lazy" = "grey", "active" = "black"))+
   facet_wrap(~target_id,  scale = "free", ncol = 3)+
   labs(x = "Sample" , y = "Gene expression (normalised counts)")+
   my_theme
 
 ggsave(file = "Figure_8/plots/barplot_top5.pdf", plot = p.barplot.top5, width = 7, height = 4)


# Boxplot per condition
 diff.top15 <- res.significant[1:15,1]
 p.boxplot.top15 =
  normalised.counts.tidy %>% filter(target_id %in% diff.top15) %>%
  ggplot(aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0))+
  geom_boxplot()+
  #geom_text_repel(aes(label = genotype)) +
  scale_fill_manual(values = c("lazy" = "grey", "active" = "black"))+
  #geom_text(aes(label = annotation), check_overlap = T, size = 2)+
  facet_wrap(~target_id, scale = "free")+
  labs(x = "Sample" , y = "Gene expression (normalised counts)")+
  my_theme
 ggsave(file = "Figure_8/plots/boxplot_top15.pdf", plot = p.boxplot.top15, width = 10, height = 10)
 

##################################################
# Plot normalised counts distribution (baseMean) #
##################################################

p.distribution =
ggplot(res_df, aes(x = res_df$baseMean))+
   geom_histogram(binwidth = 100, fill = "white", colour = "black")+
   xlim(0,10000)+
   ylim(0,3000)+
   geom_text(aes(x = 7000, y = 2000), label = "binwith = 100 counts")+
   xlab("DEseq2 baseMean")+
   ylab("number of genes")+
   my_theme
 
ggsave(file = "Figure_8/plots/baseMean_distribution.pdf", plot = p.distribution, width = 7, height = 5)
 
 ######################################
 # Principle component analysis (PCA) #
 ######################################

 # Theme for plotting
 my.pca.theme = theme_bw()+
   theme(text = element_text(),
         axis.text.x = element_text(size = 10, colour = "black"),
         axis.text.y = element_text(size = 10, colour = "black")
   )
 
 df.for.pca = rlog(DES) # normalise (log2) the data
 pca = plotPCA(df.for.pca, intgroup="Condition", returnData = TRUE, ntop = nrow(res_df)) #use the returnData = TRUE argument to enable pca to be imported in ggplot
 percentVar <- round(100 * attr(pca, "percentVar"))

 p.pca = 
 ggplot(pca, aes(x = PC1, y = PC2, colour = group))+
   geom_point()+
   geom_point(size=3) +
   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
   ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
   coord_fixed()+
   geom_label_repel(aes(label = row.names(pca), colour = group))+
   my.pca.theme
 
 
ggsave(file = "Figure_8/plots/PCA.pdf", plot = p.pca, width = 5, height = 5)

#######################
# Heatmap of DE genes #
#######################

df.for.heatmap <- 
   log(
   as.data.frame(normalised.counts) %>%
   rownames_to_column(., "target_id") %>%
   filter(target_id %in% res.significant$target_id) %>%
   column_to_rownames(., var = "target_id") + 1
   )

col_order = c("Elite_01", "PI127826_F1", "F2_151", "F2_411", "F2_445",
              "PI127826", "F2_73", "F2_127")
df.for.heatmap <- df.for.heatmap[, col_order]

pheatmap(df.for.heatmap,
         scale = "none", 
         cluster_rows = F, 
         cluster_cols = F,
         fontsize = 6,
         cellwidth = 10,
         cellheight = 5,
         annotation_colors = my_colour,
         gaps_col = 5,
         filename = "Figure_8/plots/heatmap_sig_expressed.pdf"
)
