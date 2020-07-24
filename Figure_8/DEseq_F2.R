library(ggplot2)
library(DESeq2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)
library(plyr)
library(ggrepel)
library(EnhancedVolcano)


# Load reads
counts = read.table("Figure_7/abundance_tidy.tsv", header = T, sep = "\t") %>% 
  mutate(target_id = substr(target_id, start = 1, stop = 14)) %>%
  filter(sample %in% c("Elite_01", "PI127826_F1", "F2-151", "F2-411", "F2-445",
                       "PI127826", "F2-28", "F2-73", "F2-127"))

# Make data "wide" format
counts = as.data.frame(pivot_wider(counts, names_from = sample, values_from = est_counts))

# Lazy / Active information on the samples
sampleinfo <- data.frame(cbind(colnames(counts[,2:10]),
                               c("lazy","active","lazy","active","lazy","lazy","active","active","lazy")), stringsAsFactors = F)

colnames(sampleinfo) = c("SampleName","Condition")

#shape counts dataset to make it ready for DEseq
counts4DE <- counts
counts4DE[2:10] = lapply(counts4DE[2:10], as.integer) # change floating numbers to intergers 
colnames(counts4DE)[2:10] <- sampleinfo[,1] # Add samples names 

##Clean-up dataset 
numOverTen <- function(x) {sum(x > 1)} 
ExpressionNum <- apply(counts4DE[2:10], 1, numOverTen) 
counts4DE <- counts4DE[which(ExpressionNum > 1),] #throw away data where more less than 3 samples have less than 10 counts 
counts4DE <- distinct(counts4DE, target_id, .keep_all = TRUE) # Remove rows with duplicated gene names 


###################
# DEseq2 Analysis #
###################

DES <- DESeqDataSetFromMatrix(counts4DE, sampleinfo, ~ Condition, tidy = TRUE)
head(DES)
DES <- DESeq(DES, parallel = T) #creates a normalised dataset
plotDispEsts(DES)
plotMA(DES, main = "Lazy and Active Differences in Gene Expression") #Differenitally expressed genes

# Export results
res = results(DES, contrast = c("Condition","lazy","active"))
res = lfcShrink(DES, contrast = c("Condition","lazy","active"), res=res, type = 'normal') # This is for creating the volcanoplot

res = as.data.frame(res) %>% rownames_to_column(., var = "target_id")
res.significant <- res %>% filter(padj < 0.05)

# Export result to txt
 write.table(res, file = "Figure_7/DEseq_analysis/F2_RNAseq_DEseq_resuts.tsv", sep = "\t", row.names = FALSE)
 write.table(res.significant, file = "Figure_7/DEseq_analysis/F2_RNAseq_DEseq_resuts_significant_genes.tsv", sep = "\t", row.names = FALSE)

 ###############
 # Volcanoplot #
 ###############
 EnhancedVolcano(res,
                 lab = rownames(res),
                 x = 'log2FoldChange',
                 y = 'pvalue',
                 xlim = c(-3, 3),
                 labSize = 5.0)
 
 #######
 # PCA #
 #######
 
 # Theme for plotting
 my.theme = theme_bw()+
   theme(text = element_text(),
         axis.text.x = element_text(size = 10, colour = "black"),
         axis.text.y = element_text(size = 10, colour = "black")
   )
 
 
 pca = plotPCA(vsd, intgroup=c("Condition"), returnData = TRUE)
 percentVar <- round(100 * attr(pca, "percentVar"))

 p.pca = 
 ggplot(pca, aes(x = PC1, y = PC2, colour = group))+
   geom_point()+
   geom_point(size=3) +
   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
   ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
   coord_fixed()+
   geom_label_repel(aes(label = row.names(pca), colour = group))+
   my.theme
 
 
ggsave(file = "Figure_7/DEseq_analysis/PCA.pdf", plot = p.pca, width = 5, height = 5)
 