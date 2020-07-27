library(tidyverse)
library(DESeq2)
library(pheatmap)
library(ggrepel)
library(EnhancedVolcano)


# Load reads
counts = read.table("Figure_8/raw_counts.txt", header = TRUE, sep = "\t", check.names = FALSE) %>%
  mutate(target_id = substr(Geneid, start = 6, stop = 19)) %>% select(target_id, PI127826, Elite_01,PI127826_F1, F2_28, F2_73, F2_127, F2_151, F2_411, F2_445) 
#  pivot_longer(cols = -target_id, names_to = "genotype", values_to = "value")


# Lazy / Active information on the samples
sampleinfo <- data.frame(cbind(colnames(counts[,2:10]),
                               c("active","lazy","lazy","active","active","active","lazy","lazy","lazy")), stringsAsFactors = F)

colnames(sampleinfo) = c("SampleName","Condition")

#shape counts dataset to make it ready for DEseq
counts4DE <- counts
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
res_for_volcano = lfcShrink(DES, contrast = c("Condition","lazy","active"), res=res, type = 'normal') # This is for creating the volcanoplot

res_df = as.data.frame(res) %>% rownames_to_column(., var = "target_id")
res.significant <- res_df %>% filter(padj < 0.05)

# Export result to txt
 write.table(res_df, file = "Figure_8/F2_RNAseq_DEseq_resuts.tsv", sep = "\t", row.names = FALSE)
 write.table(res.significant, file = "Figure_8/F2_RNAseq_DEseq_resuts_significant_genes.tsv", sep = "\t", row.names = FALSE)

 ###############
 # Volcanoplot #
 ###############
 EnhancedVolcano(res_for_volcano,
                 lab = rownames(res),
                 x = 'log2FoldChange',
                 y = 'pvalue',
                 labSize = 3.0)
 
 #######
 # PCA #
 #######
 # Theme for plotting
 my.theme = theme_bw()+
   theme(text = element_text(),
         axis.text.x = element_text(size = 10, colour = "black"),
         axis.text.y = element_text(size = 10, colour = "black")
   )
 
 df.for.pca = rlog(DES) # normalise (log2) the data
 pca = plotPCA(df.for.pca, intgroup="Condition", returnData = TRUE, ntop = NULL) #use the returnData = TRUE argument to enable pca to be imported in ggplot
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
 
 
ggsave(file = "Figure_8/PCA.pdf", plot = p.pca, width = 5, height = 5)
 