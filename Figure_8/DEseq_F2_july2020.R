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

# Add annotations to the results dataframe
res_df = as.data.frame(res) %>% rownames_to_column(., var = "target_id")
annotations <- read.table(file= "Figure_8/ITAG4.1_descriptions.txt", header = F, sep = "\t") 
annotations = separate(annotations, V1, sep = " ", c("target_id", "annotation"),extra = "merge") %>%
  mutate(target_id = substr(target_id, start = 1, stop = 14))

res.significant <- res_df %>% filter(padj < 0.05)
res.significant <- left_join(res.significant, annotations, by = "target_id")

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
 
 ######################
 # Barplot sig. genes #
 ######################
 
 normalised.counts <- assay(DES) 
 normalised.counts.tidy = 
   normalised.counts %>% 
   data.frame() %>%
   rownames_to_column(var="target_id") %>% 
   pivot_longer(-target_id, names_to = "genotype", values_to = "count") 

 
 con = data.frame(genotype = c("Elite_01", "PI127826_F1", "F2_151", "F2_411", "F2_445",
                             "PI127826", "F2_28", "F2_73", "F2_127"),
                  condition = c("lazy","lazy","lazy","lazy","lazy",
                                "active","active","active","active"))
 
 # Fuse the active/lazy condition with the main df
normalised.counts.tidy = left_join(normalised.counts.tidy, con, by = "genotype")
normalised.counts.tidy = left_join(normalised.counts.tidy, annotations, by = "target_id")
normalised.counts.tidy$genotype = factor(normalised.counts.tidy$genotype, 
                                         levels = c("Elite_01", "PI127826_F1", "F2_151", "F2_411", "F2_445",
                                                    "PI127826", "F2_28", "F2_73", "F2_127"),
                                         ordered = TRUE)
 
 # Determin target genes 
 res.significant =  res.significant %>% arrange(res.significant$padj)
 diff.top10 <- res.significant[1:15,1]
 
 ######################
 # Theme for plotting #
 ######################
 my_theme = theme_bw()+
   theme(text = element_text(),
         axis.text.x = element_text(size = 8, colour = "black", angle = 45, hjust = 1),
         axis.text.y = element_text(size = 8, colour = "black"),
         strip.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_rect(),
         panel.background = element_rect(fill = NA, color = "black"),
         strip.text.x = element_text(size=8, colour = "black"))

# Barplot per genotype
 
normalised.counts.tidy %>% filter(target_id %in% diff.top10) %>%
   ggplot(aes(x = genotype, y = count, fill = condition))+
   geom_bar(stat = "identity")+
   scale_fill_manual(values = c("lazy" = "grey", "active" = "black"))+
   facet_wrap(~target_id, scale = "free")+
   labs(x = "Sample" , y = "Gene expression (counts)")+
   my_theme

# Boxplot per condition
normalised.counts.tidy %>% filter(target_id %in% diff.top10) %>%
  ggplot(aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0))+
  geom_boxplot()+
  #geom_text_repel(aes(label = genotype)) +
  scale_fill_manual(values = c("lazy" = "grey", "active" = "black"))+
  facet_wrap(~target_id, scale = "free")+
  labs(x = "Sample" , y = "Gene expression (counts)")+
  my_theme
 
 
 
 
 
 
 
 
 d <- plotCounts(DES, gene= res$padj, intgroup="Condition", 
                 returnData=TRUE)
 
 ggplot(d, aes(x=Condition, y=count)) + 
   geom_point(position=position_jitter(w=0.1,h=0))+
   geom_text_repel(aes(label = rownames(d)))  
 
 
 
 
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
 