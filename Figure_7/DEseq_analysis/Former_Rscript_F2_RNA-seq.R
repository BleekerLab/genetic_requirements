library(ggplot2)
library(DESeq2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)
library(plyr)
library(ggfortify)

counts = read.table("/Users/Ruy/Documents/Zingiberene & Derivatives/Large F2 EZ/High_Low_producers/RNA-seq F2/Abundance data/F2_RNAseq_est_counts_wide.tsv", header = T, sep = "\t")
annotations = read.table("/Users/Ruy/Documents/Zingiberene & Derivatives/Large F2 EZ/High_Low_producers/RNA-seq F2/Abundance data/Annotations.tsv", header = T, sep = "\t", fill = T) %>% dplyr::rename(., Sequence_ID = target_id)

counts = select(counts, -F1_hab)

sampleinfo <- data.frame(cbind(colnames(counts[,2:10]),
                               c("lazy","active","lazy","active","lazy","lazy","active","active","lazy")), stringsAsFactors = F)


colnames(sampleinfo) = c("SampleName","Condition")

#shape counts dataset to make it ready for DEseq
counts4DE <- as.data.frame(counts)
counts4DE[2:10] = lapply(counts4DE[2:10], as.integer)
colnames(counts4DE)[2:10] <- sampleinfo[,1]
head(counts4DE)

##Clean-up dataset 
numOverTen <- function(x) {sum(x > 10)}
ExpressionNum <- apply(counts4DE[2:10], 1, numOverTen)
counts4DE <- counts4DE[which(ExpressionNum > 3),] #throw away data where more less than 3 samples have less than 10 counts 
counts4DE = dplyr::rename(counts4DE, Sequence_ID = target_id)
counts4DE = distinct(counts4DE, Sequence_ID, .keep_all = TRUE)

#DEses analysis 
str(counts4DE)
DES <- DESeqDataSetFromMatrix(counts4DE, sampleinfo, ~ Condition, tidy = TRUE)
head(DES)
DES <- DESeq(DES, parallel = T) #creates a normalised dataset
plotDispEsts(DES)
plotMA(DES, main = "Lazy and Active Differences in Gene Expression") #Differenitally expressed genes

#Export results + add functional annotation to gene identifiers + Add (non-normalised) counts from samples
res = results(DES, contrast = c("Condition","lazy","active"))
res = as.data.frame(res) %>% rownames_to_column(., var = "Sequence_ID")
head(res)
results.annotated = left_join(res, annotations, by = "Sequence_ID") %>%  select(.,c("Sequence_ID", "Functional.annotation","baseMean","log2FoldChange", "lfcSE", "stat","pvalue","padj")) %>% arrange(.,log2FoldChange, padj)
results.annotated = left_join(results.annotated, counts4DE, by = "Sequence_ID")
write.table(res, file = "F2_RNAseq_DEseq_resuts_RK.tsv", sep = "\t", row.names = FALSE)

write.table((results.annotated %>% select(., Sequence_ID, log2FoldChange)), file = "F2_results_foldChange.tsv", sep = "\t", row.names = FALSE)

#Normalise and scale data + PCA 

  #order the results by p-value (padj)
resordered <- res[order(res$padj),]
rownames(resordered) = resordered$Sequence_ID
head(resordered)
  #remove NAs
resordered <- resordered[!is.na(resordered$padj),]
dim(resordered)

  #Normalise
row.names(counts4DE) = counts4DE$Sequence_ID
  #Scale and transform data for PCA
scaled.counts <- scale(log(counts4DE[2:ncol(counts4DE)]+1))
head(scaled.counts)

  #perform and plot PCA
PCD <- prcomp(t(scaled.counts), center = T)
plot(PCD$x, pch = 20, col = as.factor(sampleinfo$Condition)
legend("top",legend = c("Lazy","Active"), pch = 20, col = c(1,2))
text(PCD$x, labels = sampleinfo$SampleName, pos = 1))

############################
# Heatmap of target groups #
############################

FCf <- scaled.counts[which(rownames(scaled.counts) %in% rownames(resordered)[1:50]),]
dim(FCf)
pheatmap(FCf, filename = "heatmap_top_50_diff.png")

### EXTRACT SPECIFIC GENE SETS ###

#ABC transporters
ABC = results.annotated[grep("ABC transporter", results.annotated$Functional.annotation),] 
ABC = ABC[order(ABC$padj),]
ABC =  ABC %>% filter(., padj < 0.05) # To display significant genes only
ABC.scaled<- scaled.counts[which(rownames(scaled.counts) %in% ABC$Sequence_ID),]
ABC.scaled = ABC.scaled[,c("PI127826","F1_hab","F2_73","F2_28","F2_127","CV","F1","F2_151","F2_411","F2_445")]
head(ABC.scaled)
pheatmap(ABC.scaled, cluster_cols = F, filename = "ABC_transporters_heatmap_scaled.png")
write.table(ABC, file = "F2_RNAseq_ABC_transporters_results.tsv", sep = "\t", row.names = FALSE)


#P450's
P450 = results.annotated[grep("P450", results.annotated$Functional.annotation),] 
P450 = P450[order(P450$padj),]
P450.scaled<- scaled.counts[which(rownames(scaled.counts) %in% P450$Sequence_ID),]
P450.scaled = P450.scaled[,c("PI127826","F1_hab","F2_73","F2_28","F2_127","CV","F1","F2_151","F2_411","F2_445")]
head(P450.scaled)
pheatmap(P450.scaled, cluster_cols = F, filename = "P450_heatmap_scaled.png")
write.table(P450, file = "F2_RNAseq_P450_results.tsv", sep = "\t", row.names = FALSE)

#Transcription Factors
TF = results.annotated[grep("transcription factor", results.annotated$Functional.annotation),] 
TF = TF[order(TF$padj),]
# TF = TF %>% filter(., padj < 0.05) # For significant genes only
TF.scaled<- scaled.counts[which(rownames(scaled.counts) %in% TF$Sequence_ID),]
TF.scaled = TF.scaled[,c("PI127826","F1_hab","F2_73","F2_28","F2_127","CV","F1","F2_151","F2_411","F2_445")]
head(TF.scaled)
pheatmap(TF.scaled, cluster_cols = F, filename = "TF_heatmap_scaled.png")
write.table(TF, file = "F2_RNAseq_TF_results.tsv", sep = "\t", row.names = FALSE)


#Precursor genes
pr.genes = counts = read.table("precursor_target_genes.tsv", header = T, sep = "\t")
PR = results.annotated[which(results.annotated$Sequence_ID %in% pr.genes$id),]
PR = PR[order(PR$padj),]
#PR = PR %>% filter(., padj < 0.05) # For significant genes only
PR.scaled = scaled.counts[which(rownames(scaled.counts) %in% PR$Sequence_ID),]
PR.scaled = PR.scaled[,c("PI127826","F2_73","F2_28","F2_127","CV","F1","F2_151","F2_411","F2_445")]
pheatmap(PR.scaled, cluster_cols = T,
         annotation_row = pr.genes)
write.table(PR, file = "F2_RNAseq_Precursor_genes_results.tsv", sep = "\t", row.names = FALSE)


###################
# Plot Precursors #
###################

precursor.counts <- counts4DE[which(counts4DE$Sequence_ID %in% rownames(pr.genes)),]
precursor.counts.long = gather(
  data = precursor.counts,
  key = "genotype",
  value = counts, -Sequence_ID)

precursor.counts.long = left_join(precursor.counts.long, sampleinfo, by = c("genotype" = "SampleName"))
precursor.counts.long = left_join(precursor.counts.long, pr.genes, by = c("Sequence_ID" = "id"))

precursor.counts.long$genotype = factor(precursor.counts.long$genotype, levels =c("PI127826","F2_73","F2_28","F2_127","CV","F1","F2_151","F2_411","F2_445"), ordered = T)

precursor.counts.long  %>% filter(., pathway == "MEP") %>% group_by(., name, genotype, Condition) %>% dplyr::summarise(counts = mean(counts), sd = sd(counts)) %>%
ggplot(., aes(x = genotype, y = counts, fill = Condition))+
  geom_bar(stat = "identity") +
  facet_wrap(~name, scale = "free") +
  theme(text = element_text(family = "Arial", color = "black", size = 8),
        axis.text.x = element_text(size = 5, angle = 90),
        strip.text = element_text(size=5, colour = "black"))

precursor.counts.long  %>% filter(., pathway == "MVA") %>% group_by(., name, genotype, Condition) %>% dplyr::summarise(counts = mean(counts), sd = sd(counts)) %>%
  ggplot(., aes(x = genotype, y = counts, fill = Condition))+
  geom_bar(stat = "identity") +
  facet_wrap(~name, scale = "free") +
  theme(text = element_text(family = "Arial", color = "black", size = 8),
        axis.text.x = element_text(size = 5, angle = 90),
        strip.text = element_text(size=5, colour = "black"))



########################
# Annotate counts file #
########################

counts4DE = counts4DE %>% rename(Sequence_ID = target_id)
counts4DE.annotated = left_join(counts4DE, annotations, by = "Sequence_ID")
counts4DE.annotated <- counts4DE.annotated[, c("Sequence_ID", "Functional.annotation","CV","F1","F1_hab","F2_127","F2_151","F2_28","F2_411","F2_445","F2_73","PI127826")]

gene.annotations = annotations[which(annotations[,1] %in% rownames(resordered)[1:50]),]
#write.table(A[,1:2], file= "top_50_differentially_expressed_genes_annotated.tsv", sep = "\t")
Top_normalised = data.frame(normalised_counts[which(rownames(normalised_counts) %in% rownames(resordered)[1:50]),]) %>% rownames_to_column(., var = "Sequence.ID")
                 

top = left_join(A, Top_normalised, by = "Sequence.ID")

#Extract top-differentially expressed genes from non-normalised table and plot them

top_genes <- counts4DE[which(counts4DE$Sequence_ID %in% resordered$Sequence_ID[1:10]),]
top_genes.long = gather(
  data = top_genes,
  key = "genotype",
  value = counts, -Sequence_ID)

top_genes.long$genotype = factor(top_genes.long$genotype, levels = c("CV","PI127826","F1","F1_hab","F2_28","F2_73","F2_127","F2_151","F2_411","F2_445"), ordered = T)
 
ggplot(top_genes.long, aes(x = genotype, y = counts))+
  geom_bar(stat = "identity") +
  facet_wrap(~Sequence_ID, scale = "free_y") +
  theme(text = element_text(family = "Arial", color = "black", size = 8),
        axis.text.x = element_text(size = 10, angle = 90))
 
ggsave("Top_10_counts.png", device = "png")
  
