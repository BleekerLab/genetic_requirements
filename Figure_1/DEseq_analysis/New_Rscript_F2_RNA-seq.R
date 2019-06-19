library(ggplot2)
library(DESeq2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)
library(plyr)

counts = read.table("F2_RNAseq_est_counts_wide.tsv", header = T, sep = "\t")
annotations = read.table("Annotations_solyc_numbers.txt", header = T, sep = "\t", fill = T)

sampleinfo <- data.frame(cbind(c("CV","F1","F1_hab","F2_127","F2_151","F2_28","F2_411","F2_445","F2_73","PI127826"),
                               c("lazy","lazy","active","active","lazy","active","lazy","lazy","active","active")), stringsAsFactors = F)

colnames(sampleinfo) = c("SampleName","Condition")

#shape counts dataset to make it ready for DEseq
counts4DE <- counts
counts4DE[2:11] = lapply(counts4DE[2:11], as.integer)
colnames(counts4DE)[2:11] <- sampleinfo[,1]
head(counts4DE)

##Clean-up dataset 
numOverTen <- function(x) {sum(x > 10)}
ExpressionNum <- apply(counts4DE[2:11], 1, numOverTen)
counts4DE <- counts4DE[which(ExpressionNum > 2),] #throw away data where cells have less that 10 cells
counts4DE = rename(counts4DE, Sequence_ID = target_id)
#DEses analysis 
str(counts4DE)
DES <- DESeqDataSetFromMatrix(counts4DE, sampleinfo, ~ Condition, tidy = TRUE)
head(DES)
DES <- DESeq(DES, parallel = T) #creates a normalised dataset
plotDispEsts(DES)
plotMA(res, main = "Lazy and Active Differences in Gene Expression") #Differenitally expressed genes

#Export results + add functional annotation to gene identifiers + Add (non-normalised) counts from samples
res = results(DES, contrast = c("Condition","lazy","active"))
res = as.data.frame(res) %>% rownames_to_column(., var = "Sequence_ID")
head(res)
results.annotated = left_join(res, annotations, by = "Sequence_ID") %>% .[, c("Sequence_ID", "Functional.annotation","baseMean","log2FoldChange", "lfcSE", "stat","pvalue","padj")] %>% arrange(.,log2FoldChange, padj)
results.annotated = left_join(results.annotated, counts4DE, by = "Sequence_ID")
write.table(results.annotated, file = "F2_RNAseq_DEseq_resuts_log2FoldChange.tsv", sep = "\t", row.names = FALSE)

#Normalise and scale data + PCA 
  #order the results by p-value (padj)
resordered <- res[order(res$padj),]
resordered = column_to_rownames(resordered, var = "Sequence_ID")
rownames(resordered) = resor$Sequence_ID
head(resordered)
  #remove NAs
resordered <- resordered[!is.na(resordered$padj),]
dim(resordered)

  #Normalise
normalised_counts <- counts(DES, normalized=TRUE)
head(normalised_counts)

  #Scale and transform data for PCA
scaled.counts <- t(scale(t(normalised_counts)))
head(scaled.counts)

  #perform and plot PCA
PCD <- prcomp(t(scaled.counts), center = T)
plot(PCD$x, pch = 20, col = as.factor(sampleinfo$Condition))
legend("topleft",legend = c("Lazy","Active"), pch = 20, col = c(1,2))

  #Heatmap top 50 diff_expressed genes
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
PR.scaled = PR.scaled[,c("PI127826","F1_hab","F2_73","F2_28","F2_127","CV","F1","F2_151","F2_411","F2_445")]
pheatmap(PR.scaled, cluster_cols = F, filename = "Precursor_genes_significant_heatmap_scaled.png")
write.table(PR, file = "F2_RNAseq_Precursor_genes_results.tsv", sep = "\t", row.names = FALSE)

  #Plot precursor genes
precursor.counts <- counts4DE[which(counts4DE$Sequence_ID %in% pr.genes$id),]
precursor.counts.long = gather(
  data = precursor.counts,
  key = "genotype",
  value = counts, -Sequence_ID)

precursor.counts.long$genotype = factor(precursor.counts.long$genotype, levels =c("PI127826","F1_hab","F2_73","F2_28","F2_127","CV","F1","F2_151","F2_411","F2_445"), ordered = T)
str(precursor.counts.long)


ggplot(precursor.counts.long, aes(x = genotype, y = counts))+
  geom_bar(stat = "identity") +
  facet_wrap(~Sequence_ID, scale = "free_y") +
  theme(text = element_text(family = "Arial", color = "black", size = 8),
        axis.text.x = element_text(size = 5, angle = 90),
        strip.text = element_text(size=5, colour = "black"))

ggsave("Precurs_genes_counts_plot.png", device = "png")




#Annotate the file
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
  
