

normalise_counts <- function(counts_to_normalise){
  
  # creates a fake design file (obligatory for dds object creation)
  col_data = data.frame(
    sample = colnames(counts_to_normalise)[2:ncol(counts_to_normalise)], 
    condition = c(rep("A",4), rep("B","5"))) # arbitrary conditions (not used)
  
  # creates DESeqDataSet
  dds <- DESeqDataSetFromMatrix(counts_to_normalise[,-1], 
                                colData = colData,
                                design = ~ condition)
  
  # estimate size factors
  dds <- estimateSizeFactors(dds)
  counts_normalised = counts(dds, normalized = TRUE)
  
  # extract the normalised counts
  return(counts_normalised)
}