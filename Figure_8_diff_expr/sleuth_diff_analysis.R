library("checkpoint")
checkpoint("2020-01-01")
library("tidyr")
library("ggplot2")
library("dplyr")
library("tibble")
library("ggrepel")
source("Figure_8_diff_expr/mypca.R")

# v0.29.0 (because 0.30.0 has bugs)
# devtools::install_github(repo = "pachterlab/sleuth",
#                         ref = "2cbb28749875eaf7dcafb30f9a795313c481fcb1", 
#                           dependencies = TRUE, 
#                           upgrade = TRUE, 
#                           quiet = TRUE,
#                           force = TRUE)
library("sleuth")  

###############################
# Creation of the Sleuth object
###############################
kallistos <- list.dirs("Supplemental_data_RNA-seq/kallisto",
                       full.names = TRUE, 
                       recursive = FALSE) 
kallistos_df <- data.frame(path = kallistos, 
                           sample = basename(kallistos)) 

samples2condition <- read.delim(file = "Supplemental_data_RNA-seq/config/samples.tsv", 
                               header = TRUE, 
                               stringsAsFactors = FALSE) %>% 
  dplyr::select("sample", "condition") %>% 
  inner_join(x = ., y = kallistos_df, by = "sample") %>% 
  mutate(path = as.character(path))
  
so <- sleuth_prep(samples2condition, 
                  full_model = ~ condition, 
                  extra_bootstrap_summary = FALSE)

###############################
# Fit model + perform Wald test
###############################

so <- sleuth_fit(obj = so, 
                 formula = ~ condition, 
                 fit_name = 'full')

so <- sleuth_wt(so, 
                which_beta = "conditionactive", 
                which_model = "full")

all_genes_wald <- sleuth_results(so, 
                            test = "conditionactive", 
                            test_type = "wt",
                            show_all = TRUE) 

significant_genes_wald <- filter(all_genes_wald, 
                            qval <= 0.05) %>% 
  rename("gene" = "target_id")

#####################################
# Fit model + perform likelihood test
#####################################

so <- sleuth_fit(so, ~ condition, fit_name = "full")
so <- sleuth_fit(so, ~ 1, fit_name = "reduced")
so <- sleuth_lrt(so, null_model = "reduced", alt_model = "full")

all_genes_lrt <- sleuth_results(so, 
                            test = "reduced:full", 
                            test_type = "lrt",
                            show_all = TRUE) 

significant_genes_lrt <- filter(all_genes_lrt, 
                            qval <= 0.05) %>% 
  rename("gene" = "target_id")

##############################################
# Annotate significant genes and save R object
##############################################

gene_annotations <- read.csv(file = "Figure_8_diff_expr/ITAG4.1_descriptions.csv", 
                              stringsAsFactors = FALSE, 
                              header = TRUE)

significant_genes_annotated <- inner_join(significant_genes_wald, 
                                          gene_annotations, 
                                          by = "gene")


save(so, file = "Figure_8_diff_expr/kallisto_sleuth_analysis.RData")

###########################
# Extract normalised counts
###########################

scaled_counts <- kallisto_table(so, use_filtered = TRUE, normalized = TRUE) %>% 
  rename("gene" = "target_id") %>% 
  select(- tpm, - eff_len) %>% 
  pivot_wider(id_cols = "gene", names_from = "sample", values_from = "est_counts") %>% 
  column_to_rownames("gene") %>% 
  t(.) # for compatibility with mypca() function

write.csv(x = scaled_counts, 
          file = "Figure_8_diff_expr/scaled_counts.csv", 
          row.names = TRUE, 
          quote = FALSE)

###########################
# PCA analysis
###########################  
pca_results <- mypca(scaled_counts)

### Scree plot
df_explained_variance <- data.frame(
  exp_var = pca_results$explained_var$exp_var
) %>% 
  rownames_to_column("PC") %>% 
  mutate(PC = factor(PC,levels = unique(PC)))

### Sample score plot
scores <- pca_results$scores %>% 
  rownames_to_column("sample") %>% 
  inner_join(., y = samples2condition, by = "sample") %>% 
  select(- path) 

pc1_vs_pc2 <- ggplot(scores, aes(label = sample)) + 
  geom_point(aes(x = PC1, y = PC2, shape = condition, col = condition), size = 3) + 
  xlab(paste0('PC1 (',df_explained_variance[1,2],'% explained variance)')) + 
  ylab(paste0('PC2 (',df_explained_variance[2,2],'% explained variance)')) + 
  ggtitle('PCA score plot: PC1 versus PC2') +
  geom_text_repel(aes(x = PC1, y = PC2))
pc1_vs_pc2 

pc2_vs_pc3 <- ggplot(scores, aes(label = sample)) + 
  geom_point(aes(x = PC2, y = PC3, shape = condition, col = condition), size = 3) + 
  xlab(paste0('PC2 (',df_explained_variance[2,2],'% explained variance)')) + 
  ylab(paste0('PC3 (',df_explained_variance[3,2],'% explained variance)')) + 
  ggtitle('PCA score plot: PC2 versus PC3') +
  geom_text_repel(aes(x = PC2, y = PC3))
pc2_vs_pc3

ggsave(filename = "Figure_8_diff_expr/PCA_pc1_vs_pc2.png", plot = pc1_vs_pc2)
ggsave(filename = "Figure_8_diff_expr/PCA_pc2_vs_pc3.png", plot = pc2_vs_pc3)

#########################################
# Sleuth live for interactive exploration
#########################################
sleuth_live(so)
