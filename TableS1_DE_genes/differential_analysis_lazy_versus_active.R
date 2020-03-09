#########
# Library
#########

## Checkpoint library manager
if (! "checkpoint" %in% installed.packages()){
  install.packages("checkpoint")
}
library("checkpoint")
checkpoint("2020-01-01")

## remotes 
library("remotes")


## Sleuth
# rdhf5 is a dependency
if (!requireNamespace("BiocManager", quietly = TRUE)) # Bioconductor package manager
  install.packages("BiocManager")
if (! "rhdf5" %in% installed.packages()){
  BiocManager::install("rhdf5")
}

if (! "sleuth" %in% installed.packages()){
  remotes::install_github("pachterlab/sleuth")
}

suppressPackageStartupMessages(library("rhdf5"))
suppressPackageStartupMessages(library("sleuth"))
suppressPackageStartupMessages(library("tidyverse"))

########################################
# Sleuth object (so) object construction
########################################

sample_to_condition <- read.delim("TableS1_DE_genes/samples.tsv", header = T, stringsAsFactors = F)


so <- sleuth_prep(sample_to_covariates = sample_to_condition, extra_bootstrap_summary = TRUE)

# fit full model: the trichome column distinguish between active from lazy trichomes
so <- sleuth_fit(so, ~ trichome, "full")

# fit the null model: the experimental design does not explain anything
so <- sleuth_fit(so, ~ 1, "null")

########################################
# Extract differentially expressed genes
########################################

so <- sleuth_wt(so, 
                which_beta = "trichomeactive", # specify a character string of denoting which grouping to test. (test against trichome lazy)
                which_model = "full")

# Col explanations: 
# https://pachterlab.github.io/sleuth/docs/sleuth_results.html
diff_results <- sleuth_results(so, 
                               test = "trichomeactive", 
                               test_type = "wt", 
                               which_model = "full", 
                               show_all = TRUE)





# add annotations
annotations <- read.delim("TableS1_DE_genes/annotations.tsv", header = T, col.names = c("target_id", "annotation"), stringsAsFactors = F)

# final result table
diff_results_with_annotations <- left_join(diff_results, annotations, by = "target_id")
write_tsv(diff_results_with_annotations, path = "TableS1_DE_genes/differentials.tsv", col_names = T)


################
# PCA analysis
################

# Extract PC variances retained by percentage 
g <- plot_pc_variance(so, use_filtered = TRUE, units = 'est_counts', pca_number = 5) 
pc1_variance <- g$data$var[1]
pc2_variance <- g$data$var[2]

# pca plot
p_pca <- plot_pca(so, 
                  color_by = "trichome",
                  use_filtered = TRUE, 
                  pc_x = 1, 
                  pc_y = 2,
                  text_labels = TRUE) + 
  labs(x = paste("PC1:",
                 round(pc1_variance,digits = 0),
                 "% of variance explained",
                 sep = " "),
       y = paste("PC2:",
                 round(pc2_variance,digits = 0),
                 "% of variance explained",
                 sep = " "))
p_pca

ggsave(filename = "TableS1_DE_genes/plot_pca.png", p_pca, width = 7, height = 5)

# Loadings
p_loadings <- plot_loadings(so, use_filtered = TRUE, pc_input = 1, pc_count = 10)
p_loadings

ggsave(filename = "TableS1_DE_genes/plot_loadings.png", p_loadings, width = 5, height = 7)

##############
# Volcano plot
#############
p_volcano <- plot_volcano(so, test = "trichomeactive", sig_color = "blue", sig_level = 0.01)
p_volcano

ggsave(filename = "TableS1_DE_genes/plot_volcano.png", p_volcano, width = 5, height = 7)