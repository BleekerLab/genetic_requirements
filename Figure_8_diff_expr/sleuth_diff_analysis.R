library("checkpoint")
checkpoint("2020-01-01")
library("tidyverse")

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
                  extra_bootstrap_summary = TRUE)

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

significant_genes <- filter(all_genes, 
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

#########################################
# Sleuth live for interactive exploration
#########################################
sleuth_live(so)
