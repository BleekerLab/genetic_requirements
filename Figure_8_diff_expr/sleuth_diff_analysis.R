library("tidyverse")
library("sleuth")  #devtools::install_github("pachterlab/sleuth", dependencies = TRUE, upgrade = TRUE, quiet = TRUE)


# Creation of the Sleuth object

kallistos <- list.dirs("results/kallisto/",full.names = TRUE, recursive = FALSE)


samples2condition <- read.delim(file = "samples.tsv", 
                               header = TRUE, 
                               stringsAsFactors = FALSE) %>% 
  dplyr::select("sample", "condition") %>% 
  dplyr::mutate(path = kallistos)


so <- sleuth_prep(samples2condition, full_model = ~ condition, extra_bootstrap_summary = TRUE)

###############################
# Fit model + perform Wald test
###############################

so <- sleuth_fit(obj = so, formula = ~ condition, fit_name = 'full')
so <- sleuth_wt(so, which_beta = "conditionactive", which_model = "full")

all_genes <- sleuth_results(so, test = "conditionactive", test_type = "wt",show_all = TRUE) 
significant_genes <- filter(all_genes, qval <= 0.05)
  
#########################################
# Sleuth live for interactive exploration
#########################################

save(so, file = "kallisto_sleuth_analysis.RData")

sleuth_live(so)