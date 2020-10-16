library("tidyverse")
library("sleuth")  #devtools::install_github("pachterlab/sleuth", dependencies = TRUE, upgrade = TRUE, quiet = TRUE)


# Creation of the Sleuth object

kallistos <- list.dirs("results/kallisto/",full.names = TRUE, recursive = FALSE)


samples2condition <- read.delim(file = "config/samples.tsv", 
                               header = TRUE, 
                               stringsAsFactors = FALSE) %>% 
  dplyr::select("sample", "condition") %>% 
  dplyr::mutate(path = kallistos)


so <- sleuth_prep(samples2condition, full_model = ~ condition, extra_bootstrap_summary = TRUE)

###########
# Fit model
###########

so <- sleuth_fit(obj = so, formula = ~ condition, fit_name = 'full')
so <- sleuth_wt(so, which_beta = "conditionactive", which_model = "full")
sleuth::plot_volcano(obj = so, test = "conditionactive", test_type = "wt",which_model = "full",sig_level = 0.01)

#########################################
# Sleuth live for interactive exploration
#########################################

sleuth_live(so)