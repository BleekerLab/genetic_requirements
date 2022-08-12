###########################
# Post-hoc test - Tukey HSD
###########################

library("tidyverse")
library("agricolae")

source("figures/Figure_01_terpenes_type6_CVxPI127826/Fig1A.R")

#########################################
# Figure 1A: 7-epizingiberene ANOVA + HSD
#########################################
p.fig1a

#### log2 transform to comply with normality assumption

anova_results_7epiZ <- 
  df.long %>% 
  filter(metabolite == "7epiZ") %>% 
  with(., aov(log2(value+1) ~ genotype))

hsd_epi <- 
  HSD.test(y = anova_results_7epiZ,
         trt = "genotype", 
         group = TRUE, 
         console = FALSE)

hsd_epi$groups

hsd_epi$groups
# log2(value + 1) groups
# Elite         0.06957234      c
# PI127826      5.01618822     ab
# F1_hab        6.84499227      a
# F1            1.07333221     bc
# LA1777        0.00000000      c

##############################
# Figure 1A: All terpenes ANOVA + HSD
##############################

anova_results_summed_terpenes <- 
  df.long %>% 
  filter(metabolite == "summed_terpenes") %>% 
  mutate(log2_value = log2(value + 1)) %>% 
  dplyr::select(genotype, log2_value) %>% 
  with(., aov(log2_value ~ genotype))

hsd_res_teprenes <- 
  HSD.test(y = anova_results_summed_terpenes,
         trt = "genotype", 
         group = TRUE, 
         console = FALSE)

hsd_res_teprenes$groups

# log2_value groups
# F1_hab     9.111309      a
# PI127826   7.283347      a
# LA1777     6.748268      a
# F1         3.499134      b
# Elite      1.331037      b
