###########################
# Post-hoc test - Tukey HSD
###########################

library("tidyverse")
library("agricolae")

source("figures/Figure_01_terpenes_type6_CVxPI127826/Fig1B.R")

#########
# Abaxial
#########

anova_trichome <- 
  df_tidy %>% 
  mutate(log2_density = log2(density_mm2 + 1)) %>% 
  filter(surface == "abaxial") %>% 
  filter(type == "type_VI") %>% 
  with(., aov(log2_density ~ genotype))


hsd_abaxial <- 
  HSD.test(y = anova_trichome,
           trt = "genotype", 
           group = TRUE, 
           console = FALSE)
hsd_abaxial$groups

###### Abaxial side ####
# Elite       0.8242796      a
# F1          1.5627404     ab
# PI127826    1.6339816     ab
# F1_hab      1.7815971     ab
# LA1777      2.1240166      b


#########
# Adaxial
#########

anova_trichome_adaxial <- 
  df_tidy %>% 
  mutate(log2_density = log2(density_mm2 + 1)) %>% 
  filter(surface == "adaxial") %>% 
  filter(type == "type_VI") %>% 
  with(., aov(log2_density ~ genotype))


hsd_aDaxial <- 
  HSD.test(y = anova_trichome,
           trt = "genotype", 
           group = TRUE, 
           console = FALSE)
hsd_aDaxial$groups

# Elite       0.8242796      a
# F1          1.5627404     ab
# PI127826    1.6339816     ab
# F1_hab      1.7815971     ab
# LA1777      2.1240166      b


