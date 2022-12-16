###########################
# Post-hoc test - Tukey HSD
###########################

library("tidyverse")
library("agricolae")

source("figures/Figure_01_terpenes_type6_CVxPI127826/Fig1B.R")


#################################
# Summed Adaxial + Abaxial side #
#################################

anova_trichome <- 
  df_tidy %>%
  group_by(genotype, plant, leafdisc, type) %>% 
  dplyr::summarise(sum_density = sum(density_mm2)) %>% 
  mutate(log2_density = log2(sum_density + 1)) %>% 
  filter(type == "type_VI") %>% 
  with(., aov(log2_density ~ genotype))


hsd <- 
  HSD.test(y = anova_trichome,
           trt = "genotype", 
           group = TRUE, 
           console = FALSE)
hsd$groups

# TYPE-VI
# log2_density groups
# LA1777       3.255455      a
# F1_hab       2.800878      a
# PI127826     2.594626     ab
# F1           2.399008     ab
# Elite        1.394708      b


# Type-I/IV
anova_trichome <- 
  df_tidy %>%
  group_by(genotype, plant, leafdisc, type) %>% 
  dplyr::summarise(sum_density = sum(density_mm2)) %>% 
  mutate(log2_density = log2(sum_density + 1)) %>% 
  filter(type == "type_I_IV") %>% 
  with(., aov(log2_density ~ genotype))


hsd <- 
  HSD.test(y = anova_trichome,
           trt = "genotype", 
           group = TRUE, 
           console = FALSE)
hsd$groups

# TYPE-I/IV
# LA1777     5.62368070      a
# PI127826   4.52725558      b
# F1_hab     3.85195444      b
# F1         1.95706983      c
# Elite      0.04005379      d


# Non-glandular trichomes
anova_trichome <- 
  df_tidy %>%
  group_by(genotype, plant, leafdisc, type) %>% 
  dplyr::summarise(sum_density = sum(density_mm2)) %>% 
  mutate(log2_density = log2(sum_density + 1)) %>% 
 # filter(type == "non_glandular") %>% 
  with(., aov(log2_density ~ genotype))


hsd <- 
  HSD.test(y = anova_trichome,
           trt = "genotype", 
           group = TRUE, 
           console = FALSE)
hsd$groups

# NON-GLANDULAR
# Elite       5.1847361      a
# F1          3.8867179      a
# F1_hab      1.6624242      b
# LA1777      0.3941611      c
# PI127826    0.3394023      c

