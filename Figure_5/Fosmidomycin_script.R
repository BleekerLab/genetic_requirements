library(ggplot2)
library(tidyverse)
library(Rmisc)
df <- read.csv(file = "Figure_5/20190115_ng_trichome_all.csv", header = TRUE, check.names = FALSE)
df$day = as.factor(df$day)
df[is.na(df)] = 0

#Filter-out F2_303 because it doesn't have the P450
df = df %>% filter(., genotype != "303")

#Make data tidy
df.long = gather(
  data = df,
  key = "metabolite",
  value = level,
  -sample, -genotype, -phenotype, -treatment,-day)
head(df.long)
df.long$metabolite = as.factor(df.long$metabolite)
df.long$level = as.numeric(df.long$level)
df.long$genotype = factor(df.long$genotype, 
                          levels = c("PI127826", "73", "CV", "411"),
                          ordered = TRUE)


#Barplot

df.long %>% 
  dplyr::group_by(genotype, phenotype, treatment, day, metabolite) %>%
  dplyr::summarise(mean_level = mean(level), se = sd(level)/sqrt(n())) %>%
  filter(metabolite %in% c("total_MVA_terpenes", "total_MEP_terpenes") & 
         day == "14") %>%
  
  ggplot(., aes(x=treatment, y=mean_level, fill = phenotype)) +
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = mean_level- se, ymax = mean_level + se), width=0.1)+
  facet_grid(genotype~metabolite, scale = "free") +
  theme_bw()

###### Relative values #######

df <- read.csv(file = "20190115_relative_values_percentage_to_control.csv", header = TRUE)
df$day = as.factor(df$day)
df[is.na(df)] = 0

#Filter-out F2_303 because it doesn't have the P450
df = df %>% filter(., genotype != "303")

#Make data tidy
df.long = gather(
  data = df,
  key = "metabolite",
  value = level,
  -sample, -genotype, -phenotype, -treatment,-day)
head(df.long)
df.long$metabolite = as.factor(df.long$metabolite)
df.long$level = as.numeric(df.long$level)
df.long$genotype = factor(df.long$genotype, 
                          levels = c("PI127826", "73", "CV", "411"),
                          ordered = TRUE)

  sum = summarySE(df.long, 
                  measurevar = "level", 
                  groupvars = c("genotype", "treatment", "day", "metabolite", "phenotype"))

sum %>% filter(., 
               sum$metabolite %in% c("total_MVA_terpenes", "total_MEP_terpenes") & 
                 sum$day == "7") %>%
  ggplot(., aes(x=genotype, y=level, fill = phenotype)) +
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = level- se, ymax = level + se), width=0.1)+
  facet_grid(treatment~metabolite) +
  theme_bw()+
  ylab("Percent change to control")






