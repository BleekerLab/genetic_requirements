library(ggplot2)
library(tidyverse)
library(Rmisc)
library(multcompView)

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


sum = summarySE(df.long, 
                measurevar = "level", 
                groupvars = c("genotype", "treatment", "day", "metabolite", "phenotype"))

#Filter data

p.volatiles = 
sum %>% filter(., 
               sum$metabolite == "total_terpenes" & 
                 sum$day == "14" &
                 sum$treatment != "mevastatin" &
                 sum$genotype == "PI127826") %>%
  ggplot(., aes(x=treatment, y=level, fill = "black")) +
  geom_bar(stat = "identity", fill = "black")+
  geom_errorbar(aes(ymin = level- se, ymax = level + se), width=0.1)+
  #facet_grid(genotype~metabolite, scale = "free") +
  theme_bw()


##################
# Cavity volumes #
##################

cavities <- read.csv(file = "Figure_5/20180808_Cavity volumes_14_days_treatment.csv", header = TRUE, check.names = FALSE)

sum.cavities = summarySE(cavities, 
                measurevar = "volume_um", 
                groupvars = c("genotype", "treatment"))

sum.cavities %>% filter(., 
                 sum.cavities$treatment != "Mevastatin" &
                 sum.cavities$genotype == "PI127826") %>%
  ggplot(., aes(x=treatment, y=volume_um, fill = "black")) +
  geom_bar(stat = "identity", fill = "black")+
  geom_errorbar(aes(ymin = volume_um- se, ymax = volume_um + se), width=0.1)+
  #facet_grid(genotype~metabolite, scale = "free") +
  theme_bw()
