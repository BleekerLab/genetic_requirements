library(ggplot2)
library(tidyverse)
library(Rmisc)
df <- read.csv(file = "Figure_5/20190115_relative_values_percentage_to_control.csv", header = TRUE)
df$day = as.factor(df$day)
df[is.na(df)] = 0

#Make data tidy
df.long = gather(
  data = df,
  key = "metabolite",
  value = percentage,
  -sample, -genotype, -phenotype, -treatment,-day)
head(df.long)
df.long$metabolite = as.factor(df.long$metabolite)
df.long$percentage = as.numeric(df.long$percentage)
df.long$genotype = factor(df.long$genotype, 
                          levels = c("PI127826", "73", "CV", "411"),
                          ordered = TRUE)


# Summarise data to get averages and sd
df.sum = summarySE(df.long, 
                measurevar = "percentage", 
                groupvars = c("genotype", "treatment", "day", "metabolite", "phenotype"))


#Subset the averaged data for MVA / MEP terpenes at day 14
df.subset = df.sum %>% filter(., 
                   df.sum$metabolite %in% c("total_MVA_terpenes", "total_MEP_terpenes") & 
                     df.sum$day == "14")

p.percentage =
ggplot(df.subset, aes(x = treatment, y = percentage, fill = phenotype))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymax = percentage + se, ymin= percentage - se), width = 0.2)+
  facet_grid(metabolite~genotype)+
  labs(x = "Inhibitor", y = "Terpene levels compared to control (%)")+
  theme_bw()

ggsave("Figure_5/MVA_MEP_terpenes_day14_percentages.pdf", plot = p.percentage, width = 9, height = 9)


  
