library(ggplot2)
library(tidyverse)
library(Rmisc)


#############################
# Costum theme for plotting #
#############################

my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 8, colour = "black"), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.text.x = element_text(size=8, colour = "black"),
        legend.position = c(1,1)
  )+
  theme_bw()


#############
# Load data #
##############

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

df.subset$metabolite = revalue(df.subset$metabolite, c("total_MVA_terpenes" = "MVA-derived terpenes", "total_MEP_terpenes" = "MEP-derived terpenes"))

p.percentage =
ggplot(df.subset, aes(x = genotype, y = percentage, fill = phenotype))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymax = percentage + se, ymin= percentage - se), width = 0.2)+
  facet_grid(treatment~metabolite)+
  labs(x = NULL, y = "Terpene levels compared to control (%)")+
  scale_fill_manual(values = c("black", "grey"))+
  my.theme

ggsave("Figure_5/MVA_MEP_terpenes_day14_percentages.svg", plot = p.percentage, width = 5, height = 7)


  
