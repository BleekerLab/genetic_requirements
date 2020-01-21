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
#Filter data
df.long %>% filter(.,
                   df.long$treatment %in% c("control", "fosmidomycin") &
                   df.long$metabolite == "total_terpenes" &
                     df.long$genotype == "PI127826" &
                     df.long$day == "14") %>%


#Boxplot
ggplot(., aes(x=treatment, y=level, fill = phenotype)) +
  geom_boxplot()+
  facet_grid(genotype~metabolite, scale = "free") +
  theme_classic() +
  theme_linedraw() +
  theme(text = element_text(family = "Times New Roman", color = "black", size = 8),
        axis.text.x = element_text(size = 10, angle = 0), 
        plot.title = element_text(size = rel(2),hjust = 0.5),
        strip.background = element_rect(colour="black", fill="grey"),
        strip.text.x = element_text(size=8, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Metabolite level (ng metabolite / type-VI trichome)") 

  ggsave("Total_cytosolic_Day7.jpg", device = "jpg", scale = 1, width = 28, height = 10, units = "cm", dpi = 300)
  
  

#Barplot

sum = summarySE(df.long, 
                measurevar = "level", 
                groupvars = c("genotype", "treatment", "day", "metabolite", "phenotype"))
write_csv(sum, "inhibitor_treatments_summary.csv")
str(sum)

sum %>% filter(., 
                   sum$metabolite %in% c("total_MVA_terpenes", "total_MEP_terpenes") & 
                     sum$day == "14") %>%
  ggplot(., aes(x=treatment, y=level, fill = phenotype)) +
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = level- se, ymax = level + se), width=0.1)+
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






