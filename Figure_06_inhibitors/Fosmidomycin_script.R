library(ggplot2)
library(tidyverse)
library(PMCMRplus)

df <- read.csv(file = "Figure_6/20190115_ng_trichome_all.csv", header = TRUE, check.names = FALSE)
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

my.theme =
  theme_bw() + 
  theme(text = element_text(),
        axis.text.x = element_text(size = 10, colour = "black", angle = 30, hjust = 1),
        axis.text.y = element_text(size = 10, colour = "black"),
        strip.background = element_blank(),
  
   
        panel.border = element_rect(),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.text.x = element_text(size=8, colour = "black")
  )
  

metabolite.labels = c("Plastidial terpenes", "Cytosolic terpenes")
names(metabolite.labels) = c("total_MEP_terpenes", "total_MVA_terpenes")

#Barplot
p1 = 
df.long %>% 
  dplyr::group_by(genotype, phenotype, treatment, day, metabolite) %>%
  dplyr::summarise(mean_level = mean(level), se = sd(level)/sqrt(n())) %>%
  filter(metabolite %in% c("total_MVA_terpenes", "total_MEP_terpenes") & 
         day == "14") %>%
  
  ggplot(aes(x=treatment, y=mean_level, fill = phenotype)) +
  geom_bar(stat = "identity", fill = "black")+
  geom_errorbar(aes(x=treatment, ymin = mean_level- se, ymax = mean_level + se), width=0.1)+
  facet_grid(genotype~metabolite, scale = "free",
             labeller = labeller(metabolite = metabolite.labels)) +
  labs(y = "Metabolite level (ng/gland)")+
  my.theme

ggsave(filename = "Figure_6/MVA_MEP terpenes day 14_barplot.pdf", plot = p1, width = 4, height = 5)

###### Statistics #######

data.summary = df.long %>% 
  filter(metabolite %in% c("total_MVA_terpenes", "total_MEP_terpenes") & 
           day == "14")
write.table(data.summary, file = "Figure_5/inhibitor_treatments_summary.tsv", row.names = FALSE, sep = "\t")

# keep only MEP and MVA derived terpenes
df.long.parsed <- df.long %>% filter(metabolite %in% c("total_MEP_terpenes", "total_MVA_terpenes")) %>% droplevels()

# STEPS to do the statistics: 
# Create a dataframe per genotype - only use the data from day 14
# Then perform an ANOVA test per metabolite, testing the effect of the treatments
# Then perform a Tukey HSD test per metabolite to see the effect of individual treatments

# PI127826
df.long.PI = df.long.parsed %>% filter(genotype == "PI127826", day == "14")
oav.PI = lapply(split(df.long.PI, 
                      df.long.PI$metabolite), 
                function(d) {aov(log(level+1) ~ treatment, data=d) })
TukeyHSD(oav.PI$total_MEP_terpenes)
TukeyHSD(oav.PI$total_MVA_terpenes)

# F2-73
df.long.73 = df.long.parsed %>% filter(genotype == "73", day == "14")
oav.73 = lapply(split(df.long.73, 
                      df.long.73$metabolite), 
                function(d) { aov(log(level+1) ~ treatment, data=d) })
TukeyHSD(oav.73$total_MEP_terpenes)
TukeyHSD(oav.73$total_MVA_terpenes)

# Elite line
df.long.CV = df.long.parsed %>% filter(genotype == "CV", day == "14")
oav.CV = lapply(split(df.long.CV, 
                      df.long.CV$metabolite), 
                function(d) { aov(log(level+1) ~ treatment, data=d) })
TukeyHSD(oav.CV$total_MEP_terpenes)
TukeyHSD(oav.CV$total_MVA_terpenes)

# F2-411
df.long.411 = df.long.parsed %>% filter(genotype == "411", day == "14")
oav.411 = lapply(split(df.long.411, 
                       df.long.411$metabolite), 
                 function(d) { aov(log(level+1) ~ treatment, data=d) })
TukeyHSD(oav.411$total_MEP_terpenes)
TukeyHSD(oav.411$total_MVA_terpenes)
