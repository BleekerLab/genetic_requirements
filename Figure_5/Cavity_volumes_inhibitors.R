library(tidyverse)
library(Rmisc)
library(multcompView)

#############################
# Costum theme for plotting #
#############################

my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, colour = "black"), 
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
# LOAD data #
#############

df = read.csv(file = "Figure_5/20180808_Cavity volumes_14_days_treatment.csv", header = TRUE, stringsAsFactors = TRUE)
df$cutting = as.factor(df$cutting)

#drop the intermediate values used to calcutate the volume of the cavity
df = df %>% select(trichome, cutting, genotype, group, treatment, volume_um) 
df$volume_um = as.numeric(df$volume_um) # set volume to numeric
df$treatment = factor(df$treatment, levels = c("Control", "Fosmidomycin", "Mevastatin"), ordered = TRUE)

sum = summarySE(df, measurevar = "volume_um", groupvars = c("genotype", "treatment", "group"))

#############
#   PLOT    #
#############

p.cavity = 
sum %>% filter(., genotype != "F2-303") %>%
  ggplot(., aes(x = treatment, y = volume_um, fill = group))+
  geom_bar(stat = "identity", aes(x = treatment, y = volume_um))+
  geom_errorbar(aes(x = treatment, ymax = volume_um + se, ymin = volume_um - se), width = 0.2)+
  facet_grid(~genotype)+
  scale_fill_manual(values = c("black", "grey"))+
  xlab(NULL)+
  ylab("Storage-cavity volume (um3)")+
  my.theme

#############
# Save plot #
#############

ggsave(file = "Figure_5/cavity_volumes_inhibitors.svg", plot = p.cavity, width = 8, height = 4)

shapiro.test(log(df$volume_um))
hist(log(df$volume_um))
##################
#   STATISTICS   #
##################

PI.aov = aov(data = df %>% filter(genotype == "PI127826"), log(volume_um) ~ treatment)
summary(PI.aov)
TK.PI = TukeyHSD(PI.aov)
multcompLetters(TK.PI$treatment[,4])

F2_73.aov = aov(data = df %>% filter(genotype == "F2-73"), log(volume_um) ~ treatment)
summary(F2_73.aov)
TK.F2_73 = TukeyHSD(F2_73.aov)
multcompLetters(TK.F2_73$treatment[,4])

CV.aov = aov(data = df %>% filter(genotype == "CV"), log(volume_um) ~ treatment)
summary(CV.aov)
TK.CV = TukeyHSD(CV.aov)
multcompLetters(TK.CV$treatment[,4])

F2_411.aov = aov(data = df %>% filter(genotype == "F2-411"), log(volume_um) ~ treatment)
summary(F2_411.aov)
TK.F2_411 = TukeyHSD(F2_411.aov)
multcompLetters(TK.F2_411$treatment[,4])



