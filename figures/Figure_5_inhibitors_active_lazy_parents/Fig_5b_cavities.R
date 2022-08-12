library(tidyverse)
library(Rmisc)
library(multcompView)

#############################
# Costum theme for plotting #
#############################

my.theme = theme_bw()+
  theme(text = element_text(),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, colour = "black"), 
        #strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.text.x = element_text(size=8, colour = "black")
  )

#############
# LOAD data #
#############

df = read.csv(file = "data/20180808_Cavity volumes_14_days_treatment.csv", header = TRUE, stringsAsFactors = TRUE)
df$cutting = as.factor(df$cutting)

#drop the intermediate values used to calcutate the volume of the cavity
df = df %>% select(trichome, cutting, genotype, group, treatment, volume_um) 

df$volume_pl = as.numeric(df$volume_um)/1000 # calculate volume to picolitres 
df$treatment = factor(df$treatment, levels = c("Control", "Fosmidomycin", "Mevastatin"), ordered = TRUE)

sum = summarySE(df, measurevar = "volume_pl", groupvars = c("genotype", "treatment", "group"))
sum$genotype <- factor(sum$genotype, levels = c("PI127826", "F2-73", "CV", "F2-411"), ordered = TRUE)
#############
#   PLOT    #
#############

p.cavity = 
sum %>% filter(., genotype != "F2-303") %>%
  ggplot(., aes(x = treatment, y = volume_pl, fill = group))+
  geom_bar(stat = "identity", aes(x = treatment, y = volume_pl))+
  geom_errorbar(aes(x = treatment, ymax = volume_pl + se, ymin = volume_pl - se), width = 0.2)+
  facet_grid(~genotype)+
  scale_fill_manual(values = c("black", "grey"))+
  xlab(NULL)+
  ylab("Storage-cavity volume (picolitre)")+
  my.theme

#############
# Save plot #
#############

ggsave(file = "figures/Figure_5_inhibitors_active_lazy_parents/cavity_volumes_inhibitors.pdf", plot = p.cavity, width = 4, height = 2.5)

shapiro.test(log(df$volume_um))
ggplot(df, aes(x = log(volume_um+1)))+
  geom_density(stat = "density")+
  #xlim(0,20000)+
  facet_wrap(genotype~treatment, scale = "free")
##################
#   STATISTICS   #
##################

PI.aov = aov(data = df %>% filter(genotype == "PI127826"), volume_um ~ treatment)
summary(PI.aov)
TK.PI = TukeyHSD(PI.aov)
multcompLetters(TK.PI$treatment[,4])

F2_73.aov = aov(data = df %>% filter(genotype == "F2-73"), volume_um ~ treatment)
summary(F2_73.aov)
TK.F2_73 = TukeyHSD(F2_73.aov)
multcompLetters(TK.F2_73$treatment[,4])

CV.aov = aov(data = df %>% filter(genotype == "CV"), volume_um ~ treatment)
summary(CV.aov)
TK.CV = TukeyHSD(CV.aov)
multcompLetters(TK.CV$treatment[,4])

F2_411.aov = aov(data = df %>% filter(genotype == "F2-411"), volume_um ~ treatment)
summary(F2_411.aov)
TK.F2_411 = TukeyHSD(F2_411.aov)
multcompLetters(TK.F2_411$treatment[,4])



