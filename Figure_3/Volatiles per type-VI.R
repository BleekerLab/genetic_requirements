library(tidyverse)
library(RColorBrewer)
library(multcompView)
library(Rmisc)
#############################
# Custom theme for plotting #
#############################

my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.text.x = element_text(size=8, colour = "black")
  )+
  theme_bw()

#############
# Load data #
#############

volatiles_VI = read.csv(file = "Figure_2/type_VI_gland_terpenes_F2.csv", header = T, stringsAsFactors = T, check.names = F) %>% filter(., !genotype %in% c("LA1777", "PI127826xLA1777", "LA1777_F1", "CV_LA1777"))
trichomes = read.csv(file = "Figure_3/trichome_density_F2_selection.csv", header = T, stringsAsFactors = T, check.names = F)

########################### 
# Zingiberene per type- VI#
###########################

volatiles_VI %>% select(., genotype, total_volatiles) %>% filter(genotype %in% c("cultivar", "PI127826","F1")) %>% dplyr::group_by(.,genotype) %>% dplyr::summarise(.,  mean = mean(total_volatiles), se = sd(total_volatiles)/3)

#p.hist = 
volatiles_VI %>% filter(!genotype %in% c("PI127826", "CV", "F1")) %>% select(., zingiberene) %>%
ggplot(., aes(log(zingiberene+1)))+
  stat_bin(binwidth = .1, colour = "grey")+
  #geom_freqpoly(binwidth = 0.5, colour = "red")+
  #xlim()+
  xlab("7-epizingiberene / type-VI trichome (log-scaled ion counts)")+
  ylab("Frequence amongst F2 progeny")+
  ggtitle("Distribution of 7-epizingiberen")+
theme_bw()

sum.volatiles = summarySE(volatiles_VI, measurevar = "total_volatiles", groupvars = c("genotype"))

ggplot(sum.volatiles, x = reorder(genotype, -total_volatiles), y = total_volatiles)+
  geom_bar(stat = "identity", aes(x =reorder(genotype, -total_volatiles), y = total_volatiles, fill = genotype))+
  geom_errorbar(aes(x = genotype, ymin = total_volatiles - se, ymax = total_volatiles + se))+
  my.theme




###### Log scaled distribution of zinigberene in trichome type-VI head cells ###

#make the density function of the data
#scaled
zingi.scaled = approxfun(density(scale(log(volatiles_VI$total_volatiles+1))))

#log
vol.log = approxfun(density((log(volatiles_VI$total_volatiles+1))))
plot(vol.log, xlim = c(-1.5,4.5), xlab = "scaled zingiberene content / type-VI trichome", ylab = "proportion of the F2 population")
abline(v = log(17.7+1)) # mean PI127826
abline(v = log(0.51+1)) # mean CV
abline(v = log(0.9 +1)) # mean F1

#Calculate the area's of the function: integrate(funtion, xmin, xmax)
integrate(zingi.scaled, -2.6,0.91)

#plot the data

pdf("/Users/Ruy/Documents/Zingiberene & Derivatives/Large F2 EZ/F2_volatiles_type_VI_density.pdf") 
p.density = plot(zingi.scaled, xlim = c(-3,3), xlab = "scaled zingiberene content / type-VI trichome", ylab = "proportion of the F2 population")
abline(v = 0.5) # line at 66%
abline(v = 0.91)#line at 75 %
abline(v=-0.85)

ggsave("/Users/Ruy/Documents/Zingiberene & Derivatives/Large F2 EZ/F2_volatiles_type_VI_histogram.pdf", plot = p.hist, height = 6, width = 6)
ggsave("/Users/Ruy/Documents/Zingiberene & Derivatives/Large F2 EZ/F2_zingiberene_full_set_histogram.pdf", plot = p.zingiberene.fullF2, height = 6, width = 6)
ggsave("/Users/Ruy/Documents/Zingiberene & Derivatives/Large F2 EZ/F2_trichome_class_full_set_histogram.pdf", plot = p.trichomes, height = 6, width = 6)
ggsave("/Users/Ruy/Documents/Zingiberene & Derivatives/Large F2 EZ/F2_trichome_class_vs_zingiberene.pdf", plot = p.density.class.zingiberene, height = 6, width = 10)

##############
# Statistics #
##############

aov_class = aov(log(df2$zingiberene) ~ as.factor(df2$sum_type_VI))
summary(aov_class)
TK = TukeyHSD(aov_class)
significant_groups = multcompLetters(TK$`as.factor(df2$sum_type_VI)`[,4])
  
