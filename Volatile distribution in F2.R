library(tidyverse)
library(RColorBrewer)
# load data
trichomes = read.csv(file= "/Users/Ruy/Documents/Zingiberene & Derivatives/Large F2 EZ/Trichomes/type_VI_trichomes_full_F2.csv", header = T, stringsAsFactors = T)
volatiles = read.csv(file = "/Users/Ruy/Documents/Zingiberene & Derivatives/Large F2 EZ/Volatiles EZ (full F2)/R analysis/GCMS_Leafwash_Full_F2_EZ.csv", header = T, stringsAsFactors = T)
homo = read.csv(file = "/Users/Ruy/Documents/Zingiberene & Derivatives/Large F2 EZ/Trichomes/homozygous_ZGB_genotypes_F2.csv", header = T, stringsAsFactors = T)
volatiles_VI = read.csv(file = "/Users/Ruy/Documents/Zingiberene & Derivatives/Large F2 EZ/Volatiles UvA/type_VI_gland_terpenes_F2.csv", header = T, stringsAsFactors = T)
volatiles.homo = left_join(homo, volatiles, by = "genotype")
volatiles.homo.trichomes = inner_join(volatiles.homo, trichomes)
volatiles.homo.trichomes$genotype = as.factor(volatiles.homo.trichomes$genotype)
df = volatiles.homo.trichomes

acylsugars = trichomes = read.csv(file= "/Users/Ruy/Documents/Zingiberene & Derivatives/Large F2 EZ/Acylsugar measurments/20170408_Corrected emission all.csv", header = T, stringsAsFactors = T)

###########
#Trichomes #
############
sum_type_VI = df %>% dplyr::group_by(genotype) %>% dplyr::summarise(., sum_type_VI = sum(type_VI))
sum_type_VI = inner_join(sum_type_VI, volatiles) 
df2 = sum_type_VI %>% filter(., !sum_type_VI %in% c("11","12", "13", "14"))

#Boxplot
df2$sum_type_VI = as.factor(df2$sum_type_VI)
df2$zingiberene = as.numeric(df2$zingiberene)
ggplot(df2, aes(x = df2$sum_type_VI, y = df2$zingiberene))+
   geom_boxplot()+
  geom_point(size = 0.75)+
  ylim(NA, 500000)

# Barplot
sum = summarySE(df2, measurevar = "zingiberene", groupvars = "sum_type_VI")

#p.density.class.zingiberene = 
ggplot(sum, aes(x = sum_type_VI, y = zingiberene))+
  geom_bar(aes(x = sum_type_VI, y = sum$zingiberene), stat= "identity", fill = "black") +
  geom_point(data = df2, aes(x = df2$sum_type_VI, y = df2$zingiberene))+
  geom_smooth(method = lm)+
  geom_errorbar(aes(x = sum$sum_type_VI, ymin = sum$zingiberene - se, ymax = sum$zingiberene + se), width = 0.2)+
  scale_x_continuous(breaks=c(1:10))+
  ylim(NA, 400000)+
  ylab("zingiberene (ion counts / leaflet)")+
  xlab("trichome-density class")+
  theme_bw()

# Class density plot
classes = approxfun(density(as.numeric(df2$sum_type_VI)))
plot(classes, xlim = c(1,11))
hist(as.numeric(df2$sum_type_VI), breaks = 10, xlim = c(1,11))

# histogram No trichomes
p.trichomes = 
df2 %>% filter(!genotype %in% c("PI127826", "CV", "F1")) %>%
  ggplot(., aes(sum_type_VI))+
  stat_bin(binwidth = 1, colour = "grey")+
  #geom_freqpoly(binwidth = 0.5, colour = "red")+
  #xlim()+
  #scale_x_continuous(breaks= 1) +
  xlab("Trichome class")+
  ylab("Frequence amongst F2 progeny")+
  ggtitle("Distribution of 7-epizingiberen in a F2 population")+
  theme_bw()


#Calculate binwidth
2 * IQR(scale(log(df2$epoxy.zingiberenol))) / length(df2$epoxy.zingiberenol)^(1/3)

#############
# Volatiles #
#############

#Barplot volatiles
volatiles %>% dplyr::group_by(., genotype, group) %>% dplyr::summarise(., mean_zingiberene = mean(zingiberene)) %>% filter(genotype != "14830-144") %>%
  ggplot(., aes(x = reorder(genotype, -mean_zingiberene),
                y = mean_zingiberene, fill= group))+
  geom_bar(stat = "identity")+
  theme_bw()+
  theme(text = element_text(),
        axis.text.x = element_text(size = 3, angle = 90, colour = "black"))


p.zingiberene.fullF2 =
ggplot(data = df2, aes(scale(log(df2$zingiberene+1))))+
  stat_bin(binwidth = .2, colour = "grey")+
  #geom_freqpoly(binwidth = 0.2, colour = "red")+
  theme_bw()+
  #xlim(-3.5, 3.5)+
  xlab("7-epizingiberene / leaflet (log-scaled ion counts)")+
  ylab("Frequence amongst F2 progeny")+
  ggtitle("Distribution of 7-epizingiberen in a F2 population")
 

zingi.lw.scaled = approxfun(density(scale(log(df2$zingiberene+1))))
integrate(zingi.lw.scaled, -4.4,0.75)

plot(zingi.lw.scaled, xlim= c(-4,4), xlab = "scaled zingiberene content / leaf", ylab = "proportion of the F2 population")
abline(v = -0.85)
abline(v = 0.75)

z.log = approxfun(density(log(df2$zingiberene+1)))
plot(z.log, xlim= c(0,20))



########################### 
# Zingiberene per type- VI#
###########################

volatiles_VI %>% select(., genotype, total_volatiles) %>% filter(genotype %in% c("cultivar", "PI127826","F1")) %>% dplyr::group_by(.,genotype) %>% dplyr::summarise(.,  mean = mean(total_volatiles), se = sd(total_volatiles)/3)

p.hist = 
volatiles_VI %>% filter(!genotype %in% c("PI127826", "CV", "F1")) %>% select(., zingiberene) %>%
ggplot(., aes(scale(log(zingiberene+1))))+
  stat_bin(binwidth = .5, colour = "grey")+
  #geom_freqpoly(binwidth = 0.5, colour = "red")+
  #xlim()+
  xlab("7-epizingiberene / type-VI trichome (log-scaled ion counts)")+
  ylab("Frequence amongst F2 progeny")+
  ggtitle("Distribution of 7-epizingiberen")+
theme_bw()

###### Log scaled distribution of zinigberene in trichome type-VI head cells ###

#make the density function of the data
#scaled
zingi.scaled = approxfun(density(scale(log(volatiles_VI$total_volatiles+1))))

#log
vol.log = approxfun(density((log(volatiles_VI$total_volatiles+1))))
plot(vol.log, xlim = c(-1,5), xlab = "scaled zingiberene content / type-VI trichome", ylab = "proportion of the F2 population")
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
# Acylsugars #
##############

plot(density((scale(log(acylsugars$Emission)))))
hist(log(acylsugars$Emission))

acylsugars %>% dplyr::group_by(., genotype, group) %>% dplyr::summarise(., emission = mean(emission)) %>%
ggplot(., aes(x = reorder(genotype, -emission), y = emission, fill = group))+
  geom_bar(stat = "identity")+
  theme_classic()+
  theme(text = element_text(),
        axis.text.x = element_text(size = 6, angle = 90, colour = "black"))
  
