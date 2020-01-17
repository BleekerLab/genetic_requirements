library(tidyverse)
library(RColorBrewer)
library(multcompView)
#############################
# Custom theme for plotting #
#############################

my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 8, colour = "black"),
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

trichomes = read.csv(file= "Figure_2/type_VI_trichomes_full_F2.csv", header = T, stringsAsFactors = T)
volatiles = read.csv(file = "Figure_2/GCMS_Leafwash_Full_F2_EZ_clean_data.csv", header = T, stringsAsFactors = T, check.names = F)
homo = read.csv(file = "Figure_2/homozygous_ZGB_genotypes_F2.csv", header = T, stringsAsFactors = T)
volatiles_VI = read.csv(file = "Figure_2/type_VI_gland_terpenes_F2.csv", header = T, stringsAsFactors = T, check.names = F)

# Take only the volatile data from plants homozygous for zFPS/ZIS + P450
volatiles.homo = left_join(homo, volatiles, by = "genotype")
volatiles.homo = volatiles.homo[complete.cases(volatiles.homo), ] # Remove rows with NA values
volatiles.parents = volatiles %>% filter(., group != "F2") #subset parents from original data
volatiles = full_join(volatiles.homo, volatiles.parents) # join homozygous F2 plants with parents

# Include the trichome couting data
df = inner_join(volatiles, trichomes)
df$genotype = as.factor(df$genotype)

acylsugars = trichomes = read.csv(file= "Figure_2/20170408_Corrected emission all.csv", header = T, stringsAsFactors = T)

###########
#Trichomes #
############

# Take the sum of abaxial / adaxial sites and join volatile data
type_VI_ab_ad = df %>% dplyr::group_by(genotype) %>% dplyr::summarise(., sum_type_VI = sum(type_VI))
type_VI_ab_ad = inner_join(type_VI_ab_ad, volatiles) 
df2 = type_VI_ab_ad %>% filter(., !sum_type_VI %in% c("11","12", "13", "14")) # remove mistakes in trichome counting (i.e. a class can not exceed 10)

##########################################
# Boxplot type-VI density vs zingiberene #
##########################################
p.box = 
ggplot(df2, aes(x = df2$sum_type_VI, y = df2$zingiberene))+
  #geom_boxplot(aes(x = as.factor(df2$sum_type_VI), y = df2$zingiberene), fill = "grey", outlier.size = 0.5)+
  geom_jitter(aes(x = as.factor(df2$sum_type_VI), y = df2$zingiberene),size = 0.5)+
  #geom_smooth(method = "lm")+
  ylim(NA, 100000)+
  xlab(NULL)+
  ylab("7-epizingiberene (ion counts / leaflet")+
  xlab("Type-VI trichome-density class")+
  my.theme

ggsave(file = "Figure_2/plots/type-VI_class_vs_zingiberene.svg", plot = p.box, width = 4, height = 3)

model= polym(formula = zingiberene~sum_type_VI, data = df2, Hess = TRUE)
summary(model)

### Supplemental:  highlighting the selected lines ####
df.parsed  = na.omit(left_join(volatiles_VI, df2, by = "genotype"))

p.box.highlight.subset = 
  ggplot(df2, aes(x = df2$sum_type_VI, y = df2$zingiberene))+
  geom_boxplot(aes(x = as.factor(df2$sum_type_VI), y = df2$zingiberene), fill = "grey", outlier.size = 0.5)+
  geom_point(aes(x = as.factor(df2$sum_type_VI), y = df2$zingiberene),size = 0.5)+
  geom_jitter(data = df.parsed, aes(x = na.omit(sum_type_VI)-1, y = zingiberene.y), color = "red", width = .05)+
    ylim(0,500000)+
  xlab(NULL)+
  ylab("7-epizingiberene (ion counts / leaflet)")+
  xlab("Type-VI trichome-density class")+
  my.theme

ggsave(file = "Figure_2/plots/type-VI_class_vs_zingiberene_highlight_subset.svg", plot = p.box.highlight.subset, width = 4, height = 4)

##############
# Statistics #
##############

aov_class = aov(log(df2$zingiberene) ~ as.factor(df2$sum_type_VI))
summary(aov_class)
TK = TukeyHSD(aov_class)
significant_groups = multcompLetters(TK$`as.factor(df2$sum_type_VI)`[,4])

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
  my.theme


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
  my.theme

ggsave(file = "Figure_2/plots/zingiberene_full_F2_histogram.svg", plot = p.zingiberene.fullF2, width = 4, height = 3)

# Frequency discribution of zgb

p.F2 = 
  volatiles %>% filter(., group == "F2") %>%
  ggplot(., aes(x = zingiberene))+
  stat_bin(binwidth = 100000, colour = "black", fill = "black") +
  #facet_grid(~group) +
  xlim(-100000,1500000)+
  my.theme  + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank()) 

p.F1 =
  volatiles %>% filter(., group == "F1") %>%
  ggplot(., aes(x = zingiberene))+
  stat_bin(binwidth = 100000, colour = "black", fill = "gray53") +
  #facet_grid(~group) +
  xlim(-100000,1500000)+
  my.theme  + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank()) 

p.cultivar =
  volatiles %>% filter(., group == "cultivar") %>%
  ggplot(., aes(x = zingiberene))+
  stat_bin(binwidth = 100000, colour = "black", fill = "gray76") +
  #facet_grid(~group) +
  xlim(-100000,1500000)+
  my.theme  + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank()) 

p.PI127826 =
  volatiles %>% filter(., group == "PI127826") %>%
  ggplot(., aes(x = zingiberene))+
  stat_bin(binwidth = 100000, colour = "black", fill = "gray30") +
  #facet_grid(~group) +
  xlim(-100000,1500000)+
  my.theme  + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank()) 

p.all = grid.arrange(p.cultivar, p.F1,p.PI127826,p.F2, ncol = 1)

p.F2.insert = 
  volatiles %>% filter(., group == "F2") %>%
  ggplot(., aes(x = zingiberene))+
  stat_bin(binwidth = 100000, colour = "black", fill = "black") +
  #facet_grid(~group) +
  xlim(100000,2500000)+
  ylim(0,10)+
  theme(text = element_text(),
        axis.text.x = element_text(size = 6, colour = "black"),
        axis.text.y = element_text(size = 4, colour = "black"),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.text.x = element_text(size=8, colour = "black")
  )+
  theme_bw()+
  theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank()) 

ggsave("Figure_2/ZGB_distribution_Full_F2.pdf", plot = p.all, height = 6, width = 3.5)
ggsave("Figure_2/F2-insert.pdf", plot = p.F2.insert, height = 1, width = 2)

##############
# Statistics #
##############

PI = volatiles %>% filter(., group == "PI127826")
summary(PI)

#######################
zingi.lw.scaled = approxfun(density(scale(log(df2$zingiberene+1))))
integrate(zingi.lw.scaled, -4.4,0.75)

plot(zingi.lw.scaled, xlim= c(-4,4), xlab = "scaled zingiberene content / leaf", ylab = "proportion of the F2 population")
abline(v = -0.85)
abline(v = 0.75)

z.log = approxfun(density(log(df2$zingiberene+1)))
plot(z.log, xlim= c(0,20))

###################
# Stacked barplot #
###################

volatiles.long = gather(volatiles,
                        key = "metabolite",
                        value = "value",
                        -sample, -genotype, -group)

volatiles.long$metabolite = factor(volatiles.long$metabolite, levels = c("B_phellandrene", "carene", "caryophyllene", "elemene", "ocymene", "pinene", "zingiberene", "zingiberenol", "epoxy-zingiberenol",
                                                                         ordered = TRUE))
p.bar = 
ggplot(volatiles.long)+
  geom_bar(stat = "identity", aes(x = reorder(genotype, -value,sum), y = value, fill = metabolite))+
  scale_fill_manual(values = c("grey", "grey","grey", "grey","grey", "grey","black", "black","black"))+
  my.theme +
  ylab("Summed volatiles (ion counts / leaflet)")+
  xlab("F2 progeny")
  
ggsave("Figure_S3/Barplot_volatiles_Full_F2.svg", plot = p.bar, height = 4, width = 12)

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
  
