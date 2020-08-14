library(ggplot2)
library(tidyr)
library(dplyr)
library(car)

df = read.table(file = "20160905_BC lines GCMS_area_mg_leaf.csv", header = T, sep = ",")
str(df)

df$cutting <- as.factor(df$cutting)
df$genotype= factor(df$genotype, levels= c("CV", "PI127826", "F1", "PHS", "PHS + 00205", "ZIS", "ZIS + 00205"))
#Change to long ('tidy') format
df_long = gather(
  data = df,
  key = "metabolite",
  value = abundance,
  -sample,-cutting,-BC1_progeny,-genotype)
head(df_long)
str(df_long)
df_long$metabolite = factor(df_long$metabolite, levels = c("zingiberene", "zingiberenol", "epoxy_zingiberenol"))
df_long$genotype= factor(df_long$genotype, levels= c("CV", "PI127826", "F1", "PHS", "PHS + 00205", "ZIS", "ZIS + 00205"))
#Filter if you want
#df_zingiberene_class3_4 <- filter(df_long, metabolite == "zingiberene" & class %in% c("3", "4"))

#Plotting
ggplot(data = df_long, aes(x = genotype, y = abundance)) +
  geom_point()+
  geom_boxplot(fill="grey")+

theme_classic() +
  theme_linedraw() +
  facet_wrap(.~metabolite, scales = 'free') +
  theme(text = element_text(family = "Arial", color = "black", size = 8),
        axis.text.x = element_text(size = 6, angle = 90), 
        plot.title = element_text(size = rel(2),hjust = 0.5),
        strip.background = element_rect(colour="black", fill="white"),
        strip.text.x = element_text(size=8, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Relative metabolite abundance (ion counts / mg fresh leaf", x = "Genotype") +
  ggsave("Zingiberene&derivatives_BC_population.svg", device = "svg", scale = 1, width = 28, height = 10, units = "cm", dpi = 300)


## ANOVA Volatiles
leveneTest(df$zingiberene ~ df$genotype, data = df) # If p > 0.05 the data is NOT normally distributed [levenetest is part of the car package]



##############  Whitefly assay  ################

df_wf = read.table(file = "20160726_WF assay BC1.csv", header = T, sep = ",")
str(df_wf)
df_wf$genotype= factor(df_wf$genotype, levels= c("CV", "PI127826", "F1", "PHS", "PHS + 00205", "ZIS", "ZIS + 00205"))
df_wf$cage=as.factor(df$cage)

ggplot(data = df_wf, aes(x= genotype, y = survival)) +
  geom_boxplot(fill = "grey") +
  #geom_jitter(width=0.05)+
  #stat_boxplot(geom= "errorbar", width = 0.1)+
  theme_classic() +
  theme_linedraw() +                                       
   theme(text = element_text (family = "Times New Roman", color = "black", size = 8),
        axis.text.x = element_text(size = 9),            
        plot.title = element_text(size = rel(2),hjust = 0.5),
        strip.background = element_rect(colour="black", fill="gray"),
        strip.text.x = element_text(size=8, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Whitefly survival (%)", x = "Genotype") +
  ggsave("whitefly survival_BC_population_per class.svg", device = "svg", scale = 1, width = 28, height = 10, units = "cm", dpi = 300)

##ANOVA Whitefly
leveneTest(df_wf$survival ~ df_wf$genotype, data = df_wf) # If p > 0.05 the data is NOT normally distributed [levenetest is part of the car package]
aov_genotype = aov(df_wf$survival ~ df_wf$genotype, data=df)
aov_table <- summary(aov_genotype) 

capture.output(summary(aov_genotype),file="ANOVA_BC1_whitefly.txt")

#TukeyHSD post-hoc
capture.output(TukeyHSD(aov_genotype),file="TukeyHSD_BC1_whitefly.txt")