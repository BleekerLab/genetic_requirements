#Create_dataframe
df <- read.delim("20170425_F1s leaf washes_peak area.txt", header = T, sep = "\t", dec = ".")
str(df)

#Change to long ('tidy') format
df_long <- gather(df, metabolite, abundance, Bpinene:Abergamotoicacid, -genotype,-plant,factor_key=TRUE)

attach(df_long)

#Define variables 
df_long$plant <- factor(df_long$plant, levels = c("CV", "PI127826", "LA1777", "CVxPI127826", "CVxLA1777", "PI127826xLA1777", 
                              ordered = TRUE))
b <- factor(metabolite)
c <- abundance
df_zingiberene <- filter(df_long, 
                         metabolite == "X7epizingiberene" & df_long$plant %in% c("PI127826", "CV", "CVxPI127826","PI127826xLA1777", ordered=TRUE))

df_zingiberene$plant <- factor(df_zingiberene$plant, levels =c("PI127826", "CV", "CVxPI127826","PI127826xLA1777",
                               odered = TRUE))
#Simple boxplot
ggplot(data = df_zingiberene, aes(x = df_zingiberene$plant, y = df_zingiberene$abundance, fill=df_zingiberene$plant)) +
  geom_boxplot() + 

#For colour plot: 
  scale_fill_brewer(palette="Dark2") +
#For bw plot: 
  #scale_fill_grey()+

    geom_point() +
     theme_classic() +
      theme_linedraw() +
        #facet_wrap(b, scales = 'free') +
          theme(text = element_text(family = "Times New Roman", color = "black", size = 10),
                axis.text.x = element_text(size = 10, angle = 0), 
                plot.title = element_text(size = rel(2),hjust = 0.5),
                strip.background = element_rect(colour="black", fill="gray"),
                strip.text.x = element_text(size=8, colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()) +
  labs(title = "7-epizingiberene in extracts", 
       x = "Genotype",
       y = "7-epizingiberene (ion counts / mg fresh leaf material)") +
ggsave("F1_leafwashes_zinigberene_boxplot.svg", device = "svg", scale = 1, width = 10, height = 3, units = "cm", dpi = 300)

#### BARPLOT ###

zing_summary <- summarySE(df_zingiberene, measurevar="abundance", groupvars=c("plant", "metabolite"))

ggplot(data = zing_summary, aes(x = zing_summary$plant, y = zing_summary$abundance, fill=zing_summary$plant)) +
  geom_bar(stat = "identity", color = "black", fill = "black") +
    geom_errorbar(aes(ymin=abundance-se, ymax=abundance+se), width=.2) +
  scale_fill_brewer(palette="Dark2") +
  theme_classic() +
  theme_linedraw() +
  theme(text = element_text(family = "Times New Roman", color = "black", size = 6),
        axis.text.x = element_text(size = 6, angle = 0),
        axis.text.y = element_text(size = 6, angle = 0),
        plot.title = element_text(size = rel(2),hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(title = "7-epizingiberene in extracts", 
       x = "Genotype",
       y = "7-epizingiberene (ion counts / mg fresh leaf material)") +
  ggsave("F1_leafwashes_zinigberene_barplot.svg", device = "svg", scale = 1, width = 10, height = 6, units = "cm", dpi = 300)
  

#ANOVA
library(car)
leveneTest(df_zingiberene$abundance ~ df_zingiberene$plant, data = df_zingiberene) # If p > 0.05 the data is NOT normally distributed [levenetest is part of the car package]
aov_genotype = aov(df_zingiberene$abundance ~ genotype, data=df_zingiberene)
summary(aov_genotype)
TK <- TukeyHSD(aov_genotype, "genotype", ordered = TRUE)
TukeyHSD_results = as.data.frame(TK$genotype)
write.csv(TukeyHSD_results, file = "F1_leafwash_TukeyHSD_results.csv")
