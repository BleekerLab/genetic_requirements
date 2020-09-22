library(Rmisc)
library(tidyverse)
library(gridExtra)
library(multcomp)


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

df = read_csv(file = "Figure_4/20190807_cavity_volumes_selected_Active_Lazy_F2.csv")

df = df %>% filter(., genotype != "F1_hab") #Remove F1_hab from the dataset

df$genotype = factor(df$genotype, levels = c("CV", "PI127826","F1", "151", "411","445","28","73", "127"),
                                             ordered = TRUE)

df$genotype = revalue(df$genotype, c("151" = "F2-151", "411" = "F2-411","445" = "F2-445","28" = "F2-28","73" = "F2-73", "127" = "F2-127"))

# Calculate gland-volume in pico-liters 
df$volume_pl = df$volume_um / 1000

sum = summarySE(df, measurevar = "volume_pl", groupvars = c("genotype", "phenotype"))
write.table(sum, file = "Figure_4/cavity_volume_F2_summary_table.txt", sep = "\t", row.names = FALSE)

p.cavity = 
ggplot(sum, aes(x = genotype, y = volume_pl, fill = phenotype)) +
  geom_bar(aes(x = genotype, y = volume_pl), stat = "identity") +
  scale_fill_manual(values = c("black", "grey"))+
  geom_errorbar(aes(ymin = volume_pl - se, ymax = volume_pl + se),
                width = 0.1)+
  xlab(NULL)+
  ylab("Storage-cavity volume (picoliter)")+
  my.theme +
  theme(legend.position = c(.8,.8), legend.background = element_rect(fill=NULL))

##############
# Statistics #
##############

aov.cavity <- aov(data = df, volume_pl ~ genotype)
summary(aov.cavity)

TK.cavity <- TukeyHSD(aov.cavity)
write.table(TK.cavity$genotype, file = "Figure_4/TukeyHSD_cavities.tsv", sep = "\t")

# Get letters from multcomp
TK.cavity2 <- glht(aov.cavity, linfct=mcp(genotype="Tukey"))
cld(TK.cavity2)


##############################
# Volatles type VI trichomes #
##############################

# load data
df2 = read.csv(file = "Figure_4/20190807_Type_VI_volatiles_High_Low_F2s.csv",
                header = TRUE)
df2 = df2 %>% filter(., genotype != "F1_hab") %>% filter(!sample %in% c("151-1","151-2","151-3"))
df2$genotype = factor(df2$genotype, levels = c("CV", "PI127826","F1", "151", "411","445","28","73", "127"),
                     ordered = TRUE)

#Remane the F2-plants
df2$genotype = revalue(df2$genotype, c("151" = "F2-151", "411" = "F2-411","445" = "F2-445","28" = "F2-28","73" = "F2-73", "127" = "F2-127"))


# BARPLOT
sum2 = summarySE(df2, measurevar = "total_terpenes", groupvars = c("genotype", "phenotype"))
write.table(sum2, file = "Total_terpenes_F2_summary_table.txt", sep = "\t", row.names = FALSE)

p.terpenes = 
  ggplot(sum2, aes(x = genotype, y = total_terpenes, fill = phenotype)) +
  geom_bar(aes(x = genotype, y = total_terpenes), stat = "identity") +
  scale_fill_manual(values = c("black", "grey"))+
  geom_errorbar(aes(ymin = total_terpenes - se, ymax = total_terpenes + se),
                width = 0.1)+
  xlab(NULL)+
  ylab("Terpenes per type-VI trichome-gland (ng)")+
  my.theme +
  theme(legend.position = "none")

##############
# Statistics #
##############

aov.terpenes <- aov(data = df2, total_terpenes ~ genotype)
summary(aov.terpenes)
TK.terpenes <- TukeyHSD(aov.terpenes)
write.table(TK.terpenes$genotype, file = "Figure_4/TukeyHSD_terpenes.txt", sep = "\t")

TK.terpenes2 <- glht(aov.terpenes, linfct=mcp(genotype="Tukey"))
cld(TK.terpenes2)


##############
# SAVE PLOTS #
##############

p.both = grid.arrange(p.cavity, p.terpenes, ncol = 1)
ggsave(file = "Figure_4/plots/caviy_volume_terpenes.pdf", plot = p.both, height = 6, width = 5)

ggsave(file = "Figure_4/plots/stoarage_cavity_volume.svg", plot = p.cavity, height = 3, width = 4)
ggsave(file = "Figure_4/plots/terpenes_per_type_VI.svg", plot = p.terpenes, height = 3, width = 4)




