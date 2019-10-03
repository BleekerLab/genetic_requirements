library(Rmisc)
library(tidyverse)
library(ggplot2)

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


sum = summarySE(df, measurevar = "volume_um", groupvars = c("genotype", "phenotype"))
write.table(sum, file = "cavity_volume_F2_summary_table.txt", sep = "\t", row.names = FALSE)

p.cavity = 
ggplot(sum, aes(x = genotype, y = volume_um, fill = phenotype)) +
  geom_bar(aes(x = genotype, y = volume_um), stat = "identity") +
  scale_fill_manual(values = c("black", "grey"))+
  geom_errorbar(aes(ymin = volume_um - se, ymax = volume_um + se),
                width = 0.1)+
  xlab(NULL)+
  ylab("Storage-cavity volume (um3)")+
  my.theme +
  theme(legend.position = c(.8,.8), legend.background = element_rect(fill=NULL))


##############################
# Volatles type VI trichomes #
##############################

# load data
df2 = read.csv(file = "Figure_4/20190807_Type_VI_volatiles_High_Low_F2s.csv",
                header = TRUE)
df2 = df2 %>% filter(., genotype != "F1_hab")
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
# SAVE PLOTS #
##############

ggsave(file = "Figure_4/plots/stoarage_cavity_volume.svg", plot = p.cavity, height = 3, width = 4)
ggsave(file = "Figure_4/plots/terpenes_per_type_VI.svg", plot = p.terpenes, height = 3, width = 4)




