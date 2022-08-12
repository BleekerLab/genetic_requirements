library(tidyverse)
library(Rmisc)
#############################
# Custom theme for plotting #
#############################

my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1, colour = "black"),
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

volatiles_VI = read.csv(file = "data/type_VI_gland_terpenes_F2.csv", header = T, stringsAsFactors = T, check.names = F) %>% filter(., !genotype %in% c("LA1777", "PI127826xLA1777", "LA1777_F1", "CV_LA1777"))
trichomes = read.table("data/20200811_F2+parents_type VI_mm2.txt",header = TRUE, check.names = FALSE)

########################### 
# Zingiberene per type- VI#
###########################


# Summarise data and make a barplot
sum.volatiles = summarySE(volatiles_VI, measurevar = "total_volatiles", groupvars = c("genotype"))

p.volatiles = 
  ggplot(sum.volatiles, x = reorder(genotype, -total_volatiles), y = total_volatiles)+
  geom_bar(stat = "identity", aes(x =reorder(genotype, -total_volatiles), y = total_volatiles))+
  #geom_errorbar(aes(x = genotype, ymin = total_volatiles - se, ymax = total_volatiles + se))+
  scale_fill_manual(values = c("gray28", "lightgrey", "gray52", "black"))+
  xlab(NULL)+
  ylab("Summed terpenes (ng / type-VI trichome)")+
  my.theme + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90))

ggsave("figures/Figure_03_Type_VI_heads/Fig3B_volatiles_type_VI.pdf", plot = p.volatiles, height = 4, width = 7)

