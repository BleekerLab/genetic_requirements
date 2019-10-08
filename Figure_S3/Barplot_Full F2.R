library(tidyverse)

#############################
# Custom theme for plotting #
#############################

my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 1, colour = "black"),
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
volatiles = read.csv(file = "Figure_2/GCMS_Leafwash_Full_F2_EZ.csv", header = T, stringsAsFactors = T, check.names = F)
homo = read.csv(file = "Figure_2/homozygous_ZGB_genotypes_F2.csv", header = T, stringsAsFactors = T)
volatiles_VI = read.csv(file = "Figure_2/type_VI_gland_terpenes_F2.csv", header = T, stringsAsFactors = T, check.names = F)

# Take only the volatile data from plants homozygous for zFPS/ZIS + P450
volatiles.homo = left_join(homo, volatiles, by = "genotype")


###################
# Stacked barplot #
###################

volatiles.long = gather(volatiles.homo,
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
  
ggsave("Figure_S3/Barplot_volatiles_Full_F2.svg", plot = p.bar, height = 4, width = 8)
