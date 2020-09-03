library(tidyverse)
library(ggrepel)
library(ggpubr)
library(Rmisc)

######################
# Theme for plotting #
######################

my.theme = theme_bw()+ 
  theme(text = element_text(),
        axis.text.x = element_text(size = 8, colour = "black", angle = 90, vjust = 0.6),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.y = element_text(size = 8, colour = "black"),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.text.x = element_text(size=8, colour = "black")
  )





df.terpenes <-read.csv(file = "Figure_2/type_VI_gland_terpenes_F2.csv", header = T, stringsAsFactors = T, check.names = F) %>% 
  filter(., !genotype %in% c("LA1777", "PI127826xLA1777", "LA1777_F1", "CV_LA1777"))

###########
# BARPLOT #
###########
sum.terpenes <- summarySE(df.terpenes, measurevar = "total_volatiles", groupvars = c("genotype", "group")) %>%
  select(genotype, total_volatiles, group) 

p.bar = 
ggplot(sum.terpenes, aes(x = reorder(genotype, -total_volatiles), y = total_volatiles, group = group))+
  geom_bar(stat = "identity", aes(fill = group, color = group)) +
  scale_fill_manual(values = c("darkgreen", "yellow", "grey", "red"))+
  scale_color_manual(values = c("black", "black", "black", "black"))+
  ylab("Total terpenes (ng/ type-VI gland)")+
  xlab(NULL)+
  my.theme

ggsave(filename = "Figure_3/type_VI_volatiles_barplot.pdf", plot = p.bar, height = 4, width = 5.5)


#######################################
# Gland volatiles vs trichome density #
#######################################

df.density <- read.delim("Figure_3/20170620_Selected_F2_type VI_mm2.txt",
                         header = TRUE, 
                         stringsAsFactors = TRUE,
                         sep = "\t",
                         check.names = FALSE)

# Calculate average density over the replicates and adaxial / abaxixal surface
sum.density <- df.density %>% dplyr::group_by(genotype) %>% 
  dplyr::summarise(avg_density = mean(type_VI_density_mm2))

# Join the two df's

df.terpenes.density <- inner_join(sum.terpenes, sum.density, by = "genotype")
df.terpenes.density$label = ""

df.names<- as.data.frame(c("Elite", "F1", "PI127826", "411", "151", "445", "28", "73", "127"))
colnames(df.names) <- "genotype"
df.names$label <- df.names$genotype
df.terpenes.density = left_join(df.terpenes.density, df.names, by = "genotype")


########
# Plot #
########

# Scatterplot trichome density versus total terpene levels

  ggplot(df.terpenes.density, aes(x = avg_density, y = total_volatiles)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, alpha  = 0.2) + 
  ylab("Total terpenes (ng / mg fresh leaf)") +
  xlab("Type-VI trichome density (trichomes / mm2)") +
 # geom_text_repel(aes(label = df.terpenes.density$label)) +
  stat_cor(
    method = "pearson",
    label.x = 10,
    label.y =30,
    label.x.npc = 0, 
    label.y.npc = 0.9) +
  my.theme


