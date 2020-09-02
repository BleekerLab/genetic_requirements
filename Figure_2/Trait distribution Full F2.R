library(tidyverse)
library(RColorBrewer)
library(multcompView)
library(Rmisc)
library(FSA)
#############################
# Custom theme for plotting #
#############################

my.theme = 
  theme_bw() + 
  theme(text = element_text(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.text.x = element_text(size=10, colour = "black")
  )

#############
# Load data #
#############

trichomes = read.csv(file= "Figure_2/type_VI_trichomes_full_F2.csv", header = T, stringsAsFactors = T)
volatiles = read.csv(file = "Figure_2/GCMS_Leafwash_Full_F2_EZ_clean_data.csv", header = T, stringsAsFactors = T, check.names = F)
homo = read.csv(file = "Figure_2/homozygous_ZGB_genotypes_F2.csv", header = T, stringsAsFactors = T) #List with plants homozygous for zFPS/ZIS
volatiles_VI = read.csv(file = "Figure_2/type_VI_gland_terpenes_F2.csv", header = T, stringsAsFactors = T, check.names = F)

# Take only the volatile data from plants homozygous for zFPS/ZIS + P450
volatiles.homo = left_join(homo, volatiles, by = "genotype")
volatiles.homo = volatiles.homo[complete.cases(volatiles.homo), ] # Remove rows with NA values
volatiles.parents = volatiles %>% filter(., group != "F2") #subset parents from original data
volatiles = full_join(volatiles.homo, volatiles.parents) # join homozygous F2 plants with parents

# Merge datasets of volatiles with trichome counting (/estimating)
df = inner_join(volatiles, trichomes)
df$genotype = as.factor(df$genotype)

########################################
# Plot distribution of ZGB over the F2 #
########################################

# Plot different genotype groups seperatly

#p.F2 = 
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

# Create an "zoom-in" figure for low frequent F2 bins

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

#############
# Main plot #
#############

p.F2.ylim = # only show the non-0 values - includes the statistics per bin 
volatiles %>% filter(., group == "F2") %>% filter(zingiberene != 0) %>%
  ggplot(., aes(x = zingiberene))+
  #stat_bin(binwidth = 100000, colour="black", fill="white") +
  geom_histogram(binwidth = 100000, colour="black", fill="lightgrey") +
  stat_bin(binwidth = 100000, aes(y=..count.., label=..count..), geom="text", vjust = -1)+
  #facet_grid(~group) +
  xlim(-100000,2000000)+
  ylim(0,45)+
  #ylim(0,30)+
  xlab("7-epizingiberene abundance (peak area)")+
  ylab("Number of F2 genotypes")+
  my.theme
ggsave(filename = "Figure_2/plots/F2_zingiberene_histogram_no_zeroes.pdf", plot = p.F2.ylim, width = 5, height = 4)

# Calculate the mean value of the parents to show in plot
sum <- volatiles.parents %>% dplyr::group_by(genotype) %>% 
  dplyr::summarise(zingi_mean = mean(zingiberene), zingi_std = sd(zingiberene))

##############################
# Zingiberene vs derivatives #
##############################

volatiles %>% filter(., group == "F2")  %>% filter(zingiberene != 0) %>% filter(zingiberene < 1000000) %>%
  ggplot(., aes(x = zingiberene, y = `epoxy-zingiberenol`))+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ x)+
  stat_cor(
    method = "pearson")



########
p.all = grid.arrange(p.cultivar, p.F1,p.PI127826,p.F2, ncol = 1)

ggsave("Figure_2/ZGB_distribution_Full_F2.pdf", plot = p.all, height = 6, width = 3.5)
ggsave("Figure_2/F2-insert.pdf", plot = p.F2.insert, height = 1, width = 2)



#######################################
# Alternativly -> Barplot of volaties #
#######################################

#Barplot volatiles
volatiles %>% dplyr::group_by(., genotype, group) %>% dplyr::summarise(., mean_zingiberene = mean(zingiberene)) %>% filter(genotype != "14830-144") %>%
  ggplot(., aes(x = reorder(genotype, -mean_zingiberene),
                y = mean_zingiberene, fill= group))+
  geom_bar(stat = "identity")+
  theme_bw()+
  theme(text = element_text(),
        axis.text.x = element_text(size = 3, angle = 90, colour = "black"))





