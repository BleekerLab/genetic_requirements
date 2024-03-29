library(tidyverse)
library(RColorBrewer)
library(multcompView)
library(Rmisc)
library(ggpubr)

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

volatiles = read.csv(file = "data/GCMS_Leafwash_Full_F2_EZ_clean_data.csv", header = T, stringsAsFactors = T, check.names = F)
homo = read.csv(file = "data/homozygous_ZGB_genotypes_F2_plus_parents.csv", header = T, stringsAsFactors = T) #List with plants homozygous for zFPS/ZIS

# Take only the volatile data from plants homozygous for zFPS/ZIS + P450
volatiles.homo = left_join(homo, volatiles, by = "genotype")
volatiles.homo = volatiles.homo[complete.cases(volatiles.homo), ] # Remove rows with NA values
volatiles.parents = volatiles %>% filter(., group != "F2") #subset parents from original data
volatiles = full_join(volatiles.homo, volatiles.parents) 



# Calculate the mean value of the parents to show in plot
# Calculate the 95% confidence interval
sum <- volatiles.parents %>% dplyr::group_by(genotype) %>% 
  dplyr::summarise(zingi_mean = mean(zingiberene), zingi_std = sd(zingiberene), n = n()) %>% 
  mutate(CI_upper = 1.96*(zingi_std/sqrt(n)))

# How many F2s
nrow(volatiles.homo)

# # How many F2 plants fall within the 95% CI of PI127826 or higher?
nrow(volatiles.homo %>% filter(zingiberene > (514142-159176)))

##########################################
# Repeat for total ZGB+derivatives       #
# used to construct Supplemental Table 1 #
##########################################

sum <- volatiles.parents %>%
  mutate(zgb_all = zingiberene + zingiberenol + `epoxy-zingiberenol`) %>%
  dplyr::group_by(genotype) %>% 
  dplyr::summarise(zingi_mean = mean(zgb_all), zingi_std = sd(zgb_all), n = n()) %>% 
  mutate(CI_upper = 1.96*(zingi_std/sqrt(n)))

# how many plants in total?
nrow(volatiles.homo)

# How many F2 plants fall within the 95% CI of PI127826 or higher?
nrow(volatiles.homo %>% 
       mutate(zgb_all = zingiberene + zingiberenol + `epoxy-zingiberenol`) %>% 
       filter(zgb_all > (558676-188767)))

########################################################
# Plot 7epiZ levels  in histrogram                     #
# 1x standard deviation of PI127826 is used as binwith #
#######################################################

p.F2.ylim = # only show the non-0 values - includes the statistics per bin 
volatiles %>% filter(., group == "F2") %>% filter(zingiberene != 0) %>%
  ggplot(., aes(x = zingiberene))+
  #stat_bin(binwidth = 100000, colour="black", fill="white") +
  geom_histogram(binwidth = 140664, colour="black", fill="lightgrey") +
  stat_bin(binwidth = 140664, aes(y=..count.., label=..count..), geom="text", vjust = -1)+
  #facet_grid(~group) +
  xlim(-100000,2000000)+
  ylim(0,45)+
  #ylim(0,30)+
  xlab("7-epizingiberene abundance (peak area)")+
  ylab("Number of F2 genotypes")+
  my.theme


  ggsave("figures/Figure_02_Full_F2_7epiZ_density/Fig2A.pdf", plot = p.F2.ylim, width = 5, height = 4)

