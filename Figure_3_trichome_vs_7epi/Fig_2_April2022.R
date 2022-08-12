library(tidyverse)
library(RColorBrewer)
library(multcompView)
library(Rmisc)
library(FSA)
library(ggpubr)

setwd(dir = "~/Documents/Github_R/genetic_requirements/")
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

trichomes = read.csv(file= "Figure_2_7epi_F2/type_VI_trichomes_full_F2.csv", header = T, stringsAsFactors = T)
volatiles = read.csv(file = "Figure_2_7epi_F2/GCMS_Leafwash_Full_F2_EZ_clean_data.csv", header = T, stringsAsFactors = T, check.names = F)
homo = read.csv(file = "Figure_2_7epi_F2/homozygous_ZGB_genotypes_F2.csv", header = T, stringsAsFactors = T) #List with plants homozygous for zFPS/ZIS
volatiles_VI = read.csv(file = "Figure_2_7epi_F2/type_VI_gland_terpenes_F2.csv", header = T, stringsAsFactors = T, check.names = F)

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

#p.F2.ylim = # only show the non-0 values - includes the statistics per bin 
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


ggsave("Figure_2_7epi_F2/plots/F2_zingiberene_histogram_stdev_zingi.svg", plot = p.F2.ylim, width = 5, height = 4)
####################################
# Distribution of trichome classes #
####################################
##############
# Figure 3B #
#############

df = read.delim("Figure_3_trichome_vs_7epi/volatiles_and_trichomes.tsv", 
                header = T, 
                stringsAsFactors = F)

# Assign the mean number of trichomes that should be in each group
df.t <- 
  df %>% 
  trichomes.parents <-
  trichomes.parents %>% 
  mutate(density = ifelse(type_VI == 1, 0,
                          ifelse(type_VI == 2, 7.5,
                                 ifelse(type_VI == 3, 22.5,
                                        ifelse(type_VI == 4, 40, 50
                                        )))))

# Filter if there is more then one sample taken
more.replicates <- df.t %>% dplyr::count(sample) %>% arrange(n) %>% filter(n != 2) %>% .$sample

# Sum adaxial and abaxial surface
df.t2 <-
  df.t %>%
  filter(!sample %in% more.replicates) %>% 
  group_by(sample) %>% 
  dplyr::summarise(ab_ad_density = sum(density))

# 7-epizingiberene content
df.zingiberene <- 
  df %>%
  group_by(sample) %>% 
  dplyr::summarise(zingiberene = mean(zingiberene))

# Fuse trichome data + 7-epizingierene data
# Assign the abaxial+adaxial summed number of trichomes to a group
df.fused <- 
  #left_join(df.zingiberene, df.t2, by = "sample") %>% 
   mutate(class = ifelse(ab_ad_density %in% c(0, 7.5), "0-10 trichomes",
                        ifelse(ab_ad_density %in% c(15), "11-20 trichomes",
                               ifelse(ab_ad_density %in% c(22.5, 30), "21-30 trichomes",
                                      ifelse(ab_ad_density %in% c(40), "31-40 trichomes",
                                             ifelse(ab_ad_density %in% c(45,47.5), "41-50 trichomes",
                                                    ifelse(ab_ad_density %in% c(50,57.5), "51-60 trichomes",
                                                           ifelse(ab_ad_density %in% c(62.5), "61-70 trichomes",
                                                                  ifelse(ab_ad_density %in% c(72.5, 80), "71-80 trichomes",
                                                                         ifelse(ab_ad_density %in% c(90),"81-90 trichomes", 
                                                                                "90-100+ trichomes"
                                                                                
                                                                         ))))))))))

# Assign class levels to trichome-density groups
df.fused$class <- factor(df.fused$class, levels = unique(df.fused$class)[c(10,9,2,8,5,7,6,3,4,1)], ordered = TRUE)
# Optional: remove "trichomes" from the labels
df.fused$class <- gsub(" trichomes", "", df.fused$class)

# Plot the distribution of the trichome classes
df.fused %>%
  dplyr::count(class) %>%
  ggplot(aes(x = class, y = n))+
  geom_col(fill = "grey", color = "black") +
  ylab("Numer of F2 plants")+
  xlab("Trichomes per leafdisc")+
  geom_text(aes(label = paste(n)), position =   position_dodge(1), vjust = -0.5, size = 3)+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))+
  ylim(0,75) +
  theme_bw()+
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 8, colour = "black", angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# ggsave(filename = "trichome_density_distribution.svg", height = 6, width = 11, units = "cm")

#######################
# Parental phenotypes #
#######################

trichomes = read.csv(file= "Figure_2_7epi_F2/type_VI_trichomes_full_F2.csv", header = T, stringsAsFactors = T)

trichomes.parents <-
trichomes  %>% 
filter(group != "F2") %>%
mutate(ID = paste0(genotype, group, name, date)) %>%
mutate(density = ifelse(type_VI == 1, 0,
                        ifelse(type_VI == 2, 7.5,
                               ifelse(type_VI == 3, 22.5,
                                      ifelse(type_VI == 4, 40, 50
))))) %>% 
  group_by(ID, group) %>% 
  dplyr::summarise(ab_ad_density = sum(density))


# Oultier analysis
trichomes.parents %>% ggplot(aes(x = group, y = ab_ad_density)) + geom_boxplot()

# Calculate average
trichomes.parents %>%
  mutate(outlier = ifelse(group == "PI127826" & ab_ad_density == 30, "yes", "no")) %>%
  filter(outlier != "yes") %>%
  group_by(group) %>% 
  dplyr::summarise(m_desity = mean(ab_ad_density),
                   SE = sd(ab_ad_density)/sqrt(n()),
                   n = n())

#############
# Figure 3C #
#############

give.n <- function(x){
  return(c(y = 18, label = length(x))) 
  # y = gives the y-position in the plot
}

# Plot trichonme density classes and 7-epizingiberene
df.fused %>% 
  mutate(zingiberene = ifelse(zingiberene ==0 & class == "0-10", 1, zingiberene)) %>%
  mutate(selected = ifelse(zingiberene != 0 ,1, ifelse(class=="0-10", 1, 0))) %>%
  mutate(chemotype = ifelse(zingiberene > (514142-159176), "PI127826", "Elite")) %>%
  # filter(selected == 1) %>%
  ggplot(aes(x = class, y = log(zingiberene)))+
  geom_boxplot(fill = "grey")+
  geom_jitter(width = 0.05, aes(color = chemotype))+
  stat_summary(fun.data = give.n, geom = "text", fun = median,
               position = position_dodge(1), size = 3)+
  ylab("7-epizingiberene (log2 ion-counts")+
  xlab("Trichomes-density class")+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
  scale_color_manual(values = c("black", "red"))+
  theme_bw()+
  #  ylim(0,40)+
  scale_fill_brewer(palette = "Greys")+
  theme(axis.text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, colour = "black", angle = 45, hjust = 1),
        axis.title = element_text(size = 10, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

# ggsave(filename = "density_vs_zingiberene.svg", height = 6, width = 11, units = "cm")