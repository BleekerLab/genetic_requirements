if (! "checkpoint" %in% installed.packages()){
  install.packages("checkpoint")
}

library("checkpoint")
checkpoint("2020-01-01")

library(tidyverse)
library(RColorBrewer)
library(multcompView)
library(Rmisc)
library(FSA)
#############################
# Custom theme for plotting #
#############################

my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 8, colour = "black"),
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

trichomes = read.csv(file= "Figure_2/type_VI_trichomes_full_F2.csv", 
                     header = T, stringsAsFactors = T)
volatiles = read.csv(file = "Figure_2/GCMS_Leafwash_Full_F2_EZ_clean_data.csv", 
                     header = T, stringsAsFactors = T, check.names = F)
homo = read.csv(file = "Figure_2/homozygous_ZGB_genotypes_F2.csv", 
                header = T, stringsAsFactors = T) #List with plants homozygous for zFPS/ZIS
volatiles_VI = read.csv(file = "Figure_2/type_VI_gland_terpenes_F2.csv", 
                        header = T, stringsAsFactors = T, check.names = F)

# Take only the volatile data from plants homozygous for zFPS/ZIS + P450
volatiles.homo = left_join(homo, volatiles, by = "genotype")

# Remove rows with NA values
volatiles.homo = volatiles.homo[complete.cases(volatiles.homo), ] 

# subset parents from original data
volatiles.parents = volatiles %>% filter(., group != "F2")        

# join homozygous F2 plants with parents
volatiles = full_join(volatiles.homo, volatiles.parents)          


# Merge datasets of volatiles with trichome counting (/estimating)
df = inner_join(volatiles, trichomes)
df$genotype = as.factor(df$genotype)


################################
# Sum trichome counts #
################################

df2 = df %>% 
  dplyr::group_by(genotype) %>% 
  # take the sum of abaxial + adaxial surface (leaf wash)
  dplyr::summarise(., sum_type_VI = sum(type_VI)) %>% 
  left_join(., volatiles, by = "genotype") %>% # re-add the volatiles info
  filter(sum_type_VI <= 10)                    # maximum class

##########################################
# Boxplot type-VI density vs zingiberene #
##########################################

p.box = 
ggplot(df2, aes(x = df2$sum_type_VI, y = df2$zingiberene)) +
  geom_boxplot(aes(x = as.factor(df2$sum_type_VI), y = df2$zingiberene), 
               fill = "grey", outlier.size = 0.5) +
  geom_jitter(aes(x = as.factor(df2$sum_type_VI), y = df2$zingiberene),
              size = 0.5, width = 0.2) +
  #geom_smooth(method = "lm") +
  ylim(NA, 2000000) +
  xlab(NULL) +
  ylab("7-epizingiberene (ion counts / leaflet)") +
  xlab("Type-VI trichome-density class") +
  my.theme

ggsave(file = "Figure_3/plot_type-VI_class_vs_zingiberene.pdf", 
       plot = p.box, 
       width = 4, 
       height = 4)


##############
# Statistics #
##############

# Test for normality
shapiro.test(df2$zingiberene)

# Kruskal-Wallis test
kruskal.test(data = df2, zingiberene~sum_type_VI)

# Dunn's test for pairwise comparison
dt = dunnTest(data = df2, zingiberene ~ sum_type_VI, method = "bh")

# show signficant comparisons from Dunn's test
sig.groups = dt$res[,c(1,4)] %>% filter(., P.adj < 0.05)



#########################################
# Optional: Distribution of the classes #
#########################################

# Barplot
sum = summarySE(df2, measurevar = "zingiberene", groupvars = "sum_type_VI")

# p.density.class.zingiberene = 
ggplot(sum, aes(x = sum_type_VI, y = zingiberene))+
  geom_bar(aes(x = sum_type_VI, y = sum$zingiberene), 
           stat= "identity", 
           fill = "black") +
  geom_point(data = df2, aes(x = df2$sum_type_VI, y = df2$zingiberene))+
  geom_smooth(method = lm)+
  geom_errorbar(aes(x = sum$sum_type_VI, ymin = sum$zingiberene - se, ymax = sum$zingiberene + se), width = 0.2)+
  scale_x_continuous(breaks=c(1:10))+
  ylim(NA, 400000)+
  ylab("7-epizingiberene (ion counts / leaflet)")+
  xlab("trichome-density class")+
  theme_bw()


#############################################################################################
# Supplemental ?                                                                            #
# highlighting the selected lines for individual type-VI gland investigation                #
#############################################################################################
df.parsed  = na.omit(left_join(volatiles_VI, df2, by = "genotype"))

p.box.highlight.subset = 
  ggplot(df2, aes(x = df2$sum_type_VI, y = df2$zingiberene))+
  geom_boxplot(aes(x = as.factor(df2$sum_type_VI), y = df2$zingiberene), fill = "grey", outlier.size = 0.5)+
  geom_point(aes(x = as.factor(df2$sum_type_VI), y = df2$zingiberene),size = 0.5)+
  geom_jitter(data = df.parsed, aes(x = na.omit(sum_type_VI)-1, y = zingiberene.y), color = "red", width = .05)+
    ylim(0,500000)+
  xlab(NULL)+
  ylab("7-epizingiberene (ion counts / leaflet)")+
  xlab("Type-VI trichome-density class")+
  my.theme

ggsave(file = "Figure_3/plot_type-VI_class_vs_zingiberene_highlight_subset.pdf", 
       plot = p.box.highlight.subset, 
       width = 4, 
       height = 4)
