if (! "checkpoint" %in% installed.packages()){
  install.packages("checkpoint")
}

library("checkpoint")
checkpoint("2020-01-01")

library(tidyverse)
library(multcompView)
library(ggrepel)
library(ggpubr)
library(gridExtra)
library(svglite)

#############################
# Custom theme for plotting #
#############################
source("theme.R")

#############
# Load data #
#############

df = read.delim("Figure_3/volatiles_and_trichomes.tsv", 
                header = T, 
                stringsAsFactors = F)

################################
# Sum trichome counts #
################################


df = df %>% 
  filter(zingiberene > 0) %>%  # keep only lines with non-null 7-epizingiberene 
  dplyr::group_by(sample, zingiberene) %>% 
  # take the sum of abaxial + adaxial surface (leaf wash)
  dplyr::summarise(., sum_type_VI = sum(type_VI)) %>% 
  filter(sum_type_VI <= 10)                    # maximum class value

 

df$sum_type_VI = as.factor(df$sum_type_VI)  # convert int to factor

# to create C2 to C10 classes of trichome counts
levels(df$sum_type_VI) = paste("C",
                               levels(df$sum_type_VI),
                               sep = "") 

# number of measurements per trichome density class
n_per_class = df %>% group_by(sum_type_VI) %>% tally() 
colnames(n_per_class) = c("class", "n_measured")
n_per_class$class = as.factor(n_per_class$class)

#####################################################
# Figure 3A: boxplot type-VI density vs zingiberene #
#####################################################

############
# statistics
############

# log2(zingiberene) follows a normal distribution (shapiro  p-value = 0.3508)
# therefore ANOVA can be used

# ANOVA and Tukey HSD
# no intercept (if trichome class = 0, then 7-epizingiberene cannot be measured)
aov_res <- aov(formula = log2(zingiberene) ~ 0 + sum_type_VI, data = df)
hsd_res <- TukeyHSD(aov_res, which = "sum_type_VI")
groups <- as.data.frame(multcompLetters(hsd_res$sum_type_VI[,4])$Letters)
colnames(groups) = "group"
groups$sum_type_VI = row.names(groups)

groups_trichome_class = 
  groups %>% 
  pivot_longer(
    data = groups,
    cols = "sum_type_VI", 
    names_to = "type", 
    values_to = "class")  

######
# plot
######
max_y_figure_3a = max(log2(df$zingiberene)) + 
  0.1 * max(log2(df$zingiberene + 1))  

labels_for_legend = c(
  "C4: 2 to 30 trichomes",
  "C5: 16 to 30 trichomes", 
  "C6: 30 to 65 trichomes",
  "C7: 45 to 90 trichomes",
  "C8: 60 to 105 trichomes",
  "C9: 80 to 125 trichomes",
  "C10:100 to 150 trichomes"
)

p_fig3a = 
  ggplot(df, aes(x = sum_type_VI, 
             y = log2(zingiberene))) +
  geom_boxplot(aes(fill = sum_type_VI), 
               outlier.size = 0.5) +
  geom_point(size = 0.5, 
              width = 0.1, 
              height = 0.1) +
  ylab("7-epizingiberene (Log2 ion counts / leaflet)") +
  xlab("Type-VI trichome-density class") +
  my_theme +
  scale_fill_brewer(name = "Class of trichome density",
                    labels = labels_for_legend,
                    palette = "Blues") 

# add post-hoc groups
p_fig3a = p_fig3a + 
  geom_text(data = groups_trichome_class,
            size = 3,
            aes(x = class, 
                y = 25, 
                label = group)) 

# add number of measurements
p_fig3a = p_fig3a + 
  geom_text(data = n_per_class,
            size = 3,
            aes(x = class, 
                y = 24, 
                label = paste("N = ", n_measured)))

p_fig3a

##################################################
# Figure 3B: trichome densities versus zingiberene
##################################################

###########
# Load data
###########
df = read.csv("Figure_3/Leafwash vs Trichome density.csv",
              header = TRUE, 
              stringsAsFactors = T, 
              check.names = FALSE)


# Change to long ('tidy') format. 
# First the volatiles (df.long) 
# and then also the trichomes (df.long2)
df.long = gather(
  data = df,
  key = "metabolite",
  value = level,
  -sample, -group, -Type_VI_mm2_Abaxial,-Type_VI_mm2_Adaxial,-Type_VI_mm)

df.long2 = gather(
  data = df.long,
  key = "trichome_position",
  value = density,
  -sample,-metabolite,-group,-level)


# Filter datafile on total_volatiles and type_VI_trichomes both abaxial and adaxial
df.parsed = df.long2 %>% 
  filter(metabolite == "total_volatiles" & trichome_position == "Type_VI_mm") %>% 
  filter(., !sample == "250")

df.parsed$metabolite = as.factor(df.parsed$metabolite)
df.parsed$trichome_position = as.factor(df.parsed$trichome_position)

############
# Statistics
############

# Scatterplot trichome density versus total terpene levels
p_fig3b  = 
  ggplot(df.parsed, aes(x = density, y = level)) +
  my_theme +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, alpha  = 0.2) + 
  ylab("Summed terpenes (ng / mg fresh leaf)") +
  xlab("Type-VI trichome density (trichomes / mm2)") +
  geom_text_repel(aes(label = df.parsed$sample), size = 3) +
  stat_regline_equation(
    label.x.npc = 0,
    label.y.npc = 1) + 
  stat_cor(
    method = "pearson",
    label.x.npc = 0, 
    label.y.npc = 0.9)
  
p_fig3b

##########################
# Figure 3: Figure 3A + 3B
##########################

# to visualise
grid.arrange(p_fig3a, p_fig3b, nrow = 1)

# to save
g <- arrangeGrob(p_fig3a, p_fig3b, nrow = 1) 
ggsave(filename = "Figure_3/Figure3.pdf", g, width = 15, height = 7)
ggsave(filename = "Figure_3/Figure3.png", g, width = 15, height = 7)
ggsave(filename = "Figure_3/Figure3.svg", g, width = 15, height = 7)


# Uncomment if you want to save individual figures
# ggsave(file = "Figure_3/Figure3A.pdf", 
#        plot = p_fig3a, 
#        width = 7, 
#        height = 5)
# 
# ggsave(file = "Figure_3/Figure3A.png", 
#        plot = p_fig3a, 
#        width = 7, 
#        height = 5)

# Save plots
# ggsave(file = "Figure_3/Figure3B.png", 
#        plot = p_fig3b, 
#        height = 5, 
#        width = 7)
# 
# 
# ggsave(file = "Figure_3/Figure3B.pdf", 
#        plot = p_fig3b, 
#        height = 5, 
#        width = 7)
=======
# Take the sum of abaxial+adaxial surface
type_VI_ab_ad = df %>% dplyr::group_by(genotype) %>% dplyr::summarise(., sum_type_VI = sum(type_VI))

# Join trichome data with volatile data
type_VI_ab_ad = inner_join(type_VI_ab_ad, volatiles) 

# remove mistakes in trichome counting (i.e. class can not exceed 10)
df2 = type_VI_ab_ad %>% filter(., !sum_type_VI %in% c("11","12", "13", "14"))

##########################################
# Boxplot type-VI density vs zingiberene #
##########################################

p.box = 
ggplot(df2, aes(x = df2$sum_type_VI, y = df2$zingiberene))+
  geom_boxplot(aes(x = as.factor(df2$sum_type_VI), y = df2$zingiberene), fill = "grey", outlier.size = 0.5)+
  geom_jitter(aes(x = as.factor(df2$sum_type_VI), y = df2$zingiberene),size = 0.5, width = 0.2)+
  #geom_smooth(method = "lm")+
  #ylim(NA, 2000000)+
  scale_y_continuous(trans='log10')+
  xlab(NULL)+
  ylab("7-epizingiberene (log10 ion counts / leaflet)")+
  xlab("Type-VI trichome-density class")+
  my.theme

ggsave("Figure_3/type-VI_class_vs_zingiberene_log10_transformed.pdf", plot = p.box, width = 4, height = 4)


##############
# Statistics #
##############

#Test for normality
shapiro.test(df2$zingiberene)

#Kruskal-Wallis test
kruskal.test(data = df2, zingiberene~sum_type_VI)

#Dunn's test for parewise comparison
dt = dunnTest(data = df2, zingiberene ~ sum_type_VI, method = "bh")
#show signficant comparisons from Dunn's test
sig.groups = dt$res[,c(1,4)] %>% filter(., P.adj < 0.05)



#########################################
# Optional: Distribution of the classes #
#########################################

# Barplot
sum = summarySE(df2, measurevar = "zingiberene", groupvars = "sum_type_VI")

#p.density.class.zingiberene = 
ggplot(sum, aes(x = sum_type_VI, y = zingiberene))+
  geom_bar(aes(x = sum_type_VI, y = sum$zingiberene), stat= "identity", fill = "black") +
  geom_point(data = df2, aes(x = df2$sum_type_VI, y = df2$zingiberene))+
  geom_smooth(method = lm)+
  geom_errorbar(aes(x = sum$sum_type_VI, ymin = sum$zingiberene - se, ymax = sum$zingiberene + se), width = 0.2)+
  scale_x_continuous(breaks=c(1:10))+
  ylim(NA, 400000)+
  ylab("zingiberene (ion counts / leaflet)")+
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


