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

setwd("~/Documents/Github_R/genetic_requirements")
#############################
# Custom theme for plotting #
#############################
my_theme <- source("theme.R")

#############
# Load data #
#############

df = read.delim("Figure_3_trichome_vs_7epi/volatiles_and_trichomes.tsv", 
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
  "C4: 2-30",
  "C5: 16-30", 
  "C6: 30-65",
  "C7: 45-90",
  "C8: 60-105",
  "C9: 80-125",
  "C10:100-150"
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
  #my_theme +
  scale_fill_brewer(name = "Class of trichome density",
                    labels = labels_for_legend,
                    palette = "Greys") 

# add post-hoc groups
p_fig3a = p_fig3a + 
  geom_text(data = groups_trichome_class,
            size = 3,
            aes(x = class, 
                y = 25, 
                label = group)) 



p_fig3a

ggsave(filename = "Figure_3/Figure3A.pdf", plot = p_fig3a, width =5, height = 3.5)


############################################################
# Calculate and plot the percentage of genotypes per class # 
# in which zingiberene is detected / undetected            #
############################################################

df = read.delim("Figure_3/volatiles_and_trichomes.tsv", 
                header = T, 
                stringsAsFactors = F)

df.zeros = df %>% 
  filter(zingiberene == 0) %>%  # keep only lines with non-null 7-epizingiberene 
  dplyr::group_by(sample, zingiberene) %>% 
  # take the sum of abaxial + adaxial surface (leaf wash)
  dplyr::summarise(., sum_type_VI = sum(type_VI)) %>% 
  filter(sum_type_VI <= 10)                    # maximum class value



df.zeros$sum_type_VI = as.factor(df.zeros$sum_type_VI)  # convert int to factor

# to create C2 to C10 classes of trichome counts
levels(df.zeros$sum_type_VI) = paste("C",
                               levels(df.zeros$sum_type_VI),
                               sep = "") 

# number of measurements per trichome density class
n_per_class = df.zeros %>% group_by(sum_type_VI) %>% tally() 
colnames(n_per_class) = c("class", "n_zeros")
n_per_class$class = as.factor(n_per_class$class)


####################################
# Now do the same for non-0 values #
####################################

df.non.zero = df %>% 
  filter(zingiberene > 0) %>%  # keep only lines with non-null 7-epizingiberene 
  dplyr::group_by(sample, zingiberene) %>% 
  # take the sum of abaxial + adaxial surface (leaf wash)
  dplyr::summarise(., sum_type_VI = sum(type_VI)) %>% 
  filter(sum_type_VI <= 10)                    # maximum class value



df.non.zero$sum_type_VI = as.factor(df.non.zero$sum_type_VI)  # convert int to factor

# to create C2 to C10 classes of trichome counts
levels(df.non.zero$sum_type_VI) = paste("C",
                                     levels(df.non.zero$sum_type_VI),
                                     sep = "")

n_non_zero = df.non.zero %>% group_by(sum_type_VI) %>% tally() %>% dplyr::rename(class = sum_type_VI) %>% dplyr::rename(n_non_zero = n)
n_per_class = left_join(n_per_class, n_non_zero, by = "class")
n_per_class[is.na(n_per_class)] <- 0

#########################
# Calculate percentages #
#########################

n_per_class$total_n <- n_per_class$n_zeros + n_per_class$n_non_zero
n_per_class$perc_zero <- (n_per_class$n_zeros / n_per_class$total_n)*100
n_per_class$perc_non_zero <- 100 - n_per_class$perc_zero

n_per_class.long <- pivot_longer(n_per_class,
                                 cols = -class,
                                 names_to = "zingiberene_measured",
                                 values_to = "number_of_genotypes")

########
# Plot #
########

n_per_class.long %>% 
  filter(zingiberene_measured %in% c("perc_non_zero", "perc_zero")) %>%
ggplot(aes(x = class, y = number_of_genotypes, group = zingiberene_measured))+
  geom_line(aes(color = zingiberene_measured))+
  geom_point(aes(color = zingiberene_measured))


#################################################
# Figure 3B: trichome densities versus zingiberene
##################################################

###########
# Load data
###########
df.density <- read.delim("Figure_3/20170620_Selected_F2_type VI_mm2.txt",
              header = TRUE, 
              stringsAsFactors = TRUE,
              sep = "\t",
              check.names = FALSE)

# Calculate average density over the replicates and adaxial / abaxixal surface
df.density.avg <- df.density %>% dplyr::group_by(genotype) %>% 
  dplyr::summarise(avg_density = mean(type_VI_density_mm2))

# load volatile data
df.terpenes <- read.delim(file = "Figure_3/20200818_leafwash_volatiles_F2.txt",
                           header = TRUE, 
                           stringsAsFactors = TRUE,
                           sep = "\t",
                           check.names = FALSE)

df.total.terpenes = df.terpenes %>% select(genotype, group, total_terpenes)

# Fuse datasets
df.density.terpenes <- inner_join(df.total.terpenes, df.density.avg, by = "genotype")

# Create labels for plotting

df.density.terpenes$labels = ""

ix_label <- c(26,27,28)
df.density.terpenes$labels[ix_label] <- df.density.terpenes$genotype[ix_label]

############
# Statistics
############
my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.y = element_text(size = 8, colour = "black"),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.text.x = element_text(size=8, colour = "black")
  )+
  theme_bw()

# Scatterplot trichome density versus total terpene levels
p_fig3b  = 
  ggplot(df.density.terpenes, aes(x = avg_density, y = total_terpenes)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, alpha  = 0.2) + 
  ylab("Total terpenes (ng / mg fresh leaf)") +
  xlab("Type-VI trichome density (trichomes / mm2)") +
  geom_text_repel(aes(label = df.density.terpenes$labels)) +
  stat_cor(
    method = "pearson",
    label.x = 7,
    label.y = 250,
    label.x.npc = 0, 
    label.y.npc = 0.9) +
    my.theme
  
p_fig3b

ggsave(filename = "Figure_3/Figure3B.pdf", plot = p_fig3b, width =3.5, height = 3.5)




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


####################################
# Distribution of trichome classes #
####################################

df.t <- 
df %>% 
  mutate(density = ifelse(type_VI == 1, 0,
                          ifelse(type_VI == 2, 7.5,
                                 ifelse(type_VI == 3, 22.5,
                                        ifelse(type_VI == 4, 40, 50
                                               )))))

more.replicates <- df.t %>% dplyr::count(sample) %>% arrange(n) %>% filter(n != 2) %>% .$sample

df.t2 <-
  df.t %>%
  filter(!sample %in% more.replicates) %>% 
  group_by(sample) %>% 
  dplyr::summarise(ab_ad_density = sum(density))


df.zingiberene <- 
  df %>%
  group_by(sample) %>% 
  dplyr::summarise(zingiberene = mean(zingiberene))

df.fused <- 
  left_join(df.t2, df.zingiberene, by = "sample") %>% 
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


df.fused$class <- factor(df.fused$class, levels = unique(df.fused$class)[c(10,9,2,8,5,7,6,3,4,1)], ordered = TRUE)
df.fused$class <- gsub(" trichomes", "", df.fused$class)

df.fused %>%
  count(class) %>%
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

ggsave(filename = "trichome_density_distribution.svg", height = 6, width = 11, units = "cm")
  
 give.n <- function(x){
  return(c(y = 18, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}
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

ggsave(filename = "density_vs_zingiberene.svg", height = 6, width = 11, units = "cm")

####
# Parental zingiberene levels

df.fused %>% 
  mutate(zingiberene = ifelse(zingiberene ==0 & class == "0-10", 1, zingiberene)) %>%
  mutate(selected = ifelse(zingiberene != 0 ,1, ifelse(class=="0-10", 1, 0))) %>%
  mutate(chemotype = ifelse(zingiberene > (514142-281327), "PI127826", "Elite")) %>%
  dplyr::count(class, chemotype) %>% 
  filter(chemotype == "PI127826")



##################################
# Selection supplemental figures #
##################################

selected_lines <- read.delim("Figure_3_trichome_vs_7epi/mean_volatiles_type_VI_glands.tsv")
selected_lines <- unique(selected_lines$genotype)

df.fused2 <- df.fused

df.fused2$sample <- gsub("14830-", "", df.fused2$sample)

df.fused2 %>% 
  mutate(zingiberene = ifelse(zingiberene ==0 & class == "0-10", 1, zingiberene)) %>%
  mutate(selected = ifelse(zingiberene != 0 ,1, ifelse(class=="0-10", 1, 0))) %>%
  mutate(chemotype = ifelse(sample %in% c(selected_lines,"073", "028"), "Selected", "Not-selected")) %>%
  mutate(label_line = ifelse(chemotype == "Selected", sample, NA)) %>%
  # filter(selected == 1) %>%
  ggplot(aes(x = class, y = log(zingiberene)))+
  geom_boxplot(fill = "grey")+
  geom_jitter(width = 0.15, aes(color = chemotype))+
  geom_text_repel(aes(label = label_line), vjust = -1)+
  ylab("7-epizingiberene (log2 ion-counts")+
  xlab("Trichomes-density class")+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
  scale_color_manual(values = c("transparent", "red"))+
  theme_bw()+
  #  ylim(0,40)+
  scale_fill_brewer(palette = "Greys")+
  theme(axis.text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, colour = "black", angle = 45, hjust = 1),
        axis.title = element_text(size = 10, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

ggsave(filename = "figure_S3A_selected_F2s_labelled.svg", height = 8, width = 13, units = "cm")
  #######
# Stats #
########

df.no.zero <- df.fused %>% 
  filter(zingiberene != 0) %>% 
  mutate(zingiberene = log(zingiberene))

#Test for normality
shapiro.test(df.no.zero$zingiberene)

anova <- aov(zingiberene ~ class, data = df.no.zero)


#Kruskal-Wallis test
kruskal.test(data = df2, zingiberene~sum_type_VI)

#Dunn's test for parewise comparison
dt = dunnTest(data = df2, zingiberene ~ sum_type_VI, method = "bh")
#show signficant comparisons from Dunn's test
sig.groups = dt$res[,c(1,4)] %>% filter(., P.adj < 0.05)


###########
# Parents #
###########

df.parents <- read.csv("Figure_2_7epi_F2/type_VI_trichomes_full_F2.csv") %>% 
  filter(group != "F2") %>% 
  mutate(density = ifelse(type_VI == 1, 0,
                          ifelse(type_VI == 2, 7.5,
                                 ifelse(type_VI == 3, 22.5,
                                        ifelse(type_VI == 4, 40, 50
                                        )))))  %>% 
  group_by(genotype, group) %>% 
  dplyr::summarise(ab_ad_density = sum(density)) %>%
  dplyr::group_by(group) %>% 
  dplyr::summarise(mean_density = mean(density))


