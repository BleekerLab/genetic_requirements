###########
# The checkpoint library checkpoint allows you to install packages 
# as they existed on CRAN on a specific snapshot.
# It ensures script reproducibility
# More info: https://rdrr.io/cran/checkpoint/
###########

library("checkpoint")
checkpoint("2020-01-01") # all package versions are from that date

library(tidyverse)
library(Rmisc)
library(multcompView)
library(stats)
source("Figure_1/full_ptable.R")
library("gridExtra")

##############
# custom theme
##############
my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1),
        axis.text.y = element_text(size = 8, colour = "black"),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.text.x = element_text(size=8, colour = "black")
  )+
  theme_bw()

#####################################################################################
# Figure 1A: 7-epizingiberene levels (norm. peak area) in F1s, parents, etc.
#####################################################################################
terpenes <- read.delim("Figure_01_terpenes_type6_CVxPI127826/20170425_leaf_washes.tsv", 
                       header = T, 
                       sep = "\t", 
                       dec = ".", 
                       check.names = F,
                       quote = "", # to avoid quotes around terpene names
                       stringsAsFactors = F)

# make a new column and calculate total terpene ion counts
terpenes$total_terpenes = rowSums(terpenes[,3:ncol(terpenes)]) 


# change to long ('tidy') format -> filter unwanted data (e.g. the LA1777-F1 measurements)
terpenes_long <- gather(data = terpenes, 
                        key = metabolite, 
                        value = abundance, 
                       -sample,
                       -genotype,
                       # metabolite converted to factor 
                       # (preserves original col order) 
                       factor_key=TRUE) %>% 
  mutate(genotype = factor(genotype, 
                           levels = c("CV", "F1", "PI127826",  "F1_hab", "LA1777"),
                           ordered = TRUE)) %>% 
  # keep only metabolites to be plotted
  filter(metabolite %in% c("total_terpenes", "7-epizingiberene")) 


##############################
# Figure 1A: metabolite levels
##############################
# summarise data for barplot
# summarySE: gives count, mean, standard deviation, 
# standard error of the mean, and confidence interval (default 95%). 
terpenes_mean_std = summarySE(terpenes_long, 
                              measurevar = "abundance", 
                              groupvars = c("genotype", "metabolite"))

#################
## Post-hoc tests
#################

# log transform the data (epizingiberene distribution not gaussian)
terpenes[terpenes == 0] <- 1 # to avoid Inf values 
terpenes_log10 = terpenes %>% 
  select(sample, genotype) %>% 
  cbind.data.frame(
    log10(terpenes[,3:ncol(terpenes)]))

# zingiberene
aov_genotype = aov(data = terpenes_log10,
                   formula = `7-epizingiberene` ~ genotype)
HSD_res <- TukeyHSD(aov_genotype, "genotype", ordered = TRUE)
groups <- as.data.frame(multcompLetters(HSD_res$genotype[,4])$Letters)
colnames(groups) = "7-epizingiberene"
groups$genotype = row.names(groups)
groups_zingi = pivot_longer(
  data = groups,
  cols = "7-epizingiberene", 
  names_to = "metabolite", 
  values_to = "group")


# total terpenes
aov_genotype_terpenes = aov(data = terpenes_log10,
                            formula = total_terpenes ~ genotype)
HSD_res2 <- TukeyHSD(aov_genotype_terpenes, "genotype", ordered = TRUE)
groups2 <- as.data.frame(multcompLetters(HSD_res2$genotype[,4])$Letters)
colnames(groups2) = "total_terpenes"
groups2$genotype = row.names(groups2)
groups_total = pivot_longer(
  data = groups2,
  cols = "total_terpenes", 
  names_to = "metabolite", 
  values_to = "group")

# collect letters in one unique dataframe
group_letters = rbind(groups_zingi,groups_total)

#######################
# Final plot: Figure 1A
#######################

# calculate maximum value + small margin to positionate the HSD letters
y_max_figure_1a <- max(terpenes_mean_std$abundance) + max(terpenes_mean_std$sd) 

p.figure1a =
  terpenes_mean_std %>% 
  ggplot(aes(x = genotype, y = abundance)) +
  geom_bar(stat = "identity", color = "black", fill = "black", alpha = 0.5) +
  geom_point(data = terpenes_long) +
  geom_errorbar(aes(x = genotype, 
                    ymin = abundance - se, 
                    ymax = abundance + se), 
                width = .2) +
  my.theme +
  facet_grid(~ metabolite) +
  ylab("Log10 of normalised metabolite abundance (ion counts / mg fresh leaf)")+
  xlab(NULL) +
#  coord_flip() + 
  geom_text(data = group_letters, 
            aes(x = genotype, y = y_max_figure_1a, label = group))

p.figure1a

###########
# Figure 1B
###########

#########################
# Import and shape data #
#########################

df = read.csv(file = "Figure_1/F1_trichome_density.csv", 
              header = T, 
              stringsAsFactors = T)

df$leafdisc = as.factor(df$leafdisc)

# Make data tidy
df_tidy = df %>% 
  gather(key = "type",
         value = "density_mm2",
         -genotype, - plant, -surface, -leafdisc, - date, -person)

#############
# Statistic #
#############

# function for easy comparisons of types / surface
test = function(x, y, z) {
  {x.sub = x %>% filter(type == y) %>% filter(surface == z)} #subsets the dataset
  {will = pairwise.wilcox.test(x.sub$density_mm2, x.sub$genotype, p.adjust.method = "none")} #wilcox test
  {letters_sig = multcompLetters(fullPTable(will$p.value), #compare the groups
                                 compare = "<",
                                 threshold = 0.05)}
  results = list(will, letters_sig)
  
  return(results)
  
}

# Calculate the statistics in a list called 'stats'
stats = list(
  type_VI_abaxial = test(df_parsed, "type_VI", "abaxial"),
  type_VI_adaxial = test(df_parsed, "type_VI", "adaxial"),
  
  type_I_IV_abaxial = test(df_parsed, "type_I_IV", "abaxial"),
  type_I_IV_adaxial = test(df_parsed, "type_I_IV", "adaxial"),
  
  type_non_glandular_abaxial = test(df_parsed, "non_glandular", "abaxial"),
  type_non_glandular_adaxial = test(df_parsed, "non_glandular", "adaxial")
)

# extract letters
groups_ab = as.data.frame(stats$type_VI_abaxial[[2]]$Letters)
colnames(groups_ab) = "abaxial"
groups_ab$genotype = row.names(groups_ab)

groups_ad = as.data.frame(stats$type_VI_adaxial[[2]]$Letters)
colnames(groups_ad) = "adaxial"
groups_ad$genotype = row.names(groups_ad)

groups_ad = pivot_longer(
  data = groups_ad,
  cols = "adaxial", 
  names_to = "surface", 
  values_to = "group")

groups_ab = pivot_longer(
  data = groups_ab,
  cols = "abaxial", 
  names_to = "surface", 
  values_to = "group")

groups_side = rbind(groups_ad, groups_ab)

###########
# Figure 1B
###########

# calculate maximum value + small margin to positionate the HSD letters
y_max_figure_1b = 
  df_tidy %>% 
  filter(type == "type_VI") %>% 
  select(density_mm2) %>% 
  max() 
y_max_figure_1b = y_max_figure_1b + 0.1 * y_max_figure_1b # 10% margin

p.figure1b = 
  df_tidy  %>% 
  filter(type == "type_VI") %>%
  ggplot(., aes(x = genotype,
                y = density_mm2)) +
  geom_boxplot() +
  geom_jitter(stat = "identity", width = 0.05) +
  facet_grid(~ surface, scales = "fixed") +
  my.theme +
  xlab(NULL) + 
  ylab(expression("Leaf trichome density, trichomes/mm"^2)) +
  scale_y_continuous(limits = c(0,15)) + 
  geom_text(data = groups_side, 
            aes(x = genotype, y = y_max_figure_1b, label = group))

p.figure1b

###############
# Save plots #
###############
ggsave(file = "Figure_1/figure1A.pdf", plot = p.figure1a)
ggsave(file = "Figure_1/figure1B.pdf", plot = p.figure1b)

grid.arrange(p.figure1a, p.figure1b, nrow = 1) # to visualise

g <- arrangeGrob(p.figure1a, p.figure1b, nrow = 1) # to save
ggsave(filename = "Figure_1/figure1.pdf", g)
ggsave(filename = "Figure_1/figure1.png", g)

