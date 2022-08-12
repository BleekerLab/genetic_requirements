###########
# The checkpoint library checkpoint allows you to install packages 
# as they existed on CRAN on a specific snapshot.
# It ensures script reproducibility
# More info: https://rdrr.io/cran/checkpoint/
###########



library(tidyverse)
library(Rmisc)
library(car)
library(multcompView)

source("Figure_1/full_ptable.R")

############################
# Custom theme for plotting#
############################
my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 8, colour = "black", angle = 90),
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
#########################
# Import and shape data #
#########################

df = read.csv(file = "data/F1_trichome_density.csv", 
              header = T, 
              stringsAsFactors = T)

df$leafdisc = as.factor(df$leafdisc)

# Make data tidy
df_tidy = df %>% 
  gather(key = "type",
         value = "density_mm2",
         -genotype, - plant, -surface, -leafdisc, - date, -person)


# Summarise different types of trichomes (average of adaxial / abaxial side)
# First over the leafdisc (Welch Two Sample t-test, p-value == 0.6539)
# leaf_disc_one = df %>% filter(leafdisc == 1)
# leaf_disc_two = df %>% filter(leafdisc == 2)
# t.test(x = leaf_disc_one$density_mm2, y = leaf_disc_two$density_mm2)

df_parsed = df_tidy %>% 
  dplyr::group_by(genotype, plant, surface, type) %>% 
  dplyr::summarise(density_mm2 = mean(density_mm2)) %>% 
  # summarise over the genotypes
  summarySE(measurevar = "density_mm2", 
            groupvars = c("genotype", "surface", "type"))
  
df_parsed$genotype = factor(df_parsed$genotype,levels = c("Elite","F1", "PI127826", "F1_hab", "LA1777"),
                               ordered = TRUE) 


# Plot only type-VI trichomes (average of abaxial / adaxial side)
p.type_VI = df_parsed  %>% 
  filter(type == "type_VI") %>%
  ggplot(., aes(x = genotype,
                y = density_mm2)) +
  geom_bar(aes(x = genotype, 
               y = density_mm2), 
           stat = "identity",
           color = "black", 
           fill = "black") +
  geom_errorbar(aes(x = genotype, 
                    ymin = density_mm2 - se, 
                    ymax = density_mm2 + se), 
                width = 0.2) +
  facet_wrap(~ surface, scales = "fixed") +
  my.theme +
  xlab(NULL) + 
  ylab(expression("Leaf trichome density, trichomes/mm^2")) 
  
ggsave(file = "Figure_1/Fig1B.pdf", plot = p.type_VI, width = 9, height = 5.5, units = "cm")


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
y_max <- max(df_tidy$density_mm2) + 0.1 * max(df_tidy$density_mm2)

p.type_VI_alternative = 
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
  geom_text(data = groups_side, aes(x = genotype, y = 15, label = group))

p.type_VI_alternative
