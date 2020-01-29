library(tidyverse)
library(ggplot2)
library(Rmisc)
library(car)
library(multcompView)
library(rcompanion)

############################
# Custom theme for plotting#
############################

my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 8, colour = "black", angle = 30, hjust = 1),
        axis.text.y = element_text(size = 8, colour = "black"),
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

df = read.csv(file = "Figure_1_F1_phenotypes/F1_trichomes_densitites_CSV.csv", header = T, stringsAsFactors = T)
str(df)
df$leafdisc = as.factor(df$leafdisc)

# Make data tidy
df = 
  gather(data = df,
         key = "type",
         value = "density_mm2",
         -genotype, - plant, -surface, -leafdisc, - date, -person)


# Summerise different types of trichomes (Average of adaxial / abaxial side)
# First over the leafdisc
df.mean.leafdisc = as.data.frame(df %>% dplyr::group_by(genotype, plant, surface, type) %>% dplyr::summarise(density_mm2 = mean(density_mm2)))

# Summarise over the genotypes
sum = summarySE(df.mean.leafdisc, measurevar = "density_mm2", groupvars = c("genotype", "surface", "type"))
sum$genotype = factor(sum$genotype, 
                          levels = c("Cultivar","F1", "PI127826",  "PI127826 x LA1777", "LA1777"),
                          ordered = TRUE)

sum$surface = factor(sum$surface, levels = c("adaxial", "abaxial"), ordered = TRUE)

# Plot only type-VI trichomes (Average of abaxial / adaxial side)
p.type_VI =
sum  %>% filter(type == "type_VI") %>%
  ggplot(., aes(x = genotype,
                y = density_mm2)) +
  geom_bar(aes(x = genotype, y = density_mm2), stat = "identity", fill = "black")+
  geom_errorbar(aes(x = genotype, ymin = density_mm2 - se, ymax = density_mm2+se), width = 0.2)+
  facet_grid(~surface, scale = "free")+
  my.theme +
  labs(x = NULL,
       y = "Leaf-trichome desity (trichomes/mm2)")

p.all =
  sum  %>% filter(!type == "type_VI") %>% filter(!type == "sum_glandular") %>% 
  ggplot(., aes(x = genotype,
                y = density_mm2)) +
  geom_bar(aes(x = genotype, y = density_mm2), stat = "identity", fill = "black")+
  geom_errorbar(aes(x = genotype, ymin = density_mm2 - se, ymax = density_mm2+se), width = 0.2)+
  facet_grid(surface~type, scale = "free")+
  my.theme +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1,  colour = "black"))+
  labs(x = NULL,
       y = "Leaf-trichome desity (trichomes/mm2)")

###############
# Print plots #
###############
ggsave(file = "Figure_1_F1_phenotypes/plots/trichome_densities_F1.svg", plot = p.all, height = 8, width = 10)
ggsave(file = "Figure_1_F1_phenotypes/plots/type_VI_densities_F1.pdf", plot = p.type_VI, width = 9, height = 5.5, units = "cm")

#############
# Statistic #
#############

#function for easy comparisons of types / surface
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
type_VI_abaxial = test(df.mean.leafdisc, "type_VI", "abaxial"),
type_VI_adaxial = test(df.mean.leafdisc, "type_VI", "adaxial"),

type_I_IV_abaxial = test(df.mean.leafdisc, "type_I_IV", "abaxial"),
type_I_IV_adaxial = test(df.mean.leafdisc, "type_I_IV", "adaxial"),

type_non_glandular_abaxial = test(df.mean.leafdisc, "non_glandular", "abaxial"),
type_non_glandular_adaxial = test(df.mean.leafdisc, "non_glandular", "adaxial")
)
print(stats)
