library(tidyverse)
library(ggplot2)
library(Rmisc)
library(car)
library(multcompView)
library(rcompanion)

############################
# Custom theme for plotting#
############################

#########################
# Import and shape data #
#########################

df = read.csv(file = "Figure_1/F1_trichomes_densitites_CSV.csv", 
              header = T, 
              stringsAsFactors = T)
str(df)
df$leafdisc = as.factor(df$leafdisc)

# Make data tidy
df_tidy = 
  gather(data = df,
         key = "type",
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
            groupvars = c("genotype", "surface", "type")) %>% 
  # fix ordering of genotypes and leaf side for plots
  mutate(genotype,
         factor(genotype,
                levels = c("Cultivar","F1", "PI127826", "PI127826 x LA1777", "LA1777"),
                ordered = TRUE)) %>% 
  mutate(surface, factor(surface, 
                         levels = c("adaxial", "abaxial"),
                         ordered = TRUE))

# Plot only type-VI trichomes (Average of abaxial / adaxial side)
p.type_VI = df_parsed  %>% 
  filter(type == "type_VI") %>%
  ggplot(., aes(x = genotype,
                y = density_mm2)) +
    geom_bar(aes(x = genotype, y = density_mm2), stat = "identity", fill = "black") +
    geom_errorbar(aes(x = genotype, 
                      ymin = density_mm2 - se, 
                      ymax = density_mm2 + se), 
                  width = 0.2) +
  facet_wrap(~ surface, scales = "fixed") +
  my.theme +
  xlab(NULL) + 
  ylab(expression("Leaf trichome density, trichomes/mm"^2)) 
  

###############
# Save plot #
###############
ggsave(file = "Figure_1/figure1B.pdf", plot = p.type_VI, width = 9, height = 5.5, units = "cm")

#############
# Statistic #
#############
# Calculate the statistics in a list called 'stats'
type_VI_abaxial = test(x = df.mean.leafdisc, 
                       y = "type_VI", 
                       z = "abaxial")

type_VI_abaxial = test(test(df.mean.leafdisc, "type_VI", "abaxial"))

print(stats)

p.type_VI
