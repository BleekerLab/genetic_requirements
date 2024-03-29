library(tidyverse)
library(Rmisc)
library(car)
library(multcompView)

############################
# Custom theme for plotting#
############################
my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 8, colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black"),
        # strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.text.x = element_text(size=10, colour = "black")
  )

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


df_tidy$genotype <- factor(df_tidy$genotype,levels = c("Elite","F1", "PI127826", "F1_hab", "LA1777"),
       ordered = TRUE) 
# Summarise different types of trichomes (average of adaxial / abaxial side)
# First over the leafdisc (Welch Two Sample t-test, p-value == 0.6539)
# leaf_disc_one = df %>% filter(leafdisc == 1)
# leaf_disc_two = df %>% filter(leafdisc == 2)
# t.test(x = leaf_disc_one$density_mm2, y = leaf_disc_two$density_mm2)

df_parsed = df_tidy %>% 
  dplyr::group_by(genotype, plant, leafdisc, type) %>% 
  dplyr::summarise(density_ab_ad_mm2 = sum(density_mm2)) %>% 
  # summarise over the genotypes
  summarySE(measurevar = "density_ab_ad_mm2", 
            groupvars = c("genotype", "type"))
  
df_parsed$genotype = factor(df_parsed$genotype,levels = c("Elite","F1", "PI127826", "F1_hab", "LA1777"),
                               ordered = TRUE) 


# Plot only type-VI trichomes (average of abaxial / adaxial side)
p.type_VI = df_parsed  %>% 
  filter(type == "type_VI") %>%
  ggplot(., aes(x = genotype,
                y = density_ab_ad_mm2)) +
  geom_bar(aes(x = genotype, 
               y = density_ab_ad_mm2), 
           stat = "identity",
           color = "black", 
           fill = "black") +
  geom_errorbar(aes(x = genotype, 
                    ymin = density_ab_ad_mm2- se, 
                    ymax = density_ab_ad_mm2+ se), 
                width = 0.2) +
  xlab(NULL) + 
  ylab(expression("Leaf trichome density, trichomes/mm^2"))+
  facet_grid(~type)+
  theme_bw()+
  my.theme
  
ggsave(file = "figures/Figure_01_terpenes_type6_CVxPI127826/Fig1B.pdf", plot = p.type_VI, width = 6.5, height = 6.5, units = "cm")

##########################
# Supplemental figure S2 #
##########################
p.type_NG_I <-
df_parsed  %>% 
  filter(!type %in% c("type_VI", "sum_glandular")) %>%
  ggplot(., aes(x = genotype,
                y = density_ab_ad_mm2)) +
  geom_bar(aes(x = genotype, 
               y = density_ab_ad_mm2), 
           stat = "identity",
           color = "black", 
           fill = "black") +
  geom_errorbar(aes(x = genotype, 
                    ymin = density_ab_ad_mm2- se, 
                    ymax = density_ab_ad_mm2+ se), 
                width = 0.2) +
  facet_grid(~type) +
  xlab(NULL) + 
  ylab(expression("Leaf trichome density, trichomes/mm^2")) +
  theme_bw()+
  my.theme

ggsave(file = "figures/Figure_01_terpenes_type6_CVxPI127826/Fig_S2.pdf", plot = p.type_NG_I, width = 9, height = 6.5, units = "cm")

# calculate maximum value + small margin to positionate the HSD letters
y_max <- max(df_tidy$density_mm2) + 0.1 * max(df_tidy$density_mm2)

#p.type_VI_alternative = 
  df_tidy  %>% 
  filter(type != "sum_glandular") %>%
  ggplot(., aes(x = genotype,
                y = density_mm2)) +
  geom_boxplot() +
  geom_jitter(stat = "identity", width = 0.05) +
  facet_grid(type~ surface, scale = "free") +
  my.theme +
  xlab(NULL) + 
  ylab(expression("Leaf trichome density, trichomes/mm"^2))

p.type_VI_alternative
