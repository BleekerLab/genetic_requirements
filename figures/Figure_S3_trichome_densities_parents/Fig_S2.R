library(tidyverse)
library(Rmisc)
library(car)
library(multcompView)


############################
# Custom theme for plotting#
############################
my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 8, colour = "black", angle = 45),
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

df_parsed = df_tidy %>% 
  dplyr::group_by(genotype, plant, surface, type) %>% 
  dplyr::summarise(density_mm2 = mean(density_mm2)) %>% 
  # summarise over the genotypes
  summarySE(measurevar = "density_mm2", 
            groupvars = c("genotype", "surface", "type"))

df_parsed$genotype = factor(df_parsed$genotype,levels = c("Elite","F1", "PI127826", "F1_hab", "LA1777"),
                            ordered = TRUE)

df_parsed$surface= factor(df_parsed$surface,levels = c("adaxial", "abaxial"),
                            ordered = TRUE) 


# Plot only type-I/IV and no-glandular trichomes (average per abaxial / adaxial side)
p.ng.typeI <- df_parsed  %>% 
  filter(!type %in% c("type_VI", "sum_glandular")) %>%
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
  facet_grid(surface~type, scales = "fixed") +
  my.theme +
  xlab(NULL) + 
  ylab("Trichome density (trichomes/mm2)")

ggsave(file = "figures/Figure_S3_trichome_densities_parents/Fig_S3_densities_NG_typeIIV.pdf", plot = p.ng.typeI, height = 7, width = 12, units = "cm")
