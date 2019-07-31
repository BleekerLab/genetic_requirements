library(tidyverse)
library(ggplot2)
library(Rmisc)
library(car)

df = read.csv(file = "F1_trichomes_densitites_CSV.csv")
str(df)
df$leafdisc = as.factor(df$leafdisc)

# Make data tidy
df = 
  gather(data = df,
         key = "type",
         value = "density_mm2",
         -genotype, - plant, -surface, -leafdisc, - date, -person)

df.wide = 
  spread(data = df,
         key = "surface",
         value = "density_mm2")

df.wide$ab_ad = (df.wide$abaxial + df.wide$adaxial)/2
df.wide$genotype = factor(df.wide$genotype, 
                             levels = c("Cultivar", "PI127826", "F1", "PI127826 x LA1777", "LA1777"),
                             ordered = TRUE)

# Plot the different types of trichomes (Average of adaxial / abaxial side)
df.wide %>% filter(., type %in% c("non_glandular", "type_I_IV", "type_VI")) %>%
ggplot(., aes(x = genotype,
              y = ab_ad)) + 
  geom_boxplot(fill = "black")+
  facet_grid(~type, scale = "free")+
  theme_bw() +
  labs(x = "Genotype",
       y = "Leaf-trichome desity (trichomes/mm2)")

# Plot only type-VI trichomes (Average of abaxial / adaxial side)
df.wide %>% filter(., type == "type_VI") %>%
  ggplot(., aes(x = genotype,
                y = ab_ad)) + 
  geom_boxplot(fill = "grey")+
  geom_point(color = "black")+
  theme_bw() +
  labs(x = "Genotype",
       y = "Leaf-trichome type-VI desity (trichomes/mm2)")


# ANOVA
leveneTest(data = df.typeVI,
            ab_ad ~ genotype) # if p < 0.05 data is NOT normally distributed


df.typeVI = df.wide %>% filter(., type == "type_VI") 
aov = aov(data = df.typeVI, ab_ad ~ genotype)
summary(aov)
TukeyHSD(aov)

#outlier detection for type VI
df.typeVI = df %>% filter(., type == "type_VI")
rosnerTest(df.typeVI$density_mm2, k =2)
