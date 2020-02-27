if (! "checkpoint" %in% installed.packages()){
  install.packages("checkpoint")
}

library("checkpoint")
#checkpoint("2020-01-01")

# nicely formatted table summaries
library(gtsummary)
library(gt)

library(tidyverse)
library(RColorBrewer)
library(multcompView)
library(Rmisc)
library(FSA)

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

df2 = df %>% 
  dplyr::group_by(genotype, zingiberene, type_VI) %>% 
  # take the sum of abaxial + adaxial surface (leaf wash)
  dplyr::summarise(., sum_type_VI = sum(type_VI)) %>% 
  filter(sum_type_VI <= 10) %>%                      # maximum class value
  mutate(sum_type_VI = as.factor(sum_type_VI))       # convert int to factor
 
# to create C2 to C10 classes of trichome counts
levels(df2$sum_type_VI) = paste("C",
                               levels(df2$sum_type_VI),
                               sep = "") 

#####################################################
# Figure 3A: boxplot type-VI density vs zingiberene #
#####################################################

# statistics
# log2(zingiberene + 1) follows a normal distribution so ANOVA can be used

# ANOVA and Tukey HSD
aov_res <- aov(formula = log2(zingiberene + 1) ~ 1 + sum_type_VI, data = df2)
hsd_res <- TukeyHSD(aov_res, which = "sum_type_VI")
groups <- as.data.frame(multcompLetters(hsd_res$sum_type_VI[,4])$Letters)
colnames(groups) = "group"
groups$sum_type_VI = row.names(groups)

groups_trichome_class = pivot_longer(
  data = groups,
  cols = "sum_type_VI", 
  names_to = "type", 
  values_to = "class")

# plot
max_y_figure_3a = max(log2(df2$zingiberene + 1)) + 
  0.1 * max(log2(df2$zingiberene + 1))  

p.box = 
  df2 %>%
#  filter(zingiberene > 0) %>% 
  ggplot(aes(x = sum_type_VI, y = log2(zingiberene + 1))) +
  geom_boxplot(fill = "grey", 
               outlier.size = 0.5) +
  geom_jitter(size = 0.5, 
              width = 0.1, 
              height = 0.1) +
  ylab("7-epizingiberene (ion counts / leaflet)") +
  xlab("Type-VI trichome density class") +
  my_theme +
  geom_text(data = groups_trichome_class, 
            mapping = aes(x = class, y = max_y_figure_3a, label = group))


p.box 


ggsave(file = "Figure_3/Figure3A.pdf", 
       plot = p.box, 
       width = 7, 
       height = 5)

ggsave(file = "Figure_3/Figure3A.png", 
       plot = p.box, 
       width = 7, 
       height = 5)



##################################################
# Figure 3B: trichome densities versus zingiberene
##################################################

df = read.csv("Figure_3/Leafwash vs Trichome density.csv",
              header = TRUE, 
              stringsAsFactors = T, 
              check.names = FALSE)


# Change to long ('tidy') format. First the volatiles (df.long) and than also the trichomes (df.long2)
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
df.parsed = df.long2 %>% filter(., df.long2$metabolite == "total_volatiles" & df.long2$trichome_position == "Type_VI_mm") %>% filter(., !sample == "250")
df.parsed$metabolite = as.factor(df.parsed$metabolite)
df.parsed$trichome_position = as.factor(df.parsed$trichome_position)

# linear model fitting the summed terpenes to the density of trichomes
model  = lm(data = df.parsed, level ~ density)
summary(model)

# Scatterplot
p.scatter = 
  ggplot(df.parsed, aes(x=density, y=level)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, alpha  = 0.2) + 
  my_theme+
  ylab("Summed terpenes (ng / mg fresh leaf)")+
  xlab("Type-VI trichome density (trichomes / mm2)")+
  geom_text(aes(label=df.parsed$sample),hjust=1.5, vjust=0.5, size = 1.5)+
  annotate(geom = "text", x = 22, y = 320, 
           label = round(summary(model)$adj.r.squared, digits= 3), 
           size = 3) +
  annotate(geom = "text", x = 19.5, y = 320, 
           label = "r^2=",
           size = 3)

# Save plot
ggsave(file = "Figure_3/trichome_density_VS_leafwash.pdf", plot = p.scatter, height = 5, width = 4)






