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
  dplyr::group_by(sample, zingiberene, type_VI) %>% 
  # take the sum of abaxial + adaxial surface (leaf wash)
  dplyr::summarise(., sum_type_VI = sum(type_VI)) %>% 
  filter(sum_type_VI <= 10)                    # maximum class value
 

df2$sum_type_VI = as.factor(df2$sum_type_VI)  # convert int to factor

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

p_fig3a = 
  df2 %>%
  ggplot(aes(x = sum_type_VI, y = log2(zingiberene + 1))) +
  geom_boxplot(fill = "grey", 
               outlier.size = 0.5) +
  geom_jitter(size = 0.5, 
              width = 0.1, 
              height = 0.1) +
  ylab("7-epizingiberene (ion counts / leaflet), log2 scale") +
  xlab("Type-VI trichome density class") +
  my_theme +
  geom_text(data = groups_trichome_class, 
            mapping = aes(x = class, y = max_y_figure_3a, label = group))


p_fig3a

ggsave(file = "Figure_3/Figure3A.pdf", 
       plot = p_fig3a, 
       width = 7, 
       height = 5)

ggsave(file = "Figure_3/Figure3A.png", 
       plot = p_fig3a, 
       width = 7, 
       height = 5)



##################################################
# Figure 3B: trichome densities versus zingiberene
##################################################

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

# linear model fitting the summed terpenes to the density of trichomes
model  = lm(data = df.parsed, level ~ density)
summary(model)

# Scatterplot
p_fig3b  = 
  ggplot(df.parsed, aes(x = density, y = level)) +
  my_theme +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, alpha  = 0.2) + 
  ylab("Summed terpenes (ng / mg fresh leaf)") +
  xlab("Type-VI trichome density (trichomes / mm2)") +
  geom_text_repel(aes(label = df.parsed$sample), size = 3) +
  stat_regline_equation(
    label.x = 3, 
    label.y = 1000) + 
  stat_cor(
    method = "pearson",
    label.x = 3, 
    label.y = max(df.parsed$level))
  

p_fig3b

# Save plots
ggsave(file = "Figure_3/Figure3B.png", 
       plot = p_fig3b, 
       height = 5, 
       width = 7)


ggsave(file = "Figure_3/Figure3B.pdf", 
       plot = p_fig3b, 
       height = 5, 
       width = 7)


##########################
# Figure 3: Figure 3A + 3B
##########################

# to visualise
grid.arrange(p_fig3a, p_fig3b, nrow = 1)

# to save
g <- arrangeGrob(p_fig3a, p_fig3b, nrow = 1) 
ggsave(filename = "Figure_3/Figure3.pdf", g)
ggsave(filename = "Figure_3/Figure3.png", g)


