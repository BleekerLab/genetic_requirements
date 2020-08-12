library(ggplot2)
library(tidyverse)
library(Rmisc)
library(gridExtra)

#############################
# Custom theme for plotting #
#############################

my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.text.x = element_text(size=8, colour = "black")
  )+
  theme_bw()

#############
# Load data #
#############

df <- read.csv(file = "Figure_5/20190115_ng_trichome_all.csv", header = TRUE, check.names = FALSE)
df$day = as.factor(df$day)
df[is.na(df)] = 0

#Filter-out F2_303 because it doesn't have the P450
df = df %>% filter(., genotype != "303")

#Make data tidy
df.long = gather(
  data = df,
  key = "metabolite",
  value = level,
  -sample, -genotype, -phenotype, -treatment,-day)
head(df.long)
df.long$metabolite = as.factor(df.long$metabolite)
df.long$level = as.numeric(df.long$level)
df.long$genotype = factor(df.long$genotype, 
                          levels = c("PI127826", "73", "CV", "411"),
                        ordered = TRUE)


sum = summarySE(df.long, 
                measurevar = "level", 
                groupvars = c("genotype", "treatment", "day", "metabolite", "phenotype"))

#Filter data

p.volatiles = 
sum %>% filter(., 
               sum$metabolite == "total_terpenes" & 
                 sum$day == "14" &
                 sum$treatment != "mevastatin" &
                 sum$genotype == "PI127826") %>%
  ggplot(., aes(x=treatment, y=level, fill = "black")) +
  geom_bar(stat = "identity", fill = "black")+
  geom_errorbar(aes(ymin = level- se, ymax = level + se), width=0.1)+
  #facet_grid(genotype~metabolite, scale = "free") +
  xlab("Volatiles")+
  ylab("Total volatiles per type-VI gland (ng)")+
  my.theme


##################
# Cavity volumes #
##################

cavities <- read.csv(file = "Figure_5/20180808_Cavity volumes_14_days_treatment.csv", header = TRUE, check.names = FALSE)

sum.cavities = summarySE(cavities, 
                measurevar = "volume_um", 
                groupvars = c("genotype", "treatment"))
p.cavities = 
sum.cavities %>% filter(., 
                 sum.cavities$treatment != "Mevastatin" &
                 sum.cavities$genotype == "PI127826") %>%
  ggplot(., aes(x=treatment, y=volume_um, fill = "black")) +
  geom_bar(stat = "identity", fill = "black")+
  geom_errorbar(aes(ymin = volume_um- se, ymax = volume_um + se), width=0.1)+
  xlab("Storage cavity")+
  ylab("Type-VI gland storage-cavity volume (um)")+
  my.theme

###################################
# Arrange both in 1 plot and save #
###################################
p.both = grid.arrange(p.volatiles, p.cavities, ncol = 2)
ggsave(file = "Figure_4/plots/fosmidomycin_PI127826_phenotypes.pdf", plot = p.both, width = 5, height = 3)
