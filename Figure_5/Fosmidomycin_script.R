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

df <- read.csv(file = "Figure_6/20190115_ng_trichome_all.csv", header = TRUE, check.names = FALSE)
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

###############
# Statisctics #
###############

# keep only MEP and MVA derived terpenes
df.long.parsed <- df.long %>% filter(metabolite %in% c("total_MEP_terpenes", "total_MVA_terpenes")) %>% droplevels()

# STEPS to do the statistics: 
# Create a dataframe per genotype - only use the data from day 14
# Then perform an ANOVA test per metabolite, testing the effect of the treatments
# Then perform a Tukey HSD test per metabolite to see the effect of individual treatments

# PI127826
df.long.PI = df.long.parsed %>% filter(genotype == "PI127826", day == "14")
oav.PI = lapply(split(df.long.PI, 
             df.long.PI$metabolite), 
             function(d) {aov(log(level+1) ~ treatment, data=d) })
TukeyHSD(oav.PI$total_MEP_terpenes)
TukeyHSD(oav.PI$total_MVA_terpenes)

# F2-73
df.long.73 = df.long.parsed %>% filter(genotype == "73", day == "14")
oav.73 = lapply(split(df.long.73, 
                     df.long.73$metabolite), 
               function(d) { aov(log(level+1) ~ treatment, data=d) })
TukeyHSD(oav.73$total_MEP_terpenes)
TukeyHSD(oav.73$total_MVA_terpenes)

# Elite line
df.long.CV = df.long.parsed %>% filter(genotype == "CV", day == "14")
oav.CV = lapply(split(df.long.CV, 
                     df.long.CV$metabolite), 
               function(d) { aov(log(level+1) ~ treatment, data=d) })
TukeyHSD(oav.CV$total_MEP_terpenes)
TukeyHSD(oav.CV$total_MVA_terpenes)

# F2-411
df.long.411 = df.long.parsed %>% filter(genotype == "411", day == "14")
oav.411 = lapply(split(df.long.411, 
                     df.long.411$metabolite), 
               function(d) { aov(log(level+1) ~ treatment, data=d) })
TukeyHSD(oav.411$total_MEP_terpenes)
TukeyHSD(oav.411$total_MVA_terpenes)


##################
# Cavity volumes #
##################

cavities <- read.csv(file = "Figure_6/20180808_Cavity volumes_14_days_treatment.csv", header = TRUE, check.names = FALSE)

# Calculatve volume of the cavitites in picoliter
cavities$volume_pl = cavities$volume_um / 1000

sum.cavities = summarySE(cavities, 
                measurevar = "volume_pl", 
                groupvars = c("genotype", "treatment"))
p.cavities = 
sum.cavities %>% filter(., 
                 sum.cavities$treatment != "Mevastatin" &
                 sum.cavities$genotype == "PI127826") %>%
  ggplot(., aes(x=treatment, y=volume_pl, fill = "black")) +
  geom_bar(stat = "identity", fill = "black")+
  geom_errorbar(aes(ymin = volume_pl- se, ymax = volume_pl + se), width=0.1)+
  xlab("Storage cavity")+
  ylab("Type-VI gland storage-cavity volume (picoliter)")+
  my.theme

###################################
# Arrange both in 1 plot and save #
###################################
p.both = grid.arrange(p.volatiles, p.cavities, ncol = 2)
ggsave(file = "Figure_5/plots/fosmidomycin_PI127826_phenotypes.pdf", plot = p.both, width = 5, height = 3)


################################################
# Plot Terpenes seperately supplemental Figure #
################################################

p.facet.metabolites =
sum %>% filter(genotype == "PI127826",
               treatment != "mevastatin",
               day == "14",
               !metabolite %in% c("RF_candidate", "TMA", "total_MEP_terpenes", "Total_Monoterpenes", "total_MVA_terpenes", "Total_sesquiterpenes", "total_terpenes","Total_zingiberene_and_derivatives",
                                  "B_pinene", "Geraniol", "Linalool", "Nerolidol", "R-curcumene")
) %>%
  ggplot(aes(x = treatment, y = level))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(x= treatment, ymin = level - se, ymax = level+se),
                width = 0.3)+
  facet_wrap(~metabolite, scale = "free", nrow = 4)+
  my.theme
ggsave(file = "Figure_5/plots/fosmidomycin_PI127826_facetted_metabolites.pdf", plot = p.facet.metabolites, width = 6, height = 7)

