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

# Filter DF 
df.long.parsed <- df.long %>% filter(genotype == "PI127826", 
                                     day == "14", 
                                     metabolite == "total_terpenes",
                                     treatment != "mevastatin") %>% droplevels()

# Check normality and perform T-test
shapiro.test(log(df.long.parsed$level))
t.terpenes <- t.test(log(df.long.parsed$level) ~ df.long.parsed$treatment)

##################
# Cavity volumes #
##################

cavities <- read.csv(file = "Figure_6/20180808_Cavity volumes_14_days_treatment.csv", header = TRUE, check.names = FALSE)

# Calculatve volume of the cavitites in picoliter
cavities$volume_pl = cavities$volume_um / 1000

# Summarise data and perform statistics 
sum.cavities = summarySE(cavities, 
                measurevar = "volume_pl", 
                groupvars = c("genotype", "treatment"))

cavities.parsed <- cavities %>% filter(genotype == "PI127826", treatment != "Mevastatin") %>% droplevels()
ggplot(cavities.parsed, aes(x = volume_pl))+
  geom_density()+
  facet_grid(~treatment)

t.cavities <- t.test(cavities.parsed$volume_pl ~ cavities.parsed$treatment)

########
# Plot #
########
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

###############
# Statisctics #
###############




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

