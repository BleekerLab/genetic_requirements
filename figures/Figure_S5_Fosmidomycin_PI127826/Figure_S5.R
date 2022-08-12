library(tidyverse)
library(Rmisc)

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

df <- read.csv(file = "data/20190115_ng_trichome_all.csv", header = TRUE, check.names = FALSE)
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

df.long$metabolite = as.factor(df.long$metabolite)
df.long$level = as.numeric(df.long$level)
df.long$genotype = factor(df.long$genotype, 
                          levels = c("PI127826", "73", "CV", "411"),
                        ordered = TRUE)


sum = summarySE(df.long, 
                measurevar = "level", 
                groupvars = c("genotype", "treatment", "day", "metabolite", "phenotype"))


#############
# Figure S5 #
#############

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

