library(tidyverse)
##############
# custom theme
##############
my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 8, colour = "black", angle = 90),
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

#######

df <- read.delim("Figure_1/20200722_F1s_leafwash_ng_mg_tissue.txt", 
                       header = T, 
                       sep = "\t", 
                       dec = ".", 
                       check.names = F,
                       quote = "", # to avoid quotes around terpene names
                       stringsAsFactors = F)

# Make data tidy -> select for the right genotypes -> order the genotypes for plotting

df.long = df %>% pivot_longer(cols = 4:ncol(df), names_to = "metabolite", values_to = "value") %>% 
  filter(genotype %in% c("Elite", "F1", "PI127826",  "F1_hab", "LA1777")) %>%
  mutate(genotype = factor(genotype, 
                           levels = c("Elite", "F1", "PI127826",  "F1_hab", "LA1777"),
                           ordered = TRUE)) 

# Optional:
#calculate to nanogram metabolite per gram (instead of milligram) fresh leaf
# df.long$value = df.long$value * 1000

# calculate average values for each metabolite

df.avg = df.long %>% dplyr::group_by(genotype, metabolite) %>% summarise(mean_value = mean(value), se = sd(value)/sqrt(3)) 


#####################################################
# Fig 1A: plot 7-epizingiberene and summed_terpenes #
#####################################################

p.fig1a = 
df.avg %>% filter(metabolite %in% c("7epiZ", "summed_terpenes")) %>% 
  ggplot(., aes(x = genotype,
                y = mean_value)) +
  geom_bar(aes(x = genotype, 
               y = mean_value), 
           stat = "identity",
           color = "black", 
           fill = "black") +
  geom_errorbar(aes(x = genotype, 
                    ymin = mean_value - se, 
                    ymax = mean_value + se), 
                width = 0.2)+
  facet_wrap(~metabolite, scale = "free")+
  my.theme +
  xlab(NULL) + 
  ylab(expression("Metabolite abundance (ng / mg fresh leaf)"))

ggsave(file = "Figure_1/Fig1A.pdf", plot = p.fig1a, width = 9, height = 6, units = "cm")

##########################################
# Supplemental Fig S1: plot all terpenes #
##########################################
p.figS1 = 
df.avg %>%
  ggplot(., aes(x = genotype,
                y = mean_value)) +
  geom_bar(aes(x = genotype, 
               y = mean_value), 
           stat = "identity",
           color = "black", 
           fill = "black") +
  geom_errorbar(aes(x = genotype, 
                    ymin = mean_value - se, 
                    ymax = mean_value + se), 
                width = 0.2)+
  facet_wrap(~metabolite, scale = "free", ncol = 4)+
  my.theme +
  xlab(NULL) + 
  ylab(expression("Metabolite abundance (ng / mg fresh leaf)"))

ggsave(file = "Figure_1/FigS1.pdf", plot = p.figS1, width = 18, height = 25, units = "cm")
