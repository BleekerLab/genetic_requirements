library(tidyverse)
library(ggrepel)

######################
# Theme for plotting #
######################

my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 8, colour = "black"),
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

#############################
#  calculated from leafwash #
#############################

# Selected F2
df.terpenes <- read.delim(file = "Figure_3/20200818_leafwash_volatiles_F2.txt",
                          header = TRUE, 
                          stringsAsFactors = TRUE,
                          sep = "\t",
                          check.names = FALSE)



df.terpenes %>% 
  ggplot(aes(x = log(df.terpenes$`B-caryophyllene`+1), y = log(df.terpenes$`7-epizingiberene`+1)))+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ x, alpha  = 0.2)+
  stat_cor(
    method = "pearson",
    label.x.npc = 0, 
    label.y.npc = 0.9) +
  xlab("B-caryophyllene (log-transformed ng / mg tissue)")+
  ylab("7-epizingiberene (log-transformed ng / mg tissue)")+
  my.theme

#####################################
#  Taken from Enza F2 measurements  #
#####################################

# Big F2 measured at EZ
df = read.delim("Figure_3/volatiles_and_trichomes.tsv", 
                header = T, 
                stringsAsFactors = F)
df %>%
ggplot(aes(x = log(caryophyllene) , y = log(zingiberene)))+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ x, alpha  = 0.2)+
  stat_cor(
    method = "pearson",
    label.x.npc = 0, 
    label.y.npc = 0.9)


#################################
# Data from trichome-head cells #
#################################

df <- read.csv(file = "Figure_2/type_VI_gland_terpenes_F2.csv", 
                header = T, 
                stringsAsFactors = F)

df.long <- pivot_longer(df, cols = 4:ncol(df),
                        names_to = "metabolite",
                        values_to = "value")

df.avg = df.long %>% dplyr::group_by(genotype, metabolite) %>% dplyr::summarise(avg.metabolite = mean(value)) %>%
  pivot_wider(names_from = metabolite, values_from = avg.metabolite)

p.z.c =
df.avg %>% 
  ggplot(aes(x = df.avg$B.caryophyllene+1 , y = df.avg$zingiberene+1))+
  geom_point()+
  #geom_text_repel(aes(label = df.avg$genotype))+
  geom_smooth(method = "lm", formula = y ~ x, alpha  = 0.2)+
  stat_cor(
    method = "pearson",
    label.x.npc = 0, 
    label.y.npc = 0.9)+
   scale_x_log10()+
   scale_y_log10()+
  xlab("Log-scaled B-caryophyllene (ng / type-VI gland)")+
  ylab("Log-scaled 7-epizingiberene (ng / type-VI gland)")+
  my.theme

ggsave(file = "Figure_6/plots/zingiberene_vs_caryophyllene.pdf", plot = p.z.c, width = 4, height = 3)
