library(tidyverse)
library(RColorBrewer)
library(multcompView)
library(Rmisc)
library(FSA)
library(ggpubr)

#############################
# Custom theme for plotting #
#############################

my.theme = 
  theme_bw() + 
  theme(text = element_text(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.text.x = element_text(size=10, colour = "black")
  )

##############
# Figure S3A #
#############
df = read.delim("Figure_03_trichome_vs_7epi/volatiles_and_trichomes.tsv", 
                header = T, 
                stringsAsFactors = F)
df.t <- 
  df %>% 
  mutate(density = ifelse(type_VI == 1, 0,
                          ifelse(type_VI == 2, 7.5,
                                 ifelse(type_VI == 3, 22.5,
                                        ifelse(type_VI == 4, 40, 50
                                        )))))

more.replicates <- df.t %>% dplyr::count(sample) %>% arrange(n) %>% filter(n != 2) %>% .$sample

df.t2 <-
  df.t %>%
  filter(!sample %in% more.replicates) %>% 
  group_by(sample) %>% 
  dplyr::summarise(ab_ad_density = sum(density))


df.zingiberene <- 
  df %>%
  group_by(sample) %>% 
  dplyr::summarise(zingiberene = mean(zingiberene))

df.fused <- 
  left_join(df.t2, df.zingiberene, by = "sample") %>% 
  mutate(class = ifelse(ab_ad_density %in% c(0, 7.5), "0-10 trichomes",
                        ifelse(ab_ad_density %in% c(15), "11-20 trichomes",
                               ifelse(ab_ad_density %in% c(22.5, 30), "21-30 trichomes",
                                      ifelse(ab_ad_density %in% c(40), "31-40 trichomes",
                                             ifelse(ab_ad_density %in% c(45,47.5), "41-50 trichomes",
                                                    ifelse(ab_ad_density %in% c(50,57.5), "51-60 trichomes",
                                                           ifelse(ab_ad_density %in% c(62.5), "61-70 trichomes",
                                                                  ifelse(ab_ad_density %in% c(72.5, 80), "71-80 trichomes",
                                                                         ifelse(ab_ad_density %in% c(90),"81-90 trichomes", 
                                                                                "90-100+ trichomes"
                                                                                
                                                                         ))))))))))


df.fused$class <- factor(df.fused$class, levels = unique(df.fused$class)[c(10,9,2,8,5,7,6,3,4,1)], ordered = TRUE)
df.fused$class <- gsub(" trichomes", "", df.fused$class)

df.fused %>%
  dplyr::count(class) %>%
  ggplot(aes(x = class, y = n))+
  geom_col(fill = "grey", color = "black") +
  ylab("Numer of F2 plants")+
  xlab("Trichomes per leafdisc")+
  geom_text(aes(label = paste(n)), position =   position_dodge(1), vjust = -0.5, size = 3)+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))+
  ylim(0,75) +
  theme_bw()+
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 8, colour = "black", angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(filename = "trichome_density_distribution.svg", height = 6, width = 11, units = "cm")



##############
# Figure S3B #
##############
# 7-epiZ versus trichome densty-class
# The genotypes with PI127826 zingiberene levels marked in red.

give.n <- function(x){
  return(c(y = 18, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}
df.fused %>% 
  mutate(zingiberene = ifelse(zingiberene ==0 & class == "0-10", 1, zingiberene)) %>%
  mutate(selected = ifelse(zingiberene != 0 ,1, ifelse(class=="0-10", 1, 0))) %>%
  mutate(chemotype = ifelse(zingiberene > (514142-159176), "PI127826", "Elite")) %>%
  # filter(selected == 1) %>%
  ggplot(aes(x = class, y = log(zingiberene)))+
  geom_boxplot(fill = "grey")+
  geom_jitter(width = 0.05, aes(color = chemotype))+
  stat_summary(fun.data = give.n, geom = "text", fun = median,
               position = position_dodge(1), size = 3)+
  ylab("7-epizingiberene (log2 ion-counts")+
  xlab("Trichomes-density class")+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
  scale_color_manual(values = c("black", "red"))+
  theme_bw()+
  #  ylim(0,40)+
  scale_fill_brewer(palette = "Greys")+
  theme(axis.text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, colour = "black", angle = 45, hjust = 1),
        axis.title = element_text(size = 10, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

ggsave(filename = "figures/Figure_S2_Full_F2_density_classes_and_7epiZ/Fig_S2_density_vs_zingiberene.svg", height = 6, width = 11, units = "cm")
