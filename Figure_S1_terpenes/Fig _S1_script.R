library(tidyverse)
library(Rmisc)

#############################
# Costum theme for plotting #
#############################

my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, colour = "black"),
        axis.text.y  = element_text(size = 11, angle = 0, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        strip.text =  element_text(size = 12, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.text.x = element_text(size=12, colour = "black")
  )


#########################
# Load and format data #
#########################

#Create_dataframe
df <- read.delim("Figure_1_terpenes_type6_CVxPI127826/20200722_F1s_leafwash_ng_mg_tissue.txt", header = T, sep = "\t", dec = ".", check.names = FALSE)
df[,3:ncol(df)] = lapply(df[,3:ncol(df)], as.numeric)


#Change to long ('tidy') format -> filter unwanted data (e.g. the LA1777-F1 measurements)
df_long <- df %>% pivot_longer(., cols = -c(genotype, sample, replicate), names_to = "metabolite", values_to = "abundance") %>% 
  filter(genotype %in% c("Elite", "PI127826", "F1", "F1_hab", "LA1777"))



#Define variables 
df_long$genotype <- factor(df_long$genotype, levels = c("Elite", "F1", "PI127826", "F1_hab", "LA1777"),
                              ordered = TRUE)

#summarise data for barplot
sum = summarySE(df_long, measurevar = "abundance", groupvars = c("genotype", "metabolite"))

# Voeg classes toe aan de metaboliten

terpene2class <- openxlsx::read.xlsx("Figure_S1_terpenes/terpene2class.xlsx")

sum <- left_join(sum, terpene2class, by = "metabolite")

sum$class <- factor(sum$class, levels = c("monoterpene", "plastidial sesquiterpene", "cytosolic sesquiterpene"),
                    ordered = T)



sum$metabolite <- gsub('_', '-', sum$metabolite)
sum$genotype <- gsub('_', '-', sum$genotype)

sum$genotype <- factor(sum$genotype, levels = c("Elite", "F1", "PI127826", "F1-hab", "LA1777"),
                       ordered = TRUE)
########### 
# BARPLOT #
###########

p.leafwash =
sum %>% 
filter(metabolite != "summed-terpenes") %>%
ggplot(., aes(x = genotype, y = abundance)) +
  geom_bar(stat = "identity", color = "black", aes(fill = class)) +
    geom_errorbar(aes(x = genotype, ymin=abundance-se, ymax=abundance+se), width=.2) +
    my.theme +
  facet_wrap(class~metabolite, ncol =4, scale = "free")+
  scale_fill_manual(values = c("Black", "Grey", "white"))+ # colors c("#002a5c", "#e9500e", "#42b8b2")
  ylab("Metabolites (ng / mg fresh leaf tissue)")+
  xlab(NULL)+
  theme_bw()+
  my.theme+
  theme(legend.position = "none")

 ggsave(filename = "Figure_S1_terpenes/Fig_S1_leafwash_3.svg", plot = p.leafwash, width = 210, height = 297, units = "mm")

p.leafwash.stacked =
  sum %>% filter(metabolite != "summed_terpenes") %>%
  ggplot(., aes(x = genotype, y = abundance)) +
  geom_bar(aes(x = genotype, y = abundance, fill = metabolite), stat = "identity", position = "stack") + 
  my.theme +
  theme(legend.position = "none")

############### 
# STATISCTICS #
###############

df.wide = sum %>% select(genotype, metabolite, abundance) %>% pivot_wider(names_from = metabolite, values_from = abundance)

df.percentage <- (df.wide[,2:ncol(df.wide)] / df.wide$summed_terpenes) * 100
rownames(df.percentage) <- df.wide$genotype
df.percentage = round(df.percentage, digits = 1)
df.percentage = rownames_to_column(df.percentage, var = "genotype")

write.table(df.percentage, file = "Figure_S1/terpene_contributions_percentage.txt", row.names = FALSE, sep = "\t")
