library(tidyverse)
library(Rmisc)

#############################
# Costum theme for plotting #
#############################

my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1, colour = "black"),
        axis.text.y  = element_text(size = 6, angle = 0, colour = "black"),
        strip.text =  element_text(size = 6, colour = "black"),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.text.x = element_text(size=8, colour = "black")
  )+
  theme_bw()


#########################
# Load and format data #
#########################

#Create_dataframe
df <- read.delim("Figure_1_F1_phenotypes/20170425_F1s leaf washes_peak area.txt", header = T, sep = "\t", dec = ".")
df$total_terpenes = rowSums(df[,3:ncol(df)]) #make a new column and calculate total terpene ion counts
df[,3:ncol(df)] = lapply(df[,3:ncol(df)], as.numeric)


#Change to long ('tidy') format -> filter unwanted data (e.g. the LA1777-F1 measurements)
df_long <- gather(df, metabolite, abundance, Bpinene:total_terpenes, -genotype,-genotype,factor_key=TRUE) %>% 
  filter(., genotype != "CV_LA1777") %>% filter(., genotype != "F1_LA1777")



#Define variables 
df_long$genotype <- factor(df_long$genotype, levels = c("CV", "PI127826", "F1", "F1_hab", "LA1777"),
                              ordered = TRUE)

#summarise data for barplot
sum = summarySE(df_long, measurevar = "abundance", groupvars = c("genotype", "metabolite"))

########### 
# BARPLOT #
###########

p.leafwash =
sum %>% filter(metabolite != "total_terpenes") %>%
ggplot(., aes(x = genotype, y = abundance)) +
  geom_bar(stat = "identity", color = "black", fill = "black") +
    geom_errorbar(aes(x = genotype, ymin=abundance-se, ymax=abundance+se), width=.2) +
    my.theme +
  facet_wrap(~metabolite, ncol =4, scale = "free")+
  ylab("Metabolite abundance (ion counts / mg fresh leaf)")+
  xlab(NULL)

ggsave(filename = "Figure_S1/Fig_S1_leafwash.svg", plot = p.leafwash, width = 12, height = 12)

############### 
# STATISCTICS #
###############
df.log = df[df == 0] <- runif(n=1,min = 1, max = 2)
df.log = df.log =  as.data.frame(c(df[1:2],log(select(df, -c(sample, genotype)))))
apply(df.log[,3:ncol(df)], 2, shapiro.test)


