library(tidyverse)
library(ggbiplot)
library(Rmisc)
library(multcompView)

#Theme for plotting
my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 8, colour = "black", angle = 30, hjust = 1),
        axis.text.y = element_text(size = 8, colour = "black"),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.text.x = element_text(size=8, colour = "black")
  )+
  theme_bw()


#Create_dataframe
df <- read.delim("Figure_1_F1_phenotypes/20170425_F1s leaf washes_peak area.txt", header = T, sep = "\t", dec = ".", check.names = F)
df$total_terpenes = rowSums(df[,3:ncol(df)]) #make a new column and calculate total terpene ion counts


#Change to long ('tidy') format -> filter unwanted data (e.g. the LA1777-F1 measurements)
df_long <- gather(df, metabolite, abundance, Bpinene:total_terpenes, -genotype,-genotype,factor_key=TRUE) %>% 
  filter(., genotype != "CV_LA1777") %>% filter(., genotype != "F1_LA1777")



#Define variables 
df_long$genotype <- factor(df_long$genotype, levels = c("CV", "F1", "PI127826",  "F1_hab", "LA1777"),
                              ordered = TRUE)

#summarise data for barplot
sum = summarySE(df_long, measurevar = "abundance", groupvars = c("genotype", "metabolite"))

#### BARPLOT ###

p.leafwash =
sum %>% filter(., metabolite %in% c("total_terpenes", "7epizingiberene")) %>%
ggplot(., aes(x = genotype, y = abundance)) +
  geom_bar(stat = "identity", color = "black", fill = "black") +
    geom_errorbar(aes(x = genotype, ymin=abundance-se, ymax=abundance+se), width=.2) +
    my.theme +
  facet_grid(~metabolite)+
  ylab("Metabolite abundance (ion counts / mg fresh leaf)")+
  xlab(NULL)

ggsave("Figure_1_F1_phenotypes/plots/F1_leafwashes_barplot.pdf", plot = p.leafwash, width = 6, height = 4)

#check for normality
shapiro.test(df$total_terpenes) 
shapiro.test(df$X7epizingiberene)

#Log transform the data
df[df == 0] <- runif(n=1,min = 1, max = 2)
df.log =  as.data.frame(c(df[1:2],log(select(df, -c(sample, genotype))))) %>% 
  filter(., genotype != "CV_LA1777") %>% filter(., genotype != "F1_LA1777")


aov_genotype = aov(X7epizingiberene ~ genotype, data=df.log)
summary(aov_genotype)
TK <- TukeyHSD(aov_genotype, "genotype", ordered = TRUE)

# Write the letters
multcompLetters(TK$genotype[,4])

#Write files with results
write.csv(TukeyHSD_results, file = "F1_leafwash_TukeyHSD_results.csv")


############
#   PCA    #
############

df2 = df %>% column_to_rownames(., var = "sample") %>% select(., -total_terpenes) %>%
  filter(., genotype != "CV_LA1777") %>% filter(., genotype != "F1-LA1777")

genotypes= df2[,1]
df2 = df2[,2:ncol(df2)]+1
df.log = log(df2)

df.pca = prcomp(df.log)
summary(df.pca)

# Plot PCA

my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 8, angle = 0, colour = "black"), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.text.x = element_text(size=8, colour = "black")
  )+
  theme_bw()

p.pca = 
ggbiplot(df.pca, ellipse = TRUE, groups = genotypes, var.axes = FALSE)+my.theme

ggsave(filename = "Figure_1_F1_phenotypes/plots/PCA_F1_leafwash.pdf", plot = p.pca, width = 12, height = 16, units = "cm")
