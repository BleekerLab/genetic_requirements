###########
# The checkpoint library checkpoint allows you to install packages 
# as they existed on CRAN on a specific snapshot.
# It ensures script reproducibility
# More info: https://rdrr.io/cran/checkpoint/
###########

library("checkpoint")
#checkpoint("2020-01-01") # all package versions are from that date

library(tidyverse)
library(Rmisc)
library(multcompView)
library(stats)

##############
# custom theme
##############
my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1),
        axis.text.y = element_text(size = 8, colour = "black"),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.text.x = element_text(size=8, colour = "black")
  )+
  theme_bw()

#####################################################################################
# Figure 1A: 7-epizingiberene levels (norm. peak area) in F1s, parents, etc.
#####################################################################################
terpenes <- read.delim("Figure_1/20170425_leaf_washes.tsv", 
                       header = T, 
                       sep = "\t", 
                       dec = ".", 
                       check.names = F,
                       quote = "", # to avoid quotes around terpene names
                       stringsAsFactors = F)

# make a new column and calculate total terpene ion counts
terpenes$total_terpenes = rowSums(terpenes[,3:ncol(terpenes)]) 


# change to long ('tidy') format -> filter unwanted data (e.g. the LA1777-F1 measurements)
terpenes_long <- gather(data = terpenes, 
                        key = metabolite, 
                        value = abundance, 
                       -sample,
                       -genotype,
                       # metabolite converted to factor 
                       # (preserves original col order) 
                       factor_key=TRUE) %>% 
  mutate(genotype = factor(genotype, 
                           levels = c("CV", "F1", "PI127826",  "F1_hab", "LA1777"),
                           ordered = TRUE)) %>% 
  # keep only metabolites to be plotted
  filter(metabolite %in% c("total_terpenes", "7-epizingiberene")) 


##############################
# Figure 1A: metabolite levels
##############################
# summarise data for barplot
# summarySE: gives count, mean, standard deviation, 
# standard error of the mean, and confidence interval (default 95%). 
terpenes_mean_std = summarySE(terpenes_long, 
                              measurevar = "abundance", 
                              groupvars = c("genotype", "metabolite"))

#################
## Post-hoc tests
#################

# log transform the data (epizingiberene distribution not gaussian)
terpenes[terpenes == 0] <- 1 # to avoid Inf values 
terpenes_log10 = terpenes %>% 
  select(sample, genotype) %>% 
  cbind.data.frame(
    log10(terpenes[,3:ncol(terpenes)]))

# zingiberene
aov_genotype = aov(data = terpenes_log10,
                   formula = `7-epizingiberene` ~ genotype)
HSD_res <- TukeyHSD(aov_genotype, "genotype", ordered = TRUE)
groups <- as.data.frame(multcompLetters(HSD_res$genotype[,4])$Letters)
colnames(groups) = "7-epizingiberene"
groups$genotype = row.names(groups)
groups_zingi = pivot_longer(
  data = groups,
  cols = "7-epizingiberene", 
  names_to = "metabolite", 
  values_to = "group")


# total terpenes
aov_genotype_terpenes = aov(data = terpenes_log10,
                            formula = total_terpenes ~ genotype)
HSD_res2 <- TukeyHSD(aov_genotype_terpenes, "genotype", ordered = TRUE)
groups2 <- as.data.frame(multcompLetters(HSD_res2$genotype[,4])$Letters)
colnames(groups2) = "total_terpenes"
groups2$genotype = row.names(groups2)
groups_total = pivot_longer(
  data = groups2,
  cols = "total_terpenes", 
  names_to = "metabolite", 
  values_to = "group")

# collect letters in one unique dataframe
group_letters = rbind(groups_zingi,groups_total)


#######################
# Final plot: Figure 1A
#######################

# calculate maximum value + small margin to positionate the HSD letters
y_max <- max(terpenes_mean_std$abundance) + max(terpenes_mean_std$sd) 

p.figure1a =
  terpenes_mean_std %>% 
  ggplot(aes(x = genotype, y = abundance)) +
  geom_bar(stat = "identity", color = "black", fill = "black", alpha = 0.5) +
  geom_point(data = terpenes_long) +
  geom_errorbar(aes(x = genotype, 
                    ymin = abundance - se, 
                    ymax = abundance + se), 
                width = .2) +
  my.theme +
  facet_grid(~ metabolite) +
  ylab("Log10 of normalised metabolite abundance (ion counts / mg fresh leaf)")+
  xlab(NULL) +
  coord_flip() + 
  geom_text(data = group_letters, 
            aes(x = genotype, y = y_max, label = group))

p.figure1a
