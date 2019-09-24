library(tidyverse)
library(Rmisc)
library(car)
library(ggpubr)
library(multcompView)
library(rcompanion)

#########################
# Import and shape data #
#########################

#Load dataset -> order genotypes -> make data "tidy"
df = read.csv(file = "Figure_1_F1_phenotypes/Whitefly_survival_F1s.csv",
              header = T, stringsAsFactors = T)
df$plant = as.factor(df$plant)
df$genotype = factor(df$genotype, 
                     levels = c("Cultivar", "PI127826", "F1", "PI127826 x LA1777", "LA1777"),
                     ordered = TRUE)

df.wide  = spread(data = df,
                  key = "clipcage",
                  value = "survival")

#take the average of the two clip cages
df.wide$survival = (df.wide$A + df.wide$B)/2

#summarise for the barplot - calculate averages per genotypem + se
sum = summarySE(df.wide, 
                measurevar = "survival",
                groupvars = "genotype")



########
# Plot #
########

# Costum theme for plotting
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

# Make barplot

p.wf = 
ggplot(sum, aes(x = genotype,
                y = survival)) +
  geom_bar(stat = "identity",
           fill = "black") +
  geom_errorbar(aes(ymin = survival-se, ymax = survival + se),
                width = 0.1)+
  my.theme +
  ylab("Whitefly survival (%)")+
  xlab(NULL)

########
# Save #
########

ggsave("Figure_1_F1_phenotypes/plots/wf_survival.svg", plot= p.wf, device = "svg")

##############
# Statistics #
##############

#check for normality
shapiro.test(df.wide$survival) 
qqPlot(df.wide$survival)

#non-parametric test
kruskal.test(data = df, survival ~ genotype)
wf.wilcox = pairwise.wilcox.test(df$survival, df$genotype, p.adjust.method = "none", paired = FALSE)
capture.output(wf.wilcox, file = "Figure_1_F1_phenotypes/plots/wf_phenotype_wilcox.txt")

letters_sig = multcompLetters(fullPTable(wf.wilcox$p.value), #compare the groups
                              compare = "<",
                              threshold = 0.05)
print(letters_sig)
