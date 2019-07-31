df = read.csv(file = "Whitefly_survival_F1s.csv",
              header = T)
df$plant = as.factor(df$plant)
df$genotype = factor(df$genotype, 
                     levels = c("Cultivar", "PI127826", "F1", "PI127826 x LA1777", "LA1777"),
                     ordered = TRUE)
str(df)

df.wide  = spread(data = df,
                  key = "clipcage",
                  value = "survival")
head(df.wide)

#take the average of the two clip cages
df.wide$survival = (df.wide$A + df.wide$B)/2

#summarise for the barplot
sum = summarySE(df.wide, 
                measurevar = "survival",
                groupvars = "genotype")
head(sum)

ggplot(sum, aes(x = genotype,
                y = survival)) +
  geom_bar(stat = "identity",
           fill = "black") +
  geom_errorbar(aes(ymin = survival-se, ymax = survival + se),
                width = 0.1)+
  theme_bw()

# Non-parametric statistics?
kruskal.test(data = df, survival ~ genotype)
pairwise.wilcox.test(df$survival, df$genotype, p.adjust.method = "none", paired = FALSE)

