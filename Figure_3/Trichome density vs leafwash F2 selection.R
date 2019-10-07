library(tidyverse)

#############################
# Custom theme for plotting #
#############################

my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 7, colour = "black"),
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

df = read.csv("Figure_3/Leafwash vs Trichome density.csv",header = TRUE, stringsAsFactors = T, check.names = FALSE)


#Change to long ('tidy') format. First the volatiles (df.long) and than also the trichomes (df.long2)
df.long = gather(
  data = df,
  key = "metabolite",
  value = level,
  -sample, -group, -Type_VI_mm2_Abaxial,-Type_VI_mm2_Adaxial,-Type_VI_mm)

df.long2 = gather(
  data = df.long,
  key = "trichome_position",
  value = density,
  -sample,-metabolite,-group,-level)


#Filter datafile on total_volatiles and type_VI_trichomes both abaxial and adaxial
df.parsed = df.long2 %>% filter(., df.long2$metabolite == "total_volatiles" & df.long2$trichome_position == "Type_VI_mm") %>% filter(., !sample == "250")
df.parsed$metabolite = as.factor(df.parsed$metabolite)
df.parsed$trichome_position = as.factor(df.parsed$trichome_position)

# linear model fitting the summed terpenes to the density of trichomes
model  = lm(data = df.parsed, level ~ density)
summary(model)

#Scatterplot
p.scatter = 
ggplot(df.parsed, aes(x=density, y=level)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, alpha  = 0.2) + 
  my.theme+
  ylab("Summed terpenes (ng / mg fresh leaf)")+
  xlab("Type-VI trichome density (trichomes / mm2)")+
  geom_text(aes(label=df.parsed$sample),hjust=1.5, vjust=0.5, size = 1.5)+
  annotate(geom = "text", x = 22, y = 320, 
           label = round(summary(model)$adj.r.squared, digits= 3), 
           size = 3) +
  annotate(geom = "text", x = 19.5, y = 320, 
           label = "r^2=",
           size = 3)

# Save plot
ggsave(file = "Figure_3/trichome_density_VS_leafwash.svg", plot = p.scatter, height = 4, width = 4)
  

