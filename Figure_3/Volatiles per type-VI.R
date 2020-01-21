library(tidyverse)
library(RColorBrewer)
library(multcompView)
library(Rmisc)
#############################
# Custom theme for plotting #
#############################

my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1, colour = "black"),
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

volatiles_VI = read.csv(file = "Figure_2/type_VI_gland_terpenes_F2.csv", header = T, stringsAsFactors = T, check.names = F) %>% filter(., !genotype %in% c("LA1777", "PI127826xLA1777", "LA1777_F1", "CV_LA1777"))
trichomes = read.csv(file = "Figure_3/trichome_density_F2_selection.csv", header = T, stringsAsFactors = T, check.names = F)

########################### 
# Zingiberene per type- VI#
###########################

# Calculate the mean / sd of the volatiles and write to file

write.table(volatiles_VI %>% dplyr::group_by(.,genotype) %>% dplyr::summarise_at(c(3,4,5,6,7,8), funs(mean, sd)),
            file = "Figure_3/mean_volatiles_type_VI_glands.tsv", row.names = F, sep = "\t", dec = ".")


# Summarise data and make a barplot
sum.volatiles = summarySE(volatiles_VI, measurevar = "total_volatiles", groupvars = c("genotype", "group"))

p.volatiles = 
ggplot(sum.volatiles, x = reorder(genotype, -total_volatiles), y = total_volatiles)+
  geom_bar(stat = "identity", aes(x =reorder(genotype, -total_volatiles), y = total_volatiles, fill = group))+
  geom_errorbar(aes(x = genotype, ymin = total_volatiles - se, ymax = total_volatiles + se))+
  scale_fill_manual(values = c("gray28", "lightgrey", "gray52", "black"))+
  xlab(NULL)+
  ylab("Summed terpenes (ng / type-VI trichome)")+
  my.theme + theme(legend.position = "none")

ggsave("Figure_3/barplot_total_volatiles_type_VI.svg", plot = p.volatiles, height = 4, width = 7)


###############################
# Distribution of zinigberene #
###############################


#make the density function of the data
#scaled
zingi.scaled = approxfun(density(scale(log(volatiles_VI$total_volatiles+1))))

#log

vol.log = approxfun(density((log(volatiles_VI$total_volatiles+1))))
pdf("Figure_3/F2_volatiles_type_VI_density_plot.pdf") 
plot(vol.log, xlim = c(-1.5,4.5), xlab = "scaled zingiberene content / type-VI trichome", ylab = "proportion of the F2 population")
abline(v = log(17.7+1)) # mean PI127826
abline(v = log(0.51+1)) # mean CV
abline(v = log(0.9 +1)) # mean F1
dev.off()

#Calculate the area's of the function: integrate(funtion, xmin, xmax)
integrate(zingi.scaled, -2.6,0.91)

#plot the data

##########################
# Density vs. volatiles #
#########################

mean.trichomes = trichomes %>% dplyr::group_by(., genotype) %>% dplyr::summarise(., mean(type_VI_mm2))
mean.volatiles = volatiles_VI %>% dplyr::group_by(., genotype) %>% dplyr::summarise(., mean(total_volatiles))

mean.trichomes.sum.volatiles = inner_join(mean.trichomes, mean.volatiles, by = "genotype")

p.density.activity = 
ggplot(mean.trichomes.sum.volatiles,
       aes(x = mean.trichomes.sum.volatiles$`mean(type_VI_mm2)`, 
           y = mean.trichomes.sum.volatiles$`mean(total_volatiles)`))+ 
  geom_point()+
  ylab("Volatiles per type-VI gland (ng / gland)")+
  xlab("Type-VI trichome density (trichomes / mm2)")+
  geom_text(aes(label=mean.trichomes.sum.volatiles$genotype),hjust=1.5, vjust=0.5, size = 2)+
  my.theme

ggsave(file = "Figure_3/desnty_VS_acitivty.pdf", plot = p.density.activity, height = 4, width = 4)





##############
# Statistics #
##############

sum.volatiles$total_volatiles = as.integer(sum.volatiles$total_volatiles)
glm(data=sum.volatiles, total_volatiles ~ genotype, family = poisson(link = "log"))

coefs = as.data.frame(summary(fit)$coefficients)
colnames(coefs)=c("coeff","stderr","zval","pval")
coefs$genotype = row.names(coefs$genotype)
coefs$genotype= sub("^genotype",replacement = "",x = row.names(coefs))
coefs$genotype[1] = "Intercept"

coefs$text = NA
coefs[which(coefs$pval > 0.05),]$text <- "ns"
coefs[which(coefs$pval < 0.05),]$text <- "***"

write.table(coefs,file = "Figure_3/GLM_type_VI_volatiles.tsv",sep = "\t",quote = F,row.names = F)

