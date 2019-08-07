library(Rmisc)
library(tidyverse)
library(ggplot2)

df = read_csv(file = "20190807_cavity_volumes_selected_Active_Lazy_F2.csv")

df = df %>% filter(., genotype != "F1_hab") #Remove F1_hab from the dataset

df$genotype = factor(df$genotype, levels = c("CV", "PI127826","F1", "151", "411","445","28","73", "127"),
                                             ordered = TRUE)

sum = summarySE(df, measurevar = "volume_um", groupvars = c("genotype", "phenotype"))
write.table(sum, file = "cavity_volume_F2_summary_table.txt", sep = "\t", row.names = FALSE)

ggplot(sum, aes(x = genotype, y = volume_um, fill = phenotype)) +
  geom_bar(aes(x = genotype, y = volume_um), stat = "identity") + 
  geom_errorbar(aes(ymin = volume_um - se, ymax = volume_um + se),
                width = 0.1)+
  theme_bw()


### Volatle in type VI trichomes
df2 = read.csv(file = "20190807_Type_VI_volatiles_High_Low_F2s.csv",
                header = TRUE)
df2 = df2 %>% filter(., genotype != "F1_hab")
df2$genotype = factor(df2$genotype, levels = c("CV", "PI127826","F1", "151", "411","445","28","73", "127"),
                     ordered = TRUE)


# Define what you want to plot / calculate statistics of (e.g. "total_terpenes")
sum2 = summarySE(df2, measurevar = "total_terpenes", groupvars = c("genotype", "phenotype"))
write.table(sum2, file = "Total_terpenes_F2_summary_table.txt", sep = "\t", row.names = FALSE)


ggplot(sum2, aes(x = genotype, y = total_terpenes, fill = phenotype)) +
  geom_bar(aes(x = genotype, y = total_terpenes), stat = "identity") + 
  geom_errorbar(aes(ymin = total_terpenes - se, ymax = total_terpenes + se),
                width = 0.1)+
  theme_bw()

## Combine
ggplot(data = NULL) +
  geom_bar(data = sum2, aes(x = genotype, y = total_terpenes), stat = "identity")+
  geom_bar(data = sum, aes(x = genotype, y = volume_um), stat = "identity")
