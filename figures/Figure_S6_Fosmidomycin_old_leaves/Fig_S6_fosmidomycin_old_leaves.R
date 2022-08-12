library(tidyverse)

# Load the data
df.raw <- read.csv("data/Fosmidomycin_14_Sofia.csv")

# The data gives cavity diameters -> Calculate volumes 
df <- df.raw %>% 
  mutate(volume = (4/3*pi*(Cavity_horizontal/2)*(Cavity_horizontal/2)*(Cavity_vertical/2))/1000)

# Summarise
df.sum <- df %>% 
  group_by(Treatment, Leaf) %>%
  dplyr::summarise(avg_volume = mean(volume),
                   se_volume = sd(volume)/sqrt(n()),
                   n = n())

# Plot
df.sum %>% 
 filter(Leaf == "Old") %>% 
  ggplot(aes(x = Treatment, y = avg_volume))+
  geom_bar(fill = "black", color = "black", stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin = avg_volume-se_volume,
                    ymax = avg_volume+se_volume,
                    x = Treatment), position = "dodge",
                width = 0.5)+
  ylab("Storage-cavity volume(picolitre)")+
  xlab(NULL)+
  theme_bw()+
  theme(axis.text = element_text(size = 12, color = "black")) +
  facet_grid(~Leaf)

ggsave("figures/Figure_S6_Fosmidomycin_old_leaves/Fig_S6.pdf")

##########
# T-test #
##########

df.test <- df %>% 
  filter(Leaf == "Old")

ggpubr::ggqqplot(df.test$volume)
shapiro.test(df.test$volume)

wilcox.test(volume ~ Treatment, data = df.test)

t.test(volume ~ Treatment, data = df.test )  

df.test %>% ggplot(aes(x = volume))+geom_density()+ facet_wrap(~Treatment, ncol = 1)+theme_bw()+xlim(0,89)
