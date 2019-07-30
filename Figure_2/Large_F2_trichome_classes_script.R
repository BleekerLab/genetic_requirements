library(ggplot2)
library(ggalt)
library(tidyverse)
library(plotly)

df = read.csv(file = "trichomes full F2 PI127826 CSV.csv", sep = ";", header = T)
str(df)

# Scatterplot type-VI vs Type-I/IV
typeVI_I_IV_scatter =
df %>% #filter(., background == "F2") %>%
  
  ggplot(., 
         aes(
         x = sum_type_VI,
         y = TypeIandIV
         )
         ) +
  
  geom_point(position = "jitter", 
             aes(colour = background),
             alpha = 1) +
  geom_smooth(method = "lm")+
  facet_wrap(~background)+
  theme_bw()
plot(typeVI_I_IV_scatter)

# 3D using plotly
plot_ly(df, x=df$TypeIandIV, 
            y=df$sum_type_VI, 
            z=df$NonGlandular, 
        type="scatter3d", 
        mode="markers", 
        color=df$background)

# Scatterplot type-VI vs Non-glandular
typeVI_non_glandular_scatter = 
  df %>% #filter(., background == "F2") %>%
  
  ggplot(., 
         aes(
           x = sum_type_VI,
           y = NonGlandular
         )
  ) +
  
  geom_point(position = "jitter", 
             aes(colour = background),
             alpha = 1) +
  geom_smooth(method = "lm")+
  facet_wrap(~background)+
  theme_bw()
plot(typeVI_non_glandular_scatter)

# Histogram of the different trichome classes
df.long = 
  gather(df,
         key = "trichome_type",
         value = "class",
         -plant, -surface, -background, -name, -date)
df.long$class = as.integer(df.long$class)

df.long %>% filter(., trichome_type %in% c("NonGlandular", "sum_type_VI", "TypeIandIV")) %>%
ggplot(., aes(x = class)) +
  geom_histogram(fill = "grey")+
  facet_grid(background~trichome_type, scale = "free")+
  #geom_density(alpha=.2, fill="#FF6666") +
  theme_bw()


