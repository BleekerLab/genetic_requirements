df = read.xlsx("GCMS_Leafwash_Full_F2_EZ.xlsx",1,header = TRUE)

#Filter datafile on zingiberene and type_VI_trichomes both abaxial and adaxial
df.parsed = filter(df, df$group == "F2")
head(df.parsed)

ggplot(df.parsed, aes(x= df.parsed$pinene, y=df.parsed$zingiberene)) +
  geom_point() +
  scale_x_continuous(trans='log10') + #Log10 scale
  scale_y_continuous(trans='log10') + #Log10 scale
  geom_smooth(method = "lm", formula = y ~ x) + 
  theme_classic() +
  theme_linedraw()+
  theme(text = element_text(family = "Times New Roman", color = "black", size = 10))+
  labs(y = "7-epizingiberene (ion counts)", x = "alpha-pinene (ion counts)") +
  #geom_text(aes(label=df.parsed$sample),hjust=1.5, vjust=0.5, size = 3) +  #Adds labels      
  ggsave("Full_F2_Zingiberen vs pinene.jpg", device = "jpg", scale = 1, width = 12, height = 8, units = "cm", dpi = 300)