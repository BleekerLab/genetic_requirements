library(pheatmap)
library(tidyverse)

df <- read.delim(file = "Figure_7/precursors_allele_expression.txt",
                 header = TRUE, row.names = 1,
                 check.names = FALSE,
                 sep = "\t",) %>% select(-annotation, -pathway) 
df <- na.omit(df)
pheatmap(mat = df, 
         scale = "none",
         color = c("black", "grey", "white"),
         breaks = c(1, 1.9, 2,2.1, 3),
         cluster_rows = F, 
         cluster_cols = F,
         fontsize = 10,
         cellwidth = 15,
         cellheight = 15,
         gaps_row = 17,
         gaps_col = 5,
         filename = "Figure_7/heatmaps/allele_expression_MEP_MVA_heatmap.pdf"
         )

