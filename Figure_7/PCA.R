library(tidyverse)
library(ggbiplot)
library(ggfortify)

#######
# PCA #
#######
# Load data 

df.for.pca <- df.wide %>% column_to_rownames(var = "target_id")

#Select for top 1000 highest expressed genes
df.for.pca$sum <- rowSums(df.for.pca)
df.for.pca = df.for.pca[order(-df.for.pca$sum),]
df.for.pca.top1000 = df.for.pca[1:10000,] %>% select(-sum)

df.for.pca <- t(df.for.pca.top1000)
df.for.pca <- log(df.for.pca +1)

pca <- prcomp(df.for.pca, scale = T, center = T)

conditions = c("lazy", "active", "lazy", "active", "lazy", "lazy", "active", "active", "lazy")

# ggbiplot(pca, groups = conditions, obs.scale = 1, var.scale = 1)
autoplot(pca, label = TRUE)
