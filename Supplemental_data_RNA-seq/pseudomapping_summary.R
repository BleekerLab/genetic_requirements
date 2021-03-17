library("checkpoint")
checkpoint("2020-01-01")
library("rjson")
suppressPackageStartupMessages(library("tidyverse"))

# Import JSON pseudomapping statistics
dirs_with_json_files <- list.dirs(path = "Supplemental_data_RNA-seq/kallisto", full.names = TRUE, recursive = FALSE)
json_files <- sapply(dirs_with_json_files, function(x){paste0(x,"/run_info.json")}) %>% 
  as.vector()
samples <- sapply(dirs_with_json_files, function(x){basename(x)})

dfs <- lapply(json_files, function(x){fromJSON(file = x)})
names(dfs) <- samples

# mapping summary df
pseudo_mapping_results <- bind_rows(dfs) %>% 
  mutate(sample = samples)
write.csv(x = pseudo_mapping_results, file = "Supplemental_data_RNA-seq/mapping_results.csv")


# Plot pseudomapping results
plot_pseudomapping <- 
  ggplot(pseudo_mapping_results, 
         aes(x = sample, y = p_pseudoaligned)) +
  geom_bar(stat = "identity") +
  labs(x = "Tomato genotype", y = "Percentage of pseudoaligned reads (%)") +
  scale_y_continuous(limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
plot_pseudomapping
ggsave(filename = "Supplemental_data_RNA-seq/mapping_summary.png", plot = plot_pseudomapping)
