library(tidyverse)

df = read.delim(file = "Figure_7/abundance_tidy.tsv", header = T)

# Shorten target_id names. Then Remove LA1777  related samples from df
df.parsed = df %>% mutate(target_id = substr(target_id, start = 1, stop = 14)) %>%
  filter(sample %in% c("Elite_01", "PI127826_F1", "F2-151", "F2-411", "F2-445",
                       "PI127826", "F2-28", "F2-73", "F2-127"))

# Custom order of the samples
df.parsed$sample = factor(df.parsed$sample, levels = c("Elite_01", "PI127826_F1", "F2-151", "F2-411", "F2-445",
                                                                   "PI127826", "F2-28", "F2-73", "F2-127"),
                          ordered = TRUE)

# Create df with active/lazy condition per samples
con = data.frame(sample = c("Elite_01", "PI127826_F1", "F2-151", "F2-411", "F2-445",
                            "PI127826", "F2-28", "F2-73", "F2-127"),
                 condition = c("lazy","lazy","lazy","lazy","lazy",
                               "active","active","active","active"))

# Fuse the active/lazy condition with the main df
df.parsed = left_join(df.parsed, con, by = "sample")

# Filter by solycnumber and plot
df.parsed %>% filter(target_id == "Solyc07g052250") %>%
  ggplot(aes(x = reorder(sample, condition), y = est_counts, fill = condition))+
  geom_bar(stat = "identity")

# Plot per condition (mean active vs. mean lazy)
df.parsed %>% filter(target_id == "Solyc09g042290") %>%
  group_by(condition) %>% dplyr::summarise(mean = mean(est_counts), sd = sd(est_counts)) %>%
  ggplot()+
  geom_bar(aes(x = condition, y = mean), stat = "identity")+
  geom_errorbar(aes(x = condition, ymin = mean-sd, ymax = mean+sd))
