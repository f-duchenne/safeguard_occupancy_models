## Table S1

library(dplyr)
library(tidyr)

species <- read.csv("database_clean_filtered.csv", stringsAsFactors = T) ## 4158540 records.

summary_table <- species %>%
  group_by(family, scientificName, taxo_group,time_period) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = time_period,
    values_from = n,
    values_fill = 0
  )

summary_table

colnames(summary_table)[3] <- "Group"

str(Trends_1921)


# Pivot the trends table so each region is a column

Trends_1921_clean <- Trends_1921 %>%
  filter(
    convergence == 0,        # keep only converged models
    !is.na(trend),           # remove missing trends
    trend >= -1 & trend <= 1 # keep plausible trends
  )

trend_summary <- Trends_1921_clean %>%
  select(Species, region_50, trend) %>%
  pivot_wider(
    names_from = region_50,
    values_from = trend
  )

library(dplyr)

# Merge with the previous summary table
final_summary <- summary_table %>%
  left_join(trend_summary, by = c("scientificName" = "Species"))

final_summary

write_csv(final_summary, "Table_S1.csv")

