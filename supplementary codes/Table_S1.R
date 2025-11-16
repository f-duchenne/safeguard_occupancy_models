############# FIRST ASSEMBLE ALL SPECIES TOGETHER
pkgs <- c("tidyr", "dplyr") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)


## Table S1
species <- read.csv("data/final_and_intermediate_outputs/database_clean_filtered.csv", stringsAsFactors = T) ## 4158540 records.

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

Trends <- read.csv("data/final_and_intermediate_outputs/all_trends.csv", stringsAsFactors = T)

Trends_1921 <- Trends %>% filter( baseline == "1921")

# Pivot the trends table so each region is a column

Trends_1921_clean <- Trends_1921 %>%
  filter(
    convergence == 0,        # keep only converged models
    !is.na(trend),           # remove missing trends
    trend >= -1 & trend <= 1 # keep plausible trends
  )

trend_summary <- Trends_1921_clean %>%
  select(species, region_50, trend) %>%
  pivot_wider(
    names_from = region_50,
    values_from = trend
  )


# Merge with the previous summary table
final_summary <- summary_table %>%
  left_join(trend_summary, by = c("scientificName" = "species"))

final_summary

#write_csv(final_summary, "data/final_and_intermediate_outputs/Table_S1.csv")

