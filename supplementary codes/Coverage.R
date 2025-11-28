############# FIRST ASSEMBLE ALL SPECIES TOGETHER
pkgs <- c("tidyverse", "iNEXT") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)


## Check coverage for richness analyses.

# Load dataset
species <- read.csv("data/final_and_intermediate_outputs/database_clean_filtered.csv", stringsAsFactors = T) ## 4158540 records.

# Filter regions of the study
species <- species %>% filter (region_100 %in% c("alpine", "mediterranean", "continental", "atlantic", "boreal"))

# Filter to keep bees only
species <- species %>% filter (taxo_group %in% "bees") #Keep only bees. # 3494628 records.

# Count occurrences of each species per bioregion and year.
species_counts <- species %>%
  group_by(region_100, endYear, scientificName) %>%
  summarise(count = n(), .groups = "drop")

# Create a unique ID for each region-year combination.
species_counts <- species_counts %>%
  mutate(region_year = paste0("R", region_100, "_Y", endYear))

# Transpose to wide to get abundance of occurrences per species.
abundance_list <- species_counts %>%
  select(region_year, scientificName, count) %>%
  pivot_wider(names_from = scientificName, values_from = count, values_fill = 0) %>%
  column_to_rownames("region_year") %>%
  as.data.frame() %>%
  split(., rownames(.)) %>%
  lapply(as.numeric)

# Combine the list into a matrix
abundance_matrix <- do.call(cbind, abundance_list)

# Transpose so rows = samples, columns = species
abundance_transposed <- t(abundance_matrix)


### Check coverage for bees at n = 100----

sample_size <- 100

abundance_transposed_filtered <- abundance_transposed[rowSums(abundance_transposed) >= sample_size, ]

## Check coverage:

# Function: estimate coverage for each region-year at fixed n
estimate_coverage_fixed <- function(abundances, n = sample_size) {
  # abundances is a numeric vector of species counts
  abundances <- abundances[abundances > 0]
  
  # Skip samples with fewer than n individuals
  if (sum(abundances) < n) return(NA)
  
  # Estimate coverage for sample size n
  out <- iNEXT(list(data = t(abundances)), datatype = "abundance", q = 0, size = sample_size)
  return(out$iNextEst[[1]])  # Extract coverage value at n individuals
}

#slow
coverage_100 <- apply(
  abundance_transposed_filtered, 
  1, 
  estimate_coverage_fixed, 
  n = sample_size
)

coverage_summary <- bind_rows(coverage_100, .id = "sample") %>%
  filter(m == 100) %>%   # keep only the rows corresponding to the selected sample size.
  mutate(
    region_year = sample
  ) %>%
  separate(region_year, into = c("region_100", "endYear"), sep = "_Y", remove = FALSE) %>%
  mutate(
    region_100 = gsub("R", "", region_100),
    endYear = as.integer(endYear)
  )

coverage_summary %>% group_by(region_100) %>% summarise(mean = mean(SC))

# Define colors for each region
region_colors <- c(
  "alpine" = "#44AA99",       # Blue
  "atlantic" = "#117733",     # Orange
  "boreal" = "#332288",       # Green
  "continental" = "#CC6677",  # Red
  "mediterranean" = "#DDCC77" # Purple
)

p_coverage_bees_100 <- ggplot(coverage_summary, aes(x = endYear, y = SC, color = region_100)) +
  geom_point() +
  geom_line() +
  facet_wrap(. ~ region_100, nrow = 1) +
  scale_color_manual(values = region_colors) +
  ylim(0, 1) +
  labs(x = "Year", y = "Coverage (100)") +
  theme_minimal() +
  theme(legend.position = "none")

p_coverage_bees_100

write.csv(coverage_summary, "data/final_and_intermediate_outputs/Coverage_bees_100.csv")

## Check coverage for bees at n = 200----

sample_size <- 200

abundance_transposed_filtered <- abundance_transposed[rowSums(abundance_transposed) >= sample_size, ]

#slow to run
coverage_200 <- apply(
  abundance_transposed_filtered, 
  1, 
  estimate_coverage_fixed, 
  n = sample_size
)

coverage_summary_200 <- bind_rows(coverage_200, .id = "sample") %>%
  filter(m == 200) %>%   # keep only the rows corresponding to the selected sample size.
  mutate(
    region_year = sample
  ) %>%
  separate(region_year, into = c("region_100", "endYear"), sep = "_Y", remove = FALSE) %>%
  mutate(
    region_100 = gsub("R", "", region_100),
    endYear = as.integer(endYear)
  )

coverage_summary_200 %>% group_by(region_100) %>% summarise(mean = mean(SC))

write.csv(coverage_summary, "data/final_and_intermediate_outputs/Coverage_bees_200.csv")

p_coverage_bees_200 <- ggplot(coverage_summary_200, aes(x = endYear, y = SC, color = region_100)) +
  geom_point() +
  geom_line() +
  facet_wrap(. ~ region_100, nrow = 1) +
  scale_color_manual(values = region_colors) +
  ylim(0, 1) +
  labs(x = "Year", y = "Coverage (200)") +
  theme_minimal() +
  theme(legend.position = "none"
  )

p_coverage_bees_200


### Coverage for hoverflies at n = 100----

species <- read.csv("data/final_and_intermediate_outputs/database_clean_filtered.csv", stringsAsFactors = T) # 4158540 records.

species <- species %>% filter (region_100 %in% c("alpine", "mediterranean", "continental", "atlantic", "boreal")) # 4135152 records.

species <- species %>% filter (taxo_group %in% "hoverflies") #Keep only hoverflies # 640524 records.

# Count occurrences of each species per bioregion and year.
species_counts <- species %>%
  group_by(region_100, endYear, scientificName) %>%
  summarise(count = n(), .groups = "drop")

# Create a unique ID poer each region-year combination.
species_counts <- species_counts %>%
  mutate(region_year = paste0("R", region_100, "_Y", endYear))

# Transpose to wide to get abundance or occurrences.
abundance_list <- species_counts %>%
  select(region_year, scientificName, count) %>%
  pivot_wider(names_from = scientificName, values_from = count, values_fill = 0) %>%
  column_to_rownames("region_year") %>%
  as.data.frame() %>%
  split(., rownames(.)) %>%
  lapply(as.numeric)


# Combine the list into a matrix
abundance_matrix <- do.call(cbind, abundance_list)

# Transpose so rows = samples, columns = species
abundance_transposed <- t(abundance_matrix)

sample_size <- 100

abundance_transposed_filtered <- abundance_transposed[rowSums(abundance_transposed) >= sample_size, ]

#slow
coverage_100 <- apply(
  abundance_transposed_filtered, 
  1, 
  estimate_coverage_fixed, 
  n = sample_size
)

coverage_summary <- bind_rows(coverage_100, .id = "sample") %>%
  filter(m == 100) %>%   # keep only the rows corresponding to the selected sample size.
  mutate(
    region_year = sample
  ) %>%
  separate(region_year, into = c("region_100", "endYear"), sep = "_Y", remove = FALSE) %>%
  mutate(
    region_100 = gsub("R", "", region_100),
    endYear = as.integer(endYear)
  )

coverage_summary %>% group_by(region_100) %>% summarise(mean = mean(SC))

write.csv(coverage_summary, "data/final_and_intermediate_outputs/Coverage_hoverflies_100.csv")

p_coverage_hoverflies_100 <- ggplot(coverage_summary, aes(x = endYear, y = SC, color = region_100)) +
  geom_point() +
  geom_line() +
  facet_wrap(. ~ region_100, nrow = 1) +
  scale_color_manual(values = region_colors) +
  ylim(0, 1) +
  labs(x = "Year", y = "Coverage (100)") +
  theme_minimal() +
  theme(legend.position = "none")

p_coverage_hoverflies_100


### Coverage for hoverflies at n = 200----

sample_size <- 200

abundance_transposed_filtered <- abundance_transposed[rowSums(abundance_transposed) >= sample_size, ]

#slow
coverage_200 <- apply(
  abundance_transposed_filtered, 
  1, 
  estimate_coverage_fixed, 
  n = sample_size
)

coverage_summary_200 <- bind_rows(coverage_200, .id = "sample") %>%
  filter(m == 200) %>%   # keep only the rows corresponding to the selected sample size.
  mutate(
    region_year = sample
  ) %>%
  separate(region_year, into = c("region_100", "endYear"), sep = "_Y", remove = FALSE) %>%
  mutate(
    region_100 = gsub("R", "", region_100),
    endYear = as.integer(endYear)
  )

coverage_summary_200 %>% group_by(region_100) %>% summarise(mean = mean(SC))

write.csv(coverage_summary_200, "data/final_and_intermediate_outputs/Coverage_hoverflies_200.csv")

p_coverage_hoverflies_200 <- ggplot(coverage_summary_200, aes(x = endYear, y = SC, color = region_100)) +
  geom_point() +
  geom_line() +
  facet_wrap(. ~ region_100, nrow = 1) +
  scale_color_manual(values = region_colors) +
  ylim(0, 1) +
  labs(x = "Year", y = "Coverage (200)") +
  theme_minimal() +
  theme(legend.position = "none")

p_coverage_hoverflies_200
