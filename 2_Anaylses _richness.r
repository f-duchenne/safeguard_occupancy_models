## Script to estimate changes in pollinator richness in the last century across European regions. 

#Install and load required packages
rm(list=ls())
pkgs <- c("vegan", "dplyr","ggplot2","mgcv", "gratia","ggeffects", "iNEXT", "textshape", "tidyr") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

select <- dplyr::select

#Load dataset
species <- read.csv("data/final_and_intermediate_outputs/database_clean_filtered.csv", stringsAsFactors = T) ## 4158540 records.

species_summary <- species %>%
  group_by(taxo_group) %>%
  distinct(scientificName) %>%
  summarise(n_species = n())

species_summary ## 1888 bee species and 784 hoverfly species (2672)

species <- species %>% filter (region_100 %in% c("alpine", "mediterranean", "continental", "atlantic", "boreal"))

# Bees----
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


## Bees: n = 100----

sample_size <- 100

abundance_transposed_filtered <- abundance_transposed[rowSums(abundance_transposed) >= sample_size, ]

nrow(abundance_transposed)-nrow(abundance_transposed_filtered) ##The rows (year*region combination) we lose.

richness <- rarefy(abundance_transposed_filtered, sample = sample_size)

richness_df <- data.frame(
  sample = rownames(abundance_transposed_filtered),
  richness = richness
)

summary_richness <- richness_df %>%
  mutate(region_year = as.character(sample)) %>%
  separate(region_year, into = c("region_100", "endYear"), sep = "_Y", remove = FALSE) %>%
  mutate(region_100 = gsub("R", "", region_100),
         endYear = as.integer(endYear))

# Estimate number of counts per region and year.
record_count <- species %>%
  group_by(region_100, endYear) %>%
  summarise(record_count = n(), .groups = "drop")

# Join with results from iNEXT.
Dataset <- left_join(summary_richness, record_count, by = c("region_100", "endYear"))

Dataset <- Dataset %>% filter (region_100 %in% c("alpine", "mediterranean", "continental", "atlantic", "boreal"))

Dataset$region_100 <- as.factor(Dataset$region_100)

## Fit GAM Model
## Model overall trend in each bioregion:

Dataset_Richness_bees <- Dataset

colnames(Dataset_Richness_bees)[2] <- "Estimator"

#write.csv(Dataset_Richness_bees, "data/final_and_intermediate_outputs/Dataset_Richness_bees_100.csv")

### GAM models for richness

gam_model_bees <- gam(
  Estimator ~ s(endYear, by = region_100, k = 3) + region_100,
  family = gaussian(),
  weights = record_count,
  data = Dataset_Richness_bees,
  method = "REML"
)

# Derivatives for each smooth (region)
deriv_bees <- derivatives(
  gam_model_bees,
  select = "endYear",        # new syntax replacing 'term'
  interval = "confidence",
  n = 100,
  level = 0.99,
  partial_match = TRUE       # allows matching e.g. s(endYear):region_100X
)

# Add significance classification
deriv_bees <- deriv_bees %>%
  mutate(
    sig = case_when(
      .lower_ci > 0 ~ "Increasing",
      .upper_ci < 0 ~ "Decreasing",
      TRUE ~ "Stable"
    )
  )

deriv_bees <- deriv_bees %>% dplyr::select (.smooth,.derivative,.lower_ci, .upper_ci,endYear,region_100, sig)


#Get fitted values for plotting

pred_full <- ggpredict(gam_model_bees, terms = c("endYear", "region_100"))

colnames(pred_full) <- c("endYear",".estimate","std_error","conf.low","conf.high","region_100")

deriv_bees$endYear <- round (deriv_bees$endYear)

# Merge derivative significance with fitted values
plot_data <- left_join(pred_full, deriv_bees, 
                       by = c("endYear", "region_100"))

#write.csv(plot_data, "data/final_and_intermediate_outputs/Plot_data_bees_100.csv")

## Bees: n = 200----

sample_size <- 200

abundance_transposed_filtered <- abundance_transposed[rowSums(abundance_transposed) >= sample_size, ]

nrow(abundance_transposed)-nrow(abundance_transposed_filtered) ##The rows (year*region combination) we lose.

richness <- rarefy(abundance_transposed_filtered, sample = sample_size)

richness_df <- data.frame(
  sample = rownames(abundance_transposed_filtered),
  richness = richness
)

summary_richness <- richness_df %>%
  mutate(region_year = as.character(sample)) %>%
  separate(region_year, into = c("region_100", "endYear"), sep = "_Y", remove = FALSE) %>%
  mutate(region_100 = gsub("R", "", region_100),
         endYear = as.integer(endYear))

# Estimate number of counts per region and year.
record_count <- species %>%
  group_by(region_100, endYear) %>%
  summarise(record_count = n(), .groups = "drop")

# Join with results from iNEXT.
Dataset <- left_join(summary_richness, record_count, by = c("region_100", "endYear"))

Dataset <- Dataset %>% filter (region_100 %in% c("alpine", "mediterranean", "continental", "atlantic", "boreal"))

Dataset$region_100 <- as.factor(Dataset$region_100)

## Fit GAM Model
## Model overall trend in each bioregion:

Dataset_Richness_bees <- Dataset

colnames(Dataset_Richness_bees)[2] <- "Estimator"

#write.csv(Dataset_Richness_bees, "data/final_and_intermediate_outputs/Dataset_Richness_bees_200.csv")

### GAM models for richness

#Fit the GAM
gam_model_bees <- gam(
  Estimator ~ s(endYear, by = region_100, k = 3) + region_100,
  family = gaussian(),
  weights = record_count,
  data = Dataset_Richness_bees,
  method = "REML"
)

# Derivatives for each smooth (region)
deriv_bees <- derivatives(
  gam_model_bees,
  select = "endYear",        # new syntax replacing 'term'
  interval = "confidence",
  n = 100,
  level = 0.99,
  partial_match = TRUE       # allows matching e.g. s(endYear):region_100X
)

# Add significance classification
deriv_bees <- deriv_bees %>%
  mutate(
    sig = case_when(
      .lower_ci > 0 ~ "Increasing",
      .upper_ci < 0 ~ "Decreasing",
      TRUE ~ "Stable"
    )
  )

deriv_bees <- deriv_bees %>% dplyr::select (.smooth,.derivative,.lower_ci, .upper_ci,endYear,region_100, sig)


#Get fitted values for plotting

pred_full <- ggpredict(gam_model_bees, terms = c("endYear", "region_100"))

colnames(pred_full) <- c("endYear",".estimate","std_error","conf.low","conf.high","region_100")

deriv_bees$endYear <- round (deriv_bees$endYear)

# Merge derivative significance with fitted values
plot_data_200 <- left_join(pred_full, deriv_bees, 
                       by = c("endYear", "region_100"))

#write.csv(plot_data_200, "data/final_and_intermediate_outputs/Plot_data_bees_200.csv")


# Hoverflies----

species <- read.csv("data/final_and_intermediate_outputs/database_clean_filtered.csv", stringsAsFactors = T) # 4158540 records.

species <- species %>% filter (region_100 %in% c("alpine", "mediterranean", "continental", "atlantic", "boreal")) # 4135152 records.

species <- species %>% filter (taxo_group %in% "hoverflies") #Keep only bees. # 640524 records.

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

nrow(abundance_transposed)# Number of region*years combinations.


## Hoverflies: n = 100----

# Compute rarefied richness
sample_size <- 100

abundance_transposed_filtered <- abundance_transposed[rowSums(abundance_transposed) >= sample_size, ]

nrow(abundance_transposed)-nrow(abundance_transposed_filtered) ##The rows (year*region combination) we lose. We lose 202 region*year combinations.

richness <- rarefy(abundance_transposed_filtered, sample = sample_size)

richness_df <- data.frame(
  sample = rownames(abundance_transposed_filtered),
  richness = richness
)

summary_richness <- richness_df %>%
  mutate(region_year = as.character(sample)) %>%
  separate(region_year, into = c("region_100", "endYear"), sep = "_Y", remove = FALSE) %>%
  mutate(region_100 = gsub("R", "", region_100),
         endYear = as.integer(endYear))

# Estimate number of counts per region and year.
record_count <- species %>%
  group_by(region_100, endYear) %>%
  summarise(record_count = n(), .groups = "drop")

# Join with results from iNEXT.
Dataset_hover <- left_join(summary_richness, record_count, by = c("region_100", "endYear"))

Dataset_hover <- Dataset_hover %>% filter (region_100 %in% c("alpine", "mediterranean", "continental", "atlantic", "boreal"))

Dataset_hover$region_100 <- as.factor(Dataset_hover$region_100)

## Model

## Model overall trend in each bioregion:

Dataset_Richness_hover <- Dataset_hover

colnames(Dataset_Richness_hover)[2] <- "Estimator"

Dataset_Richness_hover <- Dataset_Richness_hover %>%
  filter(!(region_100 == "continental" & endYear == 1938)) ## Remove temporal outlier.

#write.csv(Dataset_Richness_hover, "data/final_and_intermediate_outputs/Dataset_Richness_hover_100.csv")

#Fit the GAM
gam_model_hover <- gam(
  Estimator ~ s(endYear, by = region_100, k = 3) + region_100,
  family = gaussian(),
  weights = record_count,
  data = Dataset_Richness_hover,
  method = "REML"
)

# Derivatives for each smooth (region)
deriv_hover <- derivatives(
  gam_model_hover,
  select = "endYear",        # new syntax replacing 'term'
  interval = "confidence",
  n = 100,
  level = 0.99,
  partial_match = TRUE       # allows matching e.g. s(endYear):region_100X
)

# Add significance classification
deriv_hover <- deriv_hover %>%
  mutate(
    sig = case_when(
      .lower_ci > 0 ~ "Increasing",
      .upper_ci < 0 ~ "Decreasing",
      TRUE ~ "Stable"
    )
  )

deriv_hover <- deriv_hover %>% dplyr::select (.smooth,.derivative,.lower_ci, .upper_ci,endYear,region_100, sig)

#Get fitted values for plotting

pred_full <- ggpredict(gam_model_hover, terms = c("endYear", "region_100"))

colnames(pred_full) <- c("endYear",".estimate","std_error","conf.low","conf.high","region_100")

##Truncating years without data.
year_ranges <- Dataset_Richness_hover %>%
  group_by(region_100) %>%
  dplyr::summarise(min_year = min(endYear), max_year = max(endYear))

# Join to get year limits
pred_trimmed <- pred_full %>%
  dplyr::left_join(year_ranges, by = "region_100") %>%
  dplyr::filter(endYear >= min_year & endYear <= max_year)

deriv_hover$endYear <- round (deriv_hover$endYear)

# Merge derivative significance with fitted values
plot_data <- left_join(pred_trimmed, deriv_hover, 
                       by = c("endYear", "region_100"))

#write.csv(plot_data, "data/final_and_intermediate_outputs/Plot_data_hover_100.csv")


## Hoverflies: n=200 ----

sample_size <- 200

abundance_transposed_filtered <- abundance_transposed[rowSums(abundance_transposed) >= sample_size, ]

nrow(abundance_transposed)-nrow(abundance_transposed_filtered) ##The rows (year*region combination) we lose.

richness <- rarefy(abundance_transposed_filtered, sample = sample_size)

richness_df <- data.frame(
  sample = rownames(abundance_transposed_filtered),
  richness = richness
)

summary_richness <- richness_df %>%
  mutate(region_year = as.character(sample)) %>%
  separate(region_year, into = c("region_100", "endYear"), sep = "_Y", remove = FALSE) %>%
  mutate(region_100 = gsub("R", "", region_100),
         endYear = as.integer(endYear))

# Estimate number of counts per region and year.
record_count <- species %>%
  group_by(region_100, endYear) %>%
  summarise(record_count = n(), .groups = "drop")

# Join with results from iNEXT.
Dataset_hover <- left_join(summary_richness, record_count, by = c("region_100", "endYear"))

Dataset_hover <- Dataset_hover %>% filter (region_100 %in% c("alpine", "mediterranean", "continental", "atlantic", "boreal"))

Dataset_hover$region_100 <- as.factor(Dataset_hover$region_100)

## Model

## Model overall trend in each bioregion:

Dataset_Richness_hover <- Dataset_hover

colnames(Dataset_Richness_hover)[2] <- "Estimator"

Dataset_Richness_hover <- Dataset_Richness_hover %>%
  filter(!(region_100 == "continental" & endYear == 1938)) ## Remove temporal outlier.

#write.csv(Dataset_Richness_hover, "data/final_and_intermediate_outputs/Dataset_Richness_hover_200.csv")

#Fit the GAM
gam_model_hover <- gam(
  Estimator ~ s(endYear, by = region_100, k = 3) + region_100,
  family = gaussian(),
  weights = record_count,
  data = Dataset_Richness_hover,
  method = "REML"
)

# Derivatives for each smooth (region)
deriv_hover <- derivatives(
  gam_model_hover,
  select = "endYear",        # new syntax replacing 'term'
  interval = "confidence",
  n = 100,
  level = 0.99,
  partial_match = TRUE       # allows matching e.g. s(endYear):region_100X
)

# Add significance classification
deriv_hover <- deriv_hover %>%
  mutate(
    sig = case_when(
      .lower_ci > 0 ~ "Increasing",
      .upper_ci < 0 ~ "Decreasing",
      TRUE ~ "Stable"
    )
  )

deriv_hover <- deriv_hover %>% dplyr::select (.smooth,.derivative,.lower_ci, .upper_ci,endYear,region_100, sig)

#Get fitted values for plotting

pred_full <- ggpredict(gam_model_hover, terms = c("endYear", "region_100"))

colnames(pred_full) <- c("endYear",".estimate","std_error","conf.low","conf.high","region_100")

##Truncating years without data.
year_ranges <- Dataset_Richness_hover %>%
  group_by(region_100) %>%
  dplyr::summarise(min_year = min(endYear), max_year = max(endYear))

# Join to get year limits
pred_trimmed <- pred_full %>%
  dplyr::left_join(year_ranges, by = "region_100") %>%
  dplyr::filter(endYear >= min_year & endYear <= max_year)

deriv_hover$endYear <- round (deriv_hover$endYear)

# Merge derivative significance with fitted values
plot_data <- left_join(pred_trimmed, deriv_hover, 
                       by = c("endYear", "region_100"))

#write.csv(plot_data, "data/final_and_intermediate_outputs/Plot_data_hover_200.csv")
