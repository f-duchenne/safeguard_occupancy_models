
## Script to estimate changes in pollinator richness in the last century across European regions. 

#Install and load required packages

#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("vegan", "dplyr","ggplot2","mgcv", "gratia","ggeffects", "iNEXT") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

#defining working folder:
setwd("~/Desktop/Juan de la Cierva/Papers 2024/Preparation/STEP/Script and Data") ## Nacho Check! Change to the folder below?
#project_folder="C:/Users/Duchenne/Documents/safeguard/" ## Nacho Check!!

select <- dplyr::select

species <- read.csv("database_clean_filtered.csv", stringsAsFactors = T) ## 4158540 records.

species_summary <- species %>%
  group_by(taxo_group) %>%
  distinct(scientificName) %>%
  summarise(n_species = n())

species_summary ## 1888 bee species and 784 hoverfly species (2672)

species <- species %>% filter (region_100 %in% c("alpine", "mediterranean", "continental", "atlantic", "boreal"))

species_summary <- species %>%
  group_by(taxo_group) %>%
  distinct(scientificName) %>%
  summarise(n_species = n())

species_summary ## 1832 bee species and 774 hoverfly species (2606)

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

nrow(abundance_transposed) # 500 region*year combinations (5 regions, 100 years each).


### Check coverage at n = 100, and compute rarefied richness
## Sample size = 100
# Bees

sample_size <- 100

abundance_transposed_filtered <- abundance_transposed[rowSums(abundance_transposed) >= sample_size, ]

nrow(abundance_transposed)-nrow(abundance_transposed_filtered) ##The rows (year*region combination) we lose.


## First Check coverage:

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

### Final GAM models for richness----

# ---- Fit the GAM (you already did) ----
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


# ---- Get fitted values for plotting ----

#pred_bees <- gratia::smooth_estimates(gam_model_bees, term = "endYear")

pred_full <- ggpredict(gam_model_bees, terms = c("endYear", "region_100"))

colnames(pred_full) <- c("endYear",".estimate","std_error","conf.low","conf.high","region_100")

deriv_bees$endYear <- round (deriv_bees$endYear)

# Merge derivative significance with fitted values
plot_data <- left_join(pred_full, deriv_bees, 
                       by = c("endYear", "region_100"))

bee_plot_100 <- ggplot(plot_data, aes(x = endYear, y = .estimate)) +
  
  # Raw points
  geom_point(data = Dataset_Richness_bees,
             aes(x = endYear, y = Estimator, color = region_100, fill = region_100),
             alpha = 0.3, size = 1.5, inherit.aes = FALSE) +
  
  # Increasing / Decreasing = solid
  geom_line(
    data = subset(plot_data, sig != "Stable"),
    aes(color = region_100, group = interaction(region_100, sig)),
    linetype = "solid",
    size = 1.2
  ) +
  
  # Stable = dashed
  geom_line(
    data = subset(plot_data, sig == "Stable"),
    aes(color = region_100, group = interaction(region_100, sig)),
    linetype = "dashed",
    size = 1.2
  ) +
  
  # Confidence ribbon
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = region_100),
    alpha = 0.15, color = NA
  ) +
  coord_cartesian(ylim = c(20, NA))+
  scale_color_manual(values = region_colors) +
  scale_fill_manual(values = region_colors) +
  
  facet_wrap(~str_to_title(region_100), ncol = 1, scales = "free_y") +
  labs(x = "Year", y = "Predicted bee richness", color = "Biogeoregion") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

bee_plot_100


## sample size = 200; bees.

sample_size <- 200

abundance_transposed_filtered <- abundance_transposed[rowSums(abundance_transposed) >= sample_size, ]

nrow(abundance_transposed)-nrow(abundance_transposed_filtered) ##The rows (year*region combination) we lose.


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

p_coverage_bees_200 <- ggplot(coverage_summary_200, aes(x = endYear, y = SC, color = region_100)) +
  geom_point() +
  geom_line() +
  facet_wrap(. ~ region_100, nrow = 1) +
  scale_color_manual(values = region_colors) +
  ylim(0, 1) +
  labs(x = "Year", y = "Coverage (200)") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_blank()   # <-- hides facet labels
  )

p_coverage_bees_200


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

### Final GAM models for richness----

# ---- Fit the GAM (you already did) ----
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


# ---- Get fitted values for plotting ----

#pred_bees <- gratia::smooth_estimates(gam_model_bees, term = "endYear")

pred_full <- ggpredict(gam_model_bees, terms = c("endYear", "region_100"))

colnames(pred_full) <- c("endYear",".estimate","std_error","conf.low","conf.high","region_100")

deriv_bees$endYear <- round (deriv_bees$endYear)

# Merge derivative significance with fitted values
plot_data <- left_join(pred_full, deriv_bees, 
                       by = c("endYear", "region_100"))

bee_plot_200 <- ggplot(plot_data, aes(x = endYear, y = .estimate)) +
  
  # Raw points
  geom_point(data = Dataset_Richness_bees,
             aes(x = endYear, y = Estimator, color = region_100, fill = region_100),
             alpha = 0.3, size = 1.5, inherit.aes = FALSE) +
  
  # Increasing / Decreasing = solid
  geom_line(
    data = subset(plot_data, sig != "Stable"),
    aes(color = region_100, group = interaction(region_100, sig)),
    linetype = "solid",
    size = 1.2
  ) +
  
  # Stable = dashed
  geom_line(
    data = subset(plot_data, sig == "Stable"),
    aes(color = region_100, group = interaction(region_100, sig)),
    linetype = "dashed",
    size = 1.2
  ) +
  
  # Confidence ribbon
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = region_100),
    alpha = 0.15, color = NA
  ) +
  coord_cartesian(ylim = c(20, NA))+
  scale_color_manual(values = region_colors) +
  scale_fill_manual(values = region_colors) +
  
  facet_wrap(~str_to_title(region_100), ncol = 1, scales = "free_y") +
  labs(x = "Year", y = "Predicted bee richness", color = "Biogeoregion") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

bee_plot_200


##The same with hoverflies----

species <- read.csv("database_clean_filtered.csv", stringsAsFactors = T) # 4158540 records.

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


## Sample size = 100, hoverflies

# Compute rarefied richness
sample_size <- 100

abundance_transposed_filtered <- abundance_transposed[rowSums(abundance_transposed) >= sample_size, ]

nrow(abundance_transposed)-nrow(abundance_transposed_filtered) ##The rows (year*region combination) we lose. We lose 202 region*year combinations.


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

p_coverage_hoverflies_100 <- ggplot(coverage_summary, aes(x = endYear, y = SC, color = region_100)) +
  geom_point() +
  geom_line() +
  facet_wrap(. ~ region_100, nrow = 1) +
  scale_color_manual(values = region_colors) +
  ylim(0, 1) +
  labs(x = "Year", y = "Coverage (100)") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_blank()   # <-- hides facet labels
  )

p_coverage_hoverflies_100

# 3. Compute rarefied richness (e.g., to minimum sample size)

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


# ---- Fit the GAM (you already did) ----

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


# ---- Get fitted values for plotting ----

#pred_bees <- gratia::smooth_estimates(gam_model_bees, term = "endYear")

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

hoverfly_plot_100 <- ggplot(plot_data, aes(x = endYear, y = .estimate)) +
  # Raw points
  geom_point(data = Dataset_Richness_hover,
             aes(x = endYear, y = Estimator, color = region_100, fill = region_100),
             alpha = 0.3, size = 1.5, inherit.aes = FALSE) +
  # Increasing / Decreasing = solid
  geom_line(
    data = subset(plot_data, sig != "Stable"),
    aes(color = region_100, group = interaction(region_100, sig)),
    linetype = "solid",
    size = 1.2
  ) +
  
  # Stable = dashed
  geom_line(
    data = subset(plot_data, sig == "Stable"),
    aes(color = region_100, group = interaction(region_100, sig)),
    linetype = "dashed",
    size = 1.2
  ) +
  
  # Confidence ribbon
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = region_100),
    alpha = 0.15, color = NA
  ) +
  scale_color_manual(values = region_colors) +
  scale_fill_manual(values = region_colors) +
  facet_wrap(~str_to_title(region_100), ncol = 1, scales = "free_y") +
  labs(x = "Year", y = "Predicted hoverfly richness", color = "Biogeoregion") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

hoverfly_plot_100


## Now using 200 redords, hoverflies----

sample_size <- 200

abundance_transposed_filtered <- abundance_transposed[rowSums(abundance_transposed) >= sample_size, ]

nrow(abundance_transposed)-nrow(abundance_transposed_filtered) ##The rows (year*region combination) we lose.


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

p_coverage_hoverflies_200 <- ggplot(coverage_summary_200, aes(x = endYear, y = SC, color = region_100)) +
  geom_point() +
  geom_line() +
  facet_wrap(. ~ region_100, nrow = 1) +
  scale_color_manual(values = region_colors) +
  ylim(0, 1) +
  labs(x = "Year", y = "Coverage (200)") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_blank()   # <-- hides facet labels
  )


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

# ---- Fit the GAM (you already did) ----
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

# ---- Get fitted values for plotting ----

#pred_bees <- gratia::smooth_estimates(gam_model_bees, term = "endYear")

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


hoverfly_plot_200 <- ggplot(plot_data, aes(x = endYear, y = .estimate)) +
  # Raw points
  geom_point(data = Dataset_Richness_hover,
             aes(x = endYear, y = Estimator, color = region_100, fill = region_100),
             alpha = 0.3, size = 1.5, inherit.aes = FALSE) +
  # Increasing / Decreasing = solid
  geom_line(
    data = subset(plot_data, sig != "Stable"),
    aes(color = region_100, group = interaction(region_100, sig)),
    linetype = "solid",
    size = 1.2
  ) +

  # Stable = dashed
  geom_line(
    data = subset(plot_data, sig == "Stable"),
    aes(color = region_100, group = interaction(region_100, sig)),
    linetype = "dashed",
    size = 1.2
  ) +
  
  # Confidence ribbon
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = region_100),
    alpha = 0.15, color = NA
  ) +
  scale_color_manual(values = region_colors) +
  scale_fill_manual(values = region_colors) +
  facet_wrap(~str_to_title(region_100), ncol = 1, scales = "free_y") +
  labs(x = "Year", y = "Predicted hoverfly richness", color = "Biogeoregion") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

hoverfly_plot_200


bee_plot_100 + hoverfly_plot_100

bee_plot_200 + hoverfly_plot_200


## Combined Plot for coverage----

p_coverage_bees_100 + 
  p_coverage_bees_200 + 
  p_coverage_hoverflies_100 +
  p_coverage_hoverflies_200 + plot_layout (nrow = 4)


