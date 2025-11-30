############# FIRST ASSEMBLE ALL SPECIES TOGETHER
pkgs <- c("emmeans", "tidyverse", "lme4", "MuMIn") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

## Logistic regressions----

species <- read.csv("data/final_and_intermediate_outputs/database_clean_filtered.csv", stringsAsFactors = T) # 4158540 records.

data_filtered <- species %>%
  group_by(endYear, scientificName, region_50) %>%
  summarise(Freq = n(), .groups = "drop")

species_by_region <- as.data.frame(spread (data_filtered, key=scientificName, value = Freq))

# total records per year in alpine (from raw observations)
obs_by_year <- species %>%
  group_by(region_50,endYear) %>%
  summarise(total = n(), .groups = "drop")

# Identify species columns
species_cols <- setdiff(names(species_by_region), c("endYear", "region_50"))

species_matrix <- species_by_region %>%
  group_by(endYear, region_50) %>%
  mutate(total_abundance = rowSums(across(all_of(species_cols)), na.rm = TRUE)) %>%
  mutate(across(all_of(species_cols), ~ . / total_abundance)) %>%
  select(-total_abundance) %>%
  ungroup()

species_matrix$endYear <- as.factor(as.character(species_matrix$endYear))

species_matrix_alpine <- species_matrix %>% filter (region_50 == "alpine")

species_matrix_atlantic <- species_matrix %>% filter (region_50 == "atlantic")

species_matrix_boreal <- species_matrix %>% filter (region_50 == "boreal")

species_matrix_mediterranean <- species_matrix %>% filter (region_50 == "mediterranean")

species_matrix_continental <- species_matrix %>% filter (region_50 == "continental")


##Clean datasets (mutate NAs to 0s, remove species that appear fewer than 20 years)

#Mutate Nas to 0s.

species_matrix_alpine <- species_matrix_alpine %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

species_matrix_atlantic <- species_matrix_atlantic %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

species_matrix_boreal <- species_matrix_boreal %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

species_matrix_mediterranean <- species_matrix_mediterranean %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

species_matrix_continental <- species_matrix_continental %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))


##Alpine
# Identify numeric columns (species)
species_cols <- species_matrix_alpine %>%
  select(where(is.numeric)) %>%
  names()

# Count in how many rows each species appears
species_presence <- colSums(species_matrix_alpine[species_cols] > 0)

# Keep only those appearing in >= 20 years
species_to_keep <- names(species_presence[species_presence >= 20]) ## More than 20 years present.

# Rebuild the dataset with kept species
species_matrix_alpine <- species_matrix_alpine %>%
  select(endYear, region_50, all_of(species_to_keep))

# Atlantic
species_cols <- species_matrix_atlantic %>%
  select(where(is.numeric)) %>%
  names()

species_presence <- colSums(species_matrix_atlantic[species_cols] > 0)

species_to_keep <- names(species_presence[species_presence >= 20]) ## More than 20 years present.

species_matrix_atlantic <- species_matrix_atlantic %>%
  select(endYear, region_50, all_of(species_to_keep))


# Continental
species_cols <- species_matrix_continental %>%
  select(where(is.numeric)) %>%
  names()

species_presence <- colSums(species_matrix_continental[species_cols] > 0)

species_to_keep <- names(species_presence[species_presence >= 20]) ## More than 20 years present.

species_matrix_continental <- species_matrix_continental %>%
  select(endYear, region_50, all_of(species_to_keep))

# Boreal
species_cols <- species_matrix_boreal %>%
  select(where(is.numeric)) %>%
  names()

species_presence <- colSums(species_matrix_boreal[species_cols] > 0)

species_to_keep <- names(species_presence[species_presence >= 20]) ## More than 20 years present.

species_matrix_boreal <- species_matrix_boreal %>%
  select(endYear, region_50, all_of(species_to_keep))

# Mediterranean
species_cols <- species_matrix_mediterranean %>%
  select(where(is.numeric)) %>%
  names()

species_presence <- colSums(species_matrix_mediterranean[species_cols] > 0)

species_to_keep <- names(species_presence[species_presence >= 20]) ## More than 20 years present.

species_matrix_mediterranean <- species_matrix_mediterranean %>%
  select(endYear, region_50, all_of(species_to_keep))


#Re-stablish year to numeric.

species_matrix_alpine$endYear <- as.numeric(as.character(species_matrix_alpine$endYear))
species_matrix_atlantic$endYear <- as.numeric(as.character(species_matrix_atlantic$endYear))
species_matrix_boreal$endYear <- as.numeric(as.character(species_matrix_boreal$endYear))
species_matrix_mediterranean$endYear <- as.numeric(as.character(species_matrix_mediterranean$endYear))
species_matrix_continental$endYear <- as.numeric(as.character(species_matrix_continental$endYear))

##Alpine region----
## Iterate across all species.

# join totals into your species_matrix_alpine (which has one row per year)
species_matrix_alpine <- species_matrix_alpine %>%
  left_join(obs_by_year, by = c("endYear","region_50"))

Results <- data.frame(Estimate= NA, Std_error = NA, Estimate_emm = NA, 
                      LCI = NA, UCI = NA, Species = NA)

# Iterate through each response variable column and fit models
for (i in 3:(ncol(species_matrix_alpine)-1)) {  
  i
  # Fit the model for each response variable
  model <- glm(formula = species_matrix_alpine[[i]] ~ endYear, data = species_matrix_alpine, family = binomial, weights = total)
  
  # Print the summary of each model
  cat("\nSummary for model with", names(species_matrix_alpine)[i], "as response variable:\n")
  print(summary(model))
  Estimate_emm <-emmeans(model,specs="endYear")
  temp <- data.frame(Estimate = summary(model)$coefficients[2,1], Std_error = summary(model)$coefficients[2,2], Estimate_emm = summary(Estimate_emm)[[2]],LCI = summary(Estimate_emm)[[5]], UCI = summary(Estimate_emm)[[6]], Species = names(species_matrix_alpine)[i])
  Results <- rbind(Results, temp)
  
}

Results_alpine <- Results[-1,]

Results_alpine <- Results_alpine[order(Results_alpine$Estimate, decreasing = TRUE), ]

Results_alpine$region_50 <- "alpine"


# Atlantic

Results <- data.frame(Estimate= NA, Std_error = NA, Estimate_emm = NA, 
                      LCI = NA, UCI = NA, Species = NA)

species_matrix_atlantic <- species_matrix_atlantic %>%
  left_join(obs_by_year, by = c("endYear","region_50"))

# Iterate through each response variable column and fit models
for (i in 3:(ncol(species_matrix_atlantic)-1)) {  
  i
  # Fit the model for each response variable
  model <- glm(formula = species_matrix_atlantic[[i]] ~ endYear, data = species_matrix_atlantic, family = binomial, weights = total)
  
  # Print the summary of each model
  cat("\nSummary for model with", names(species_matrix_atlantic)[i], "as response variable:\n")
  print(summary(model))
  Estimate_emm <-emmeans(model,specs="endYear")
  temp <- data.frame(Estimate = summary(model)$coefficients[2,1], Std_error = summary(model)$coefficients[2,2], Estimate_emm = summary(Estimate_emm)[[2]],LCI = summary(Estimate_emm)[[5]], UCI = summary(Estimate_emm)[[6]], Species = names(species_matrix_atlantic)[i])
  Results <- rbind(Results, temp)
  
}

Results_atlantic <- Results[-1,]

Results_atlantic <- Results_atlantic[order(Results_atlantic$Estimate, decreasing = TRUE), ]

Results_atlantic$region_50 <- "atlantic"


# Boreal

Results <- data.frame(Estimate= NA, Std_error = NA, Estimate_emm = NA, 
                      LCI = NA, UCI = NA, Species = NA)

species_matrix_boreal <- species_matrix_boreal %>%
  left_join(obs_by_year, by = c("endYear","region_50"))

# Iterate through each response variable column and fit models
for (i in 3:(ncol(species_matrix_boreal)-1)) {  
  i
  # Fit the model for each response variable
  model <- glm(formula = species_matrix_boreal[[i]] ~ endYear, data = species_matrix_boreal, family = binomial, weights = total)
  
  # Print the summary of each model
  cat("\nSummary for model with", names(species_matrix_boreal)[i], "as response variable:\n")
  print(summary(model))
  Estimate_emm <-emmeans(model,specs="endYear")
  temp <- data.frame(Estimate = summary(model)$coefficients[2,1], Std_error = summary(model)$coefficients[2,2], Estimate_emm = summary(Estimate_emm)[[2]],LCI = summary(Estimate_emm)[[5]], UCI = summary(Estimate_emm)[[6]], Species = names(species_matrix_boreal)[i])
  Results <- rbind(Results, temp)
}

Results_boreal <- Results[-1,]

Results_boreal <- Results_boreal[order(Results_boreal$Estimate, decreasing = TRUE), ]

Results_boreal$region_50 <- "boreal"


# Mediterranean

Results <- data.frame(Estimate= NA, Std_error = NA, Estimate_emm = NA, 
                      LCI = NA, UCI = NA, Species = NA)

species_matrix_mediterranean <- species_matrix_mediterranean %>%
  left_join(obs_by_year, by = c("endYear","region_50"))

# Iterate through each response variable column and fit models
for (i in 3:(ncol(species_matrix_mediterranean)-1)) {  
  i
  # Fit the model for each response variable
  model <- glm(formula = species_matrix_mediterranean[[i]] ~ endYear, data = species_matrix_mediterranean, family = binomial, weights = total)
  
  # Print the summary of each model
  cat("\nSummary for model with", names(species_matrix_mediterranean)[i], "as response variable:\n")
  print(summary(model))
  Estimate_emm <-emmeans(model,specs="endYear")
  temp <- data.frame(Estimate = summary(model)$coefficients[2,1], Std_error = summary(model)$coefficients[2,2], Estimate_emm = summary(Estimate_emm)[[2]],LCI = summary(Estimate_emm)[[5]], UCI = summary(Estimate_emm)[[6]], Species = names(species_matrix_mediterranean)[i])
  Results <- rbind(Results, temp)
}

Results_mediterranean <- Results[-1,]

Results_mediterranean <- Results_mediterranean[order(Results_mediterranean$Estimate, decreasing = TRUE), ]

Results_mediterranean$region_50 <- "mediterranean"


# Continental

Results <- data.frame(Estimate= NA, Std_error = NA, Estimate_emm = NA, 
                      LCI = NA, UCI = NA, Species = NA)

species_matrix_continental <- species_matrix_continental %>%
  left_join(obs_by_year, by = c("endYear","region_50"))

# Iterate through each response variable column and fit models
for (i in 3:(ncol(species_matrix_continental)-1)) {  
  i
  # Fit the model for each response variable
  model <- glm(formula = species_matrix_continental[[i]] ~ endYear, data = species_matrix_continental, family = binomial, weights = total)
  
  # Print the summary of each model
  cat("\nSummary for model with", names(species_matrix_continental)[i], "as response variable:\n")
  print(summary(model))
  Estimate_emm <-emmeans(model,specs="endYear")
  temp <- data.frame(Estimate = summary(model)$coefficients[2,1], Std_error = summary(model)$coefficients[2,2], Estimate_emm = summary(Estimate_emm)[[2]],LCI = summary(Estimate_emm)[[5]], UCI = summary(Estimate_emm)[[6]], Species = names(species_matrix_continental)[i])
  Results <- rbind(Results, temp)
}

Results_continental <- Results[-1,]

Results_continental <- Results_continental[order(Results_continental$Estimate, decreasing = TRUE), ]

Results_continental$region_50 <- "continental"


Logistic_models <- rbind (Results_alpine, Results_atlantic, 
                          Results_boreal, Results_mediterranean, Results_continental)

#Load_trends

Trends <- read.csv("all_trends.csv")

colnames(Trends)[1]<- "Species"

Trends_1921 <- Trends %>% filter (baseline == "1921")

Trends_1921 <- Trends_1921 %>% filter (trend >=-1, trend <=1)

# Merge occupancy trends with trends estimated from logistic models by species and bioregion.

Trends_logistic_occupancy_1921 <- left_join (Logistic_models, Trends_1921, by = c("Species", "region_50"))

#write.csv(Trends_logistic_occupancy_1921, "Trends_logistic_occupancy_1921.csv")

#Trends_logistic_occupancy_1921 <- read.csv("Trends_logistic_occupancy_1921.csv")


pearson <- cor.test(scale(Trends_logistic_occupancy_1921$trend), scale(Trends_logistic_occupancy_1921$Estimate),
    use = "complete.obs", method = "pearson") 

pearson

spearman <- cor.test(scale(Trends_logistic_occupancy_1921$trend), scale(Trends_logistic_occupancy_1921$Estimate),
         use = "complete.obs", method = "spearman")

spearman


p1 <- ggplot(Trends_logistic_occupancy_1921, aes(x = scale(trend), y = scale(Estimate))) +
  geom_point(alpha=0.22) +
  xlim(values = c(-5,5))+
  geom_smooth(method="lm", fill = "lightblue") +
  labs(
    x = "Standardized raw occupancy trend",
    y = "Standardized raw logistic trend",
  )+
  theme_minimal() +
  theme(
    axis.text.x = element_text(),
    plot.title = element_text(hjust = 0.5, margin = margin(b = 25)),
    plot.title.position = "panel"
  )

p1

## Confusion matrix.

Trends_logistic_occupancy_1921 <- Trends_logistic_occupancy_1921 %>%
  mutate(
    trend_class = case_when(
      Estimate > 0 ~ "log. positive",
      Estimate < 0 ~ "log. negative",
      TRUE ~ "occ. nonsignificant"
      ),
    model_class = case_when(
      trend > 0  ~ "occ. positive",
      trend < 0  ~ "occ. negative",
      TRUE ~ "occ. nonsignificant"
    )
  )

confusion_full <- table(Trends_logistic_occupancy_1921$trend_class,
                        Trends_logistic_occupancy_1921$model_class)

confusion_full
