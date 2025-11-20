############# FIRST ASSEMBLE ALL SPECIES TOGETHER
pkgs <- c("patchwork", "tidyverse", "brms", "stringr") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

### Fig_4.
#### Traits modelled with brms----
###############warning#############
#The code used traits_models_brms_phylo.RData
load("data/final_and_intermediate_outputs/traits_models_brms.RData") 
model_bees=models[[1]]
model_hovs=models[[2]]

## Plots for bees.
# Create marginal effects for each predictor
me_ITD <- conditional_effects(model_bees, effects = "ITD_F_mm", method = "posterior_epred")
me_STI <- conditional_effects(model_bees, effects = "STI_Species_temperature_index", method = "posterior_epred")
me_Soc <- conditional_effects(model_bees, effects = "Sociality_simplified", method = "posterior_epred")
me_LDB <- conditional_effects(model_bees, effects = "Larval_diet_breadth", method = "posterior_epred")

pred_ITD <- as.data.frame(me_ITD[[1]])
pred_STI <- as.data.frame(me_STI[[1]])
pred_Soc <- as.data.frame(me_Soc[[1]])
pred_LDB <- as.data.frame(me_LDB[[1]])

#Load traits to plot raw points.
Traits <- read.csv("data/final_and_intermediate_outputs/traits_table.csv")
Trends <- read.csv("data/final_and_intermediate_outputs/all_trends.csv")
colnames(Trends)[which(names(Trends) == "species")] <- "Species"
Trends <- Trends %>% filter (baseline == "1921")

Bee_Traits <- Traits %>% filter (taxo_group == "bees", Species %in% Trends$Species)
Bee_Traits <- Bee_Traits %>% left_join(Trends, by = "Species")

Hoverfly_Traits <- Traits %>% filter (taxo_group == "hoverflies", Species %in% Trends$Species)
Hoverfly_Traits <- Hoverfly_Traits %>% left_join(Trends, by = "Species")

plot_ITD_bee <- ggplot(pred_ITD, aes(x = effect1__, y = estimate__)) +
  geom_point(data = Bee_Traits %>% filter(!is.na(ITD_F_mm), !is.na(trend)),
             aes(x = ITD_F_mm, y = trend),
             alpha = 0.1, color = "gray50", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", size = 0.5) +
  geom_line(color = "#2ca02c", size = 1.5, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "#2ca02c", alpha = 0.2) +
  ylim(-0.1, 0.1) +
  labs(title = "Body Size", x = "Female ITD (mm)", y = "Predicted Trend") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(),
    plot.title = element_text(hjust = 0.5, margin = margin(b = 25)),
    plot.title.position = "panel"
  )

plot_STI_bee <- ggplot(pred_STI, aes(x = effect1__, y = estimate__)) +
  geom_point(data = Bee_Traits %>% filter(!is.na(STI_Species_temperature_index), !is.na(trend)),
             aes(x = STI_Species_temperature_index, y = trend),
             alpha = 0.1, color = "gray50", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", size = 0.5) +
  geom_line(color = "#2ca02c", size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "#2ca02c", alpha = 0.2) +
  ylim(-0.1, 0.1) +
  labs(title = "Temperature", x = "STI (ºC)", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(),
        plot.title = element_text(hjust = 0.5, margin = margin(b = 25)),
        plot.title.position = "panel")


plot_Soc_bee <- ggplot(pred_Soc, aes(x = effect1__, y = estimate__)) +
  geom_jitter(data = Bee_Traits %>% filter(!is.na(Sociality_simplified), !is.na(trend)),
              aes(x = Sociality_simplified, y = trend),
              width = 0.2, alpha = 0.1, color = "gray50", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", size = 0.5) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__), width = 0.2, size=1, color = "#1f77b4") +
  geom_point(size = 4, color = "#1f77b4") +
  geom_point(size = 2, color = "white") + geom_point(size = 1, color = "#1f77b4") +
  labs(title = "Life History", x = "Sociality", y = "") +
  ylim(-0.08, 0.08) +
  scale_x_discrete(labels = function(x) stringr::str_to_title(x))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, margin = margin(b = 25)),
        plot.title.position = "panel")


plot_LDB_bee <- ggplot(pred_LDB, aes(x = effect1__, y = estimate__)) +
  geom_jitter(data = Bee_Traits %>% filter(Larval_diet_breadth != "na", !is.na(trend)),
              aes(x = Larval_diet_breadth, y = trend),
              width = 0.2, alpha = 0.1, color = "gray50", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", size = 0.5) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__), width = 0.2, size=1, color = "#1f77b4") +
  geom_point(size = 4, color = "#1f77b4") +
  geom_point(size = 2, color = "white") + geom_point(size = 1, color = "#1f77b4") +
  labs(title = "Diet", x = "Larval diet breadth", y = "") +
  ylim(-0.08, 0.08) +
  scale_x_discrete(labels = function(x) stringr::str_to_title(x)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, margin = margin(b = 25)),
        plot.title.position = "panel")


##Plots for hoverflies

# Load model
model_hovs <- models[[2]]

# Generate predictions for each variable ---

# For continuous predictors
pred_BS  <- conditional_effects(model_hovs, effects = "Adult_Body_size_num", method = "posterior_epred")
pred_STI <- conditional_effects(model_hovs, effects = "STI_Species_temperature_index", method = "posterior_epred")

# For categorical predictors
pred_FH  <- conditional_effects(model_hovs, effects = "Flight_height", method = "posterior_epred")
pred_LN  <- conditional_effects(model_hovs, effects = "Larval_nutrition", method = "posterior_epred")


pred_BS  <- as.data.frame(pred_BS[[1]])
pred_STI <- as.data.frame(pred_STI[[1]])
pred_FH  <- as.data.frame(pred_FH[[1]])
pred_LN  <- as.data.frame(pred_LN[[1]])


plot_BS_hover <- ggplot(pred_BS, aes(x = effect1__, y = estimate__)) +
  geom_point(data = Hoverfly_Traits %>% filter(!is.na(Adult_Body_size_num), !is.na(trend)),
             aes(x = Adult_Body_size_num, y = trend),
             alpha = 0.1, color = "gray50", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black",size=0.5) +
  geom_line(color = "#2ca02c", size = 1.5) +
  ylim (-0.2,0.2)+
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "#2ca02c", alpha = 0.2) +
  labs( x = "Body Size (ordinal)", y = "Predicted Trend") +
  theme_minimal()

plot_BS_hover


plot_STI_hover <- ggplot(pred_STI, aes(x = effect1__, y = estimate__)) +
  geom_point(data = Hoverfly_Traits %>% filter(!is.na(STI_Species_temperature_index), !is.na(trend)),
             aes(x = STI_Species_temperature_index, y = trend),
             alpha = 0.1, color = "gray50", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black",size=0.5) +
  geom_line(color = "#2ca02c", size = 1.5, linetype="dashed") +
  ylim (-0.2,0.2)+
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "#2ca02c", alpha = 0.2) +
  labs(x = "STI (ºC)", y = "") +
  theme_minimal()

plot_STI_hover

plot_FH_hover <- ggplot(pred_FH, aes(x = effect1__, y = estimate__)) +
  geom_jitter(data = Hoverfly_Traits %>% filter(!is.na(Flight_height), !is.na(trend)), 
              aes(x = Flight_height, y = trend),
              width = 0.2, alpha = 0.1, color = "gray50", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black",size=0.5) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__), width = 0.2, size=1, color = "#1f77b4") +
  geom_point(size = 4, color = "#1f77b4") +
  geom_point(size = 2, color = "white") + geom_point(size = 1, color = "#1f77b4") +
  labs(x = "Flight height", y = "") +
  ylim (-0.2,0.2)+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_FH_hover

plot_LN_hover <- ggplot(pred_LN, aes(x = effect1__, y = estimate__)) +
  geom_jitter(data = Hoverfly_Traits %>% filter(!is.na(Larval_nutrition), !is.na(trend)), 
              aes(x = Larval_nutrition, y = trend),
              width = 0.2, alpha = 0.1, color = "gray50", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black",size=0.5) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__), width = 0.2, size=1, color = "#1f77b4") +
  geom_point(size = 4, color = "#1f77b4") +
  geom_point(size = 2, color = "white") + geom_point(size = 1, color = "#1f77b4") +
  labs(x = "Larval Diet", y = "Predicted Trend") +
  ylim (-0.2,0.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_LN_hover

# Combine plots
plot_BS_hover + plot_STI_hover + plot_FH_hover + plot_LN_hover


# top row of 4 plots
top_row <- plot_ITD_bee + plot_STI_bee + plot_Soc_bee + plot_LDB_bee + plot_layout(ncol=4)

# bottom row of 4 plots
bottom_row <- plot_BS_hover + plot_STI_hover + plot_FH_hover + plot_LN_hover + plot_layout(ncol=4)

# combine with 4 columns and add vertical space

plot_traits <- top_row /
  plot_spacer() /
  bottom_row +
  plot_layout(heights = c(1, 0.2, 1))  # 0.05 is the spacer height

plot_traits


### Fig_5. Trends by IUCN Status and Genera.

Bee_IUCN <- Bee_Traits %>%
  filter(European_Category_2025 != "" & !is.na(European_Category_2025))

Bee_IUCN <- Bee_IUCN %>%
  mutate(IUCN_2025 = factor(European_Category_2025,
                            levels = c("DD","LC","NT","VU","EN","CR")))

plot_IUCN_bee <- ggplot(Bee_IUCN, aes(x = IUCN_2025, y = trend)) +
  geom_jitter(width = 0.2, alpha = 0.2) +
  geom_boxplot(fill = "lightblue", alpha = 0.7, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black",linewidth=0.5) +
  labs(x = "Bee IUCN Category", y = "Trend") +
  coord_cartesian(ylim = c(-0.33, 0.2)) +   # Adjust these limits to your data
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_IUCN_bee


# Boxplot for selected genera

genus_trend_summary <- Bee_Traits %>%
  group_by(Genus) %>%
  dplyr::summarise(
    median_trend = median(trend, na.rm = TRUE),
    sd_trend = sd(trend, na.rm = TRUE),
    n = n()) %>%
  arrange(desc(median_trend))  # Optional: Sort by trend

target_genera <- c("Hoplitis", "Tetralonia", "Osmia","Anthidium", "Amegilla",
                   "Sphecodes", "Seladonia","Dufourea", "Panurgus", "Melecta")

Bee_subset <- Bee_Traits %>% 
  filter(Genus %in% target_genera)

plot_genus_bee <- ggplot(Bee_subset, aes(x = fct_reorder(Genus, trend, .fun = median, .desc = TRUE), y = trend))+ 
  geom_jitter(width = 0.2, alpha = 0.2) +
  geom_boxplot(fill = "#fc8d62", alpha = 0.7, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black",size=0.5) +
  coord_cartesian(ylim = c(-0.2, 0.2)) +   # Adjust these limits to your data
  labs(x = "Bee genus", y = "Trend") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="italic"))

plot_genus_bee


# Same for hoverflies

Hoverfly_IUCN <- Hoverfly_Traits %>%
  filter(European_Category != "" & !is.na(European_Category))

Hoverfly_IUCN <- Hoverfly_IUCN %>%
  mutate(European_Category = factor(European_Category,
                            levels = c("DD","LC","NT","VU","EN","CR")))

plot_IUCN_hover <- ggplot(Hoverfly_IUCN, aes(x = European_Category, y = trend)) +
  geom_jitter(width = 0.2, alpha = 0.2) +
  geom_boxplot(fill = "lightblue", alpha = 0.7, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black",size=0.5) +
  labs(x = "Hoverfly IUCN Category", y = "") +
  coord_cartesian(ylim = c(-0.25, 0.3)) +   # Adjust these limits to your data
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_IUCN_hover


# Boxplot for selected genera

#Identify genera with more information.
Hoverfly_Traits %>%
  filter(!is.na(trend)) %>%
  count(Genus, sort = TRUE) %>%
  slice_max(n, n = 20)

top_bottom_genera <- Hoverfly_Traits %>%
  filter(!is.na(trend)) %>%
  group_by(Genus) %>%
  summarise(
    mean_trend = median(trend, na.rm = TRUE),
    n = n()
  ) %>%
  arrange(desc(mean_trend))

target_genera <- c("Syrphus", "Philhelius", "Dasysyrphus", "Melanostoma", "Chalcosyrphus",
                   "Melangyna", "Merodon", "Pipiza", "Parasyrphus", "Brachyopa")

Hoverfly_subset <- Hoverfly_Traits %>% 
  filter(Genus %in% target_genera)

plot_genus_hover <- ggplot(Hoverfly_subset, aes(x = fct_reorder(Genus, trend, .fun = median, .desc = TRUE), y = trend)) +
  geom_jitter(width = 0.2, alpha = 0.2) +
  geom_boxplot(fill = "#fc8d62", alpha = 0.7, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black",size=0.5) +
  coord_cartesian(ylim = c(-0.15, 0.25)) +
  labs(x = "Hoverfly genus", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,face="italic"))

plot_genus_hover


###Fig:5. Plot trends by IUCN Status and genera.

plot_iucn_and_genera <- (plot_IUCN_bee + plot_IUCN_hover) /
  (plot_genus_bee + plot_genus_hover) +
  plot_layout(heights = c(1, 1.2))

plot_iucn_and_genera

dev.off()




