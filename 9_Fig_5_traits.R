### Fig_5.

library(patchwork)
library(tidyverse)

#### Traits modelled with brms----

load("traits_models_brms_phylo.RData")
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
