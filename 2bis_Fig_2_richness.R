## Fig 2_changes in richness.

#Install and load required packages
rm(list=ls())
pkgs <- c("ggplot2", "stringr") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

### Plot change in bee richness at n = 100----

plot_data <- read.csv("data/final_and_intermediate_outputs/Plot_data_bees_100.csv")

Dataset_Richness_bees <- read.csv("data/final_and_intermediate_outputs/Dataset_Richness_bees_100.csv")

## Plot richness bees n = 100.

# Define colors for each region
region_colors <- c(
  "alpine" = "#44AA99",       # Blue
  "atlantic" = "#117733",     # Orange
  "boreal" = "#332288",       # Green
  "continental" = "#CC6677",  # Red
  "mediterranean" = "#DDCC77" # Purple
)

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
  labs(x = "Year", y = "Predicted bee richness", color = "Bioregion") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

bee_plot_100


### Plot change in hoverfly richness at n = 100----

plot_data <- read.csv("data/final_and_intermediate_outputs/Plot_data_hover_100.csv")
Dataset_Richness_hover <- read.csv("data/final_and_intermediate_outputs/Dataset_Richness_hover_100.csv")

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
  labs(x = "Year", y = "Predicted hoverfly richness", color = "Bioregion") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

hoverfly_plot_100


### Plots for supplementary material (n = 200).

### Plot of bees richness at n = 200----

plot_data_bee_200 <- read.csv("data/final_and_intermediate_outputs/Plot_data_bees_200.csv")
Dataset_Richness_bees <- read.csv("data/final_and_intermediate_outputs/Dataset_Richness_bees_200.csv")

bee_plot_200 <- ggplot(plot_data_bee_200, aes(x = endYear, y = .estimate)) +
  
  # Raw points
  geom_point(data = Dataset_Richness_bees,
             aes(x = endYear, y = Estimator, color = region_100, fill = region_100),
             alpha = 0.3, size = 1.5, inherit.aes = FALSE) +
  
  # Increasing / Decreasing = solid
  geom_line(
    data = subset(plot_data_bee_200, sig != "Stable"),
    aes(color = region_100, group = interaction(region_100, sig)),
    linetype = "solid",
    size = 1.2
  ) +
  
  # Stable = dashed
  geom_line(
    data = subset(plot_data_bee_200, sig == "Stable"),
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
  labs(x = "Year", y = "Predicted bee richness", color = "Bioregion") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

bee_plot_200


### Plot change in hoverfly richness at n = 200----

plot_data <- read.csv("data/final_and_intermediate_outputs/Plot_data_hover_200.csv")
Dataset_Richness_hover <- read.csv("data/final_and_intermediate_outputs/Dataset_Richness_hover_200.csv")

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
  labs(x = "Year", y = "Predicted hoverfly richness", color = "Bioregion") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

hoverfly_plot_200



##Combined plots

bee_plot_100 + hoverfly_plot_100

bee_plot_200 + hoverfly_plot_200

