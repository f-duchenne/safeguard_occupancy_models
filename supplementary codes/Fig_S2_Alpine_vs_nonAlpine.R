
library(tidyverse)

### Fig_S2_Test trends of alpine vs. non-alpine species in the alpine region.

#Load trends for species and rename first column.
Trends <- read.csv("all_trends.csv")
colnames(Trends)[1]<- "Species"

library(tidyverse)
#Join list of alpine species, designated by experts.

#Hoverflies
Alpine_expert_hoverflies <- read.csv ("data/raw_data/Summary_records_and_trends_AV.csv", sep=";")

#Create a lookup for expert alpine designation
Alpine_expert_hoverflies <- Alpine_expert_hoverflies %>%
  mutate(Alp_exp = ifelse(X.2 == "Alp", "Yes", "No")) %>%
  select(TAXON, Alp_exp)

colnames(Alpine_expert_hoverflies)[1]<- "Species"

#Join with trends data
Trends <- Trends %>%
  left_join(Alpine_expert_hoverflies, by = "Species")

Trends %>% filter (Alp_exp == "Yes") %>% distinct(Species) %>% count()

#Add alpine Bees
Alpine_expert_bees <- read.csv("data/raw_data/Alpine_bee_species.csv")
nrow(Alpine_expert_bees)  # How many Alpine species

Trends$Alp_exp[Trends$Species %in% Alpine_expert_bees$scientificName] <- "Yes"


## Subset trends only for the alpine region, selceting two different baselines.
Trends_alpine_1981 <- Trends %>% filter (baseline == "1981" & region_50 == "alpine")
Trends_alpine_1981 <- Trends_alpine_1981 %>% filter(trend > -1, trend < 1)


#Plot
p1 <- ggplot(Trends_alpine_1981, aes(y = trend, x = Alp_exp, fill = Alp_exp,show.legend = F)) +
  geom_jitter(alpha = 0.1, aes(fill = Alp_exp))+
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  #scale_fill_manual(values = c("No" = "#F8766D", "Yes" = "#00BFC4")) +
  ylim(-0.5, 0.5) +
  labs(
    x = "Alpine species",
    y = "Trend",
    fill = "Alpine",
    title = "Baseline 1981"
  ) +
  facet_wrap(~ taxo_group) +  # <<–– add this line
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold")
  )

p1

Trends_alpine_1981_bees <- Trends_alpine_1981 %>% filter (taxo_group == "bees")
Trends_alpine_1981_hoverflies <- Trends_alpine_1981 %>% filter (taxo_group == "hoverflies")


model1 <- lm ( data = Trends_alpine_1981_bees,
                 formula = trend ~ Alp_exp,
                 weights = (1/Trends_alpine_1981_bees$sde)^2)

hist(resid(model1), breaks = 100)
summary(model1)

emmeans(model1, "Alp_exp")

model2 <- lm ( data = Trends_alpine_1981_hoverflies,
               formula = trend ~ Alp_exp,
               weights = (1/Trends_alpine_1981_hoverflies$sde)^2)

hist(resid(model2), breaks = 100)
summary(model2)

emmeans(model2, "Alp_exp")
