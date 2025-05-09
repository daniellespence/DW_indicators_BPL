###############################################################################
####### BPWTP Data -- GAMs ##########
####### Danielle Spence ##########
####### Created 5/4/2023 ########
###############################################################################
### Clear memory
rm(list = ls())


setwd("C:/Users/danis/OneDrive/R/DW indicators")

library(pacman)
p_load(tidyverse, ggplot2, mgcv, gratia, readr, ggExtra, dplyr, lubridate,
       cowplot, tibble, patchwork, install = TRUE)


effects_table<- read_csv("C:/Users/danis/OneDrive/R/DW Indicators/Data/Effects table no temp.csv")

papertheme <- theme_bw(base_size=14, base_family = 'Arial') 


# Ensure your data is in wide format
heatmap_data <- effects_table %>%
  select(Predictor, Response, Rank) %>%
  spread(key = Response, value = Rank)


# Convert back to long format for ggplot
heatmap_long <- gather(heatmap_data, key = "Response", value = "Rank", -Predictor)

heatmap_long$Predictor<- factor(heatmap_long$Predictor, levels = c("PDO_SOI", "Lake Diefenbaker flow", "Soluble reactive phosphorus",  "Dissolved inorganic nitrogen"))

ggplot(heatmap_long, aes(x = Predictor, y = Response, fill = Rank)) +
  geom_tile(color = "white") +
  papertheme+
  theme_minimal(base_size = 14) +
  scale_fill_gradientn(colors = c("darkred", "lightyellow"), 
                       values = scales::rescale(c(1, 2, 3)), 
                       breaks = c(1, 2, 3),  # Ensure only whole numbers appear
                       guide = guide_colorbar(reverse = TRUE)) +
  labs(x = "Predictor",
       y = "Indicator",
       fill = "Rank")+
scale_x_discrete(labels = c("Soluble reactive phosphorus" = "SRP", 
                            "Dissolved inorganic nitrogen" = "DIN", 
                            #"Water temperature" = "W_Temp",
                            "Lake Diefenbaker flow" = "QLD",
                            "PDO_SOI" = "PDO_SOI"))


ggsave("Output/Predictor_Strength_Heatmap TDS no temp.png", width = 8, height = 6)

