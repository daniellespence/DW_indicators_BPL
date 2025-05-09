###############################################################################
####### BPWTP Data -- GAMs ##########
####### Danielle Spence ##########
####### Created 5/4/2023 ########
###############################################################################
### Clear memory
rm(list = ls())

library(pacman)
p_load(tidyverse, ggplot2, mgcv, gratia, readr, GGally, dplyr, lubridate,
       cowplot, tibble,viridis, install = TRUE)

setwd("C:/Users/danis/OneDrive/R/DW indicators")

papertheme <- theme_bw(base_size=14, base_family = 'Arial') 

##----------------------------------------------------------------------------##
## 1. Read in data 
##----------------------------------------------------------------------------##

df <- read_csv("data/bpgamdataCLEAN_Conductivity.csv")

df$Date<- as.Date(df$Date, "%m/%d/%Y")

# add in Year and nMonth for numeric month and a proper Date class

df <- mutate(df,
             year = as.numeric(format(Date,'%Y')),
             DOY = as.numeric(format(Date,'%j')),
             nMonth = as.numeric(format(Date,'%m')),
             week = as.numeric(format(Date, '%W')))%>% 
  filter(year %in% c(1990:2022))
# 1039 obs


## order the data
df <- df[with(df, order(Date)),]


# rename for ease of use
df <- df%>%
  dplyr::rename(
    date = "Date",
    Conductivity = "Conductivity",
    NO3 = "NO3_mg.L",
    NH3 = "NH3_mg.L",
    SRP = "SRP_ug.L",
    TP = "TP_ug.L",
    W_temp = "Temp_C",
    QBP = "combined_05JG004.cms",
    QLD = "SK05JG006.cms",
    QWS = "RC_IC_cms",
    orgN = "Org_N_mg.L"
  )


#Remove rows without Conductivity data
df <- df[complete.cases(df$Conductivity),]
# 1022 obs



tds <- read_csv("data/bpgamdataCLEAN_TDS.csv")

# 305 obs

tds$Date<- as.Date(tds$Date, "%m/%d/%Y")

# add in Year and nMonth for numeric month and a proper Date class


## order the data
tds <- tds[with(tds, order(Date)),]


# rename for ease of use
tds <- tds%>%
  dplyr::rename(
    date = "Date",
    NO3 = "NO3_mg.L",
    NH3 = "NH3_mg.L",
    SRP = "SRP_ug.L",
    TP = "TP_ug.L",
    W_temp = "Temp_C",
    QBP = "combined_05JG004.cms",
    QLD = "SK05JG006.cms",
    QWS = "RC_IC_cms",
    orgN = "Org_N_mg.L",
  )

#Remove rows without Chla data
tds <- tds[complete.cases(tds$TDS),]

tds<- tds %>% select(date, TDS)
# 284 obs


df1<- left_join(df, tds, by="date")

# Fit linear model
lm_model <- lm(TDS ~ Conductivity, data = df1)
summary(lm_model)

# Extract the coefficients from the model
intercept <- coef(lm_model)[1]
slope <- coef(lm_model)[2]


# Create the equation string
equation <- paste0("y = ", format(slope, digits = 2), "x + ", format(intercept, digits = 2))

# Extract R² value
r_squared <- summary(lm_model)$r.squared
r_squared_text <- paste("R² =", round(r_squared, 3))


p1<- ggplot(df1, aes(x = Conductivity, y = TDS)) +
  geom_point(color = "black", alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  scale_x_continuous(limits=c(400, 1200))+
  #annotate("text", x = max(df1$Conductivity) * 0.7, y = max(df1$TDS) * 0.9,
   #        label = r_squared_text, size = 5, color = "black") +
  annotate("text", x = 500, y = 750, label = paste0(equation, "\n ", format(r_squared_text))) + # Add equation and R-squared
  labs(x = "Conductivity (µS/cm)", 
       y = "TDS (mg/L)") +
  theme_bw()
p1

ggsave('output/TDS vs conductivity.png', p1, height = 6, width  = 8)

# are they correlated?
cor.test(df1$Conductivity, df1$TDS, method = "spearman")

cor.test(df1$Conductivity, df1$TDS, method = "kendall")


# Predict TDS where it's missing
merged_data <- df1 %>%
  mutate(Predicted_TDS = ifelse(is.na(TDS), predict(lm_model, newdata = df1), TDS))


p2<- ggplot(merged_data, aes(x = Predicted_TDS, y = TDS)) +
  geom_point(color = "black", alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  annotate("text", x = max(df1$Conductivity) * 0.7, y = max(df1$TDS) * 0.9,
           label = r_squared_text, size = 5, color = "black") +
  labs(x = "Predicted TDS (mg/L)", 
       y = "TDS (mg/L)") +
  theme_minimal()
p2

ggsave('output/predicted vs measured TDS.png', p2, height = 6, width  = 8)


merged_data <- merged_data %>%
  mutate(TDS_Type = ifelse(is.na(TDS), "Predicted", "Measured"))


p3 <- ggplot(merged_data, aes(x = year)) +
  # Measured TDS points and line
  geom_point(aes(y = TDS, color = "Measured", shape = "Measured"), size = 3, show.legend = TRUE) +  # Show legend for points
  geom_line(aes(y = TDS), color = "blue", linetype = "solid", show.legend = FALSE) +  # Don't show legend for line
  # Predicted TDS points and line
  geom_point(aes(y = Predicted_TDS, color = "Predicted", shape = "Predicted"), size = 2, shape = 17, show.legend = TRUE) +  # Show legend for points
  geom_line(aes(y = Predicted_TDS), color = "red", linetype = "dashed", show.legend = FALSE) +  # Don't show legend for line
  labs(title = "Comparison of Measured and Predicted TDS Over Time",
       y = "TDS (mg/L)", x = "Year",
       color = "TDS Type",  # Legend title for color
       shape = "TDS Type") +  # Legend title for shape
  scale_color_manual(values = c("Measured" = "blue", "Predicted" = "red")) +
  scale_shape_manual(values = c("Measured" = 16, "Predicted" = 17)) +
  theme_minimal()

p3

ggsave('output/predicted vs measured TDS timeseries.png', p3, height = 6, width  = 8)


ggplot(merged_data, aes(x = Predicted_TDS, y = TDS, color = TDS_Type, shape = TDS_Type)) +
  geom_point(size = 3, alpha = 0.7) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # 1:1 line
  scale_color_manual(values = c("red", "blue")) +
  scale_shape_manual(values = c(17, 16)) +
  labs(title = "Predicted vs. Measured TDS",
       x = "Predicted TDS (mg/L)",
       y = "Measured TDS (mg/L)",
       color = "TDS Type",
       shape = "TDS Type") +
  theme_minimal()


# Blank plot with two geom_point layers
ggplot() +
  # Measured TDS points (blue circles)
  geom_point(data = merged_data %>% filter(TDS_Type == "Measured"), 
             aes(x = year, y = TDS), color = "blue", size = 3) +
  # Predicted TDS points (red triangles)
  geom_point(data = merged_data %>% filter(TDS_Type == "Predicted"), 
             aes(x = year, y = Predicted_TDS), color = "red", size = 3, shape = 17) +
  # Labels and title
  labs(title = "Measured vs. Predicted TDS",
       x = "Year",
       y = "TDS (mg/L)") +
  theme_minimal()
