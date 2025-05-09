###############################################################################
####### BPWTP Data -- GAMs ##########
####### Danielle Spence ##########
####### Created 5/4/2023 ########
###############################################################################
### Clear memory
rm(list = ls())

library(pacman)
p_load(tidyverse, ggplot2, mgcv, gratia, readr, GGally, dplyr, lubridate, readr,
       cowplot, tibble,viridis, install = TRUE)

setwd("C:/Users/danis/OneDrive/R/DW indicators")

papertheme <- theme_bw(base_size=12, base_family = 'Arial') 

##----------------------------------------------------------------------------##
## 1. Read in data 
##----------------------------------------------------------------------------##

df <- read_csv("data/bpgamdataCLEAN_Turb.csv")

df$Date<- as.Date(df$Date, "%m/%d/%Y")

# add in Year and nMonth for numeric month and a proper Date class

df <- mutate(df,
             year = as.numeric(format(Date,'%Y')),
             DOY = as.numeric(format(Date,'%j')),
             nMonth = as.numeric(format(Date,'%m')),
             week = as.numeric(format(Date, '%W')))%>% 
  filter(year %in% c(1990:2022))
# 974 obs


## order the data
df <- df[with(df, order(Date)),]


# rename for ease of use
df <- df%>%
  dplyr::rename(
    date = "Date",
    turb = "Turbidity",
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


#Remove rows without turb data
df <- df[complete.cases(df$turb),]
# 974 obs


chla <- read_csv("data/bpgamdataCLEAN_chla.csv")
chla$Date<- as.Date(chla$Date, "%m/%d/%Y")
# 1244 obs


# rename for ease of use
chla <- chla%>%
  dplyr::rename(
    date = "Date",
    chla = "Chla_ug.L")


#Remove rows without Chla data
chla <- chla[complete.cases(chla$chla),]
# 1238 obs

df1<- left_join(df, chla, by="date")

# select only ice-free months

#df1 <- df1%>% 
 # filter(nMonth.x %in% c(5, 6, 7, 8, 9, 10))

# Fit linear model
lm_model <- lm(turb ~ chla, data = df1)
summary(lm_model)

# Extract the coefficients from the model
intercept <- coef(lm_model)[1]
slope <- coef(lm_model)[2]


# Create the equation string
equation <- paste0("y = ", format(slope, digits = 2), "x + ", format(intercept, digits = 2))

# Extract R² value
r_squared <- summary(lm_model)$r.squared
r_squared_text <- paste("R² =", round(r_squared, 3))


# Extract R² value
r_squared <- summary(lm_model)$r.squared
r_squared_text <- paste("R² =", round(r_squared, 3))


p1<- ggplot(df1, aes(x = turb, y = chla)) +
  geom_point(color = "black", alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  #annotate("text", x = max(df1$Conductivity) * 0.7, y = max(df1$TDS) * 0.9,
  #        label = r_squared_text, size = 5, color = "black") +
  annotate("text", x = 40, y = 150, label = paste0(equation, "\n ", format(r_squared_text))) + # Add equation and R-squared
  labs(x = "Turbidity (NTU)", 
       y = expression(paste("Chlorophyll ", italic("a"), " (", mu, "g L"^-1*")"))) +
  theme_bw()
p1

ggsave('output/chla vs turbidity.png', p1, height = 6, width  = 8)


?cor.te
cor.test(df1$turb, df1$chla, method = "spearman")

cor.test(df1$turb, df1$chla, method = "kendall")

# assess correlations

dfx<-df1%>%
  select(turb, chla)
## visualize correlations
p<- ggpairs(dfx[])+ theme_bw()  # can see that turb and turbidity, temperature, Org_N and TP are correlated. Ammonia and conductivity are not. (NH3, NO3, and TP correlated; TP and temp). Org N and TP are similarly correlated with turb (0.640 and 0.667, respectively)
p


ggsave('output/turb vs chla Correlations.png', p, height = 8, width  = 10)

