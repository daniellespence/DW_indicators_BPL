###############################################################################
####### cH.1 GAMs -- data processing ##########
####### Danielle Spence ##########
####### Created 3/30/2021 ########
###############################################################################
### Clear memory

rm(list = ls())

library(pacman)
p_load(tidyverse, reshape2, dplyr, tidyr, readxl,  janitor,  install = TRUE)
library(dplyr)

setwd("C:/Users/danis/OneDrive/R/DW Indicators")

## Read in data 

df <- read_csv("Data/Raw data/bpwtp_labdat_db.csv")
## pull out response variable and potential explanatory variables 

## select only raw data

df <- df%>%
  dplyr::select(date_ymd, parameter, result, station, year)

df<- df%>%
  filter(station == "Raw",
         year %in% c(2006:2022))


df<- df%>%
  filter(parameter %in% c("Nitrate","Ammonia N","Phosphate (ortho)","Phosphate (total)","Temperature","Organic N", "TDS"))

df<- df%>% dplyr::select(-station)

df$date_ymd<- as.Date(df$date_ymd, "%m/%d/%Y")

 # parameter_summary <- df %>%
 #       dplyr::filter(parameter == "TDS") %>%
 #   dplyr::summarise(
 #             mean = mean(result, na.rm = TRUE),
 #             median = median(result, na.rm = TRUE),
 #             sd = sd(result, na.rm = TRUE),
 #             min = min(result, na.rm = TRUE),
 #             max = max(result, na.rm = TRUE),
 #             count = n()
 #         )
 # 
 #  print(parameter_summary)

  # mean median    sd   min   max count
  # <dbl>  <dbl> <dbl> <dbl> <dbl> <int>
  # 421.    392  115.   250   810   284


## Pivot wider to make parameters into columns


d2<- df %>%
  tidyr::pivot_wider(names_from = "parameter", values_from = "result",
                     values_fn = mean)


df <- d2%>%
  dplyr::rename(
    Date = "date_ymd",
    
    NO3_mg.L = "Nitrate",
    NH3_mg.L = "Ammonia N",
    SRP_ug.L = "Phosphate (ortho)",
    TP_ug.L = "Phosphate (total)",
    Temp_C = "Temperature",
    Org_N_mg.L = "Organic N"
    )


df$Date<- as.Date(df$Date, "%m/%d/%Y")
df <- mutate(df,
             nMonth = as.numeric(format(Date,'%m')))

df <- df[rowSums(is.na(df)) <= 5, ]  


df %>% filter(is.na(NH3_mg.L)) #53 NAs
df %>% filter(is.na(NO3_mg.L)) #33 NAs
df %>% filter(is.na(SRP_ug.L)) #34 NAs
df %>% filter(is.na(TDS)) #21 NAs

## Add flow d## Add flow d## Add flow data

flow <- read_csv("Data/daily flow.csv")

flow<- flow %>%
  select(date, combined_05JG004.cms, SK05JG006.cms, RC_IC_cms)%>%
  dplyr::rename("Date" = date)


full_df<- merge(df, flow, by="Date", all=TRUE)


full_df <- full_df[rowSums(is.na(full_df)) <= 4, ]       # Apply rowSums & is.na


full_df <- mutate(full_df,
                  DOY = as.numeric(format(Date,'%j')),
                  nMonth = as.numeric(format(Date,'%m')))


# add anthony imputations from 2003-2004 for SRP and ammonia

#impute from nearest neighbours

# # ammonia
 full_df %>% filter(year == 2012) # c(0.02, 0.19) for months 10,11,12
 full_df %>% filter(year == 2013)
 full_df %>% filter(year == 2022)

full_df <- full_df %>% 
  mutate(NH3_mg.L = case_when(
    is.na(NH3_mg.L) & year == c(2012) & nMonth %in% c(10,11,12)~ mean(c(0.02, 0.19)),
    is.na(NH3_mg.L) & year == c(2022) & nMonth %in% c(5,6,7,8,9,10)~ mean(c(0.0, 0.00)),
    TRUE ~ as.numeric(NH3_mg.L)
  ))




sum(full_df$NH3_mg.L == 0, na.rm = TRUE) #73
sum(full_df$NO3_mg.L == 0, na.rm = TRUE) #96
sum(full_df$SRP_ug.L == 0, na.rm = TRUE) #34


## change 0s for N and P to limits of detection

full_df$NO3_mg.L[full_df$NO3_mg.L == 0] <- (0.057/2)
full_df$NH3_mg.L[full_df$NH3_mg.L == 0] <- (0.086/2)
full_df$SRP_ug.L[full_df$SRP_ug.L == 0] <- (3/2)


# 
# parameter_summary <- full_df %>%
#   dplyr::summarise(
#     mean = mean(TDS, na.rm = TRUE),
#     median = median(TDS, na.rm = TRUE),
#     sd = sd(TDS, na.rm = TRUE),
#     min = min(TDS, na.rm = TRUE),
#     max = max(TDS, na.rm = TRUE),
#     count = n()
#   )
# 
# print(parameter_summary)

# mean        median       sd     min  max count
# 420.9437    392      114.7851    250 810   305


write.csv(full_df, file="data/bpgamdataCLEAN_TDS.csv", row.names=F)

