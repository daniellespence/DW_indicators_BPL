###############################################################################
####### cH.1 GAMs -- data processing ##########
####### Danielle Spence ##########
###############################################################################
### Clear memory

rm(list = ls())

library(pacman)
p_load(tidyverse, reshape2, dplyr, tidyr, readxl,  janitor,  install = TRUE)
library(dplyr)

setwd("C:/Users/danis/OneDrive/R/DW Indicators")

## Read in data 

df <- read_csv("Data/Raw data/bpwtp_labdat_db.csv")%>%
  filter(year %in% c(1990:2022))
## pull out response variable and potential explanatory variables 

## select only raw data

df <- df%>%
  dplyr::select(date_ymd, parameter, result, station, year)

df<- df%>%
  filter(station == "Raw",
         year %in% c(1984:2022))


df<- df%>%
  filter(parameter %in% c("Nitrate","Ammonia N","Phosphate (ortho)","Phosphate (total)","Temperature","Organic N", "Odour"))

df<- df%>% dplyr::select(-station)

df$date_ymd<- as.Date(df$date_ymd, "%m/%d/%Y")

#  parameter_summary <- df %>%
#        filter(parameter == "Odour") %>%
#        summarise(
#              mean = mean(result, na.rm = TRUE),
#              median = median(result, na.rm = TRUE),
#              sd = sd(result, na.rm = TRUE),
#              min = min(result, na.rm = TRUE),
#              max = max(result, na.rm = TRUE),
#              count = n()
#          )
#  
#   print(parameter_summary)
# # A tibble: 1 × 6
# mean median    sd   min   max count
# <dbl>  <dbl> <dbl> <dbl> <dbl> <int>
# 71.1     43  88.8     2   800  1938

# add missing  data


missing_2001 <- read_excel("Data/Raw data/ROUTINE LAB DATA 2001(v2).xlsx",
                           sheet = "Weekly Data", col_names = TRUE, range = "A8:BD67") %>% 
  filter(Parameters %in% c("Nitrate","Ammonia N","Phosphate (ortho)","Phosphate (total)","Temperature","Organic N", "Odour")) %>% 
  select(-c("...3":"...4")) %>% 
  pivot_longer(cols = -c(Parameters, Units), names_to = "date_ymd", values_to = "result") %>% 
  mutate_at(c("date_ymd", "result"), as.numeric) %>% 
  mutate(date_ymd = excel_numeric_to_date(date_ymd),
         year = year(date_ymd)) %>% 
  select(parameter = Parameters, date_ymd, year, result) 

bp_longterm_sans2001 <- df %>% filter(!year == 2001) 



bp_longterm_infill <- bind_rows(bp_longterm_sans2001, missing_2001) %>% 
  arrange(date_ymd)



## Pivot wider to make parameters into columns


d2<- bp_longterm_infill %>%
  tidyr::pivot_wider(names_from = "parameter", values_from = "result",
                     values_fn = mean)


 # parameter_summary <- d2 %>%
 #       dplyr::summarise(
 #             mean = mean(Odour, na.rm = TRUE),
 #             median = median(Odour, na.rm = TRUE),
 #             sd = sd(Odour, na.rm = TRUE),
 #             min = min(Odour, na.rm = TRUE),
 #             max = max(Odour, na.rm = TRUE),
 #             count = n()
 #         )
 # 
 #  print(parameter_summary)

  # # A tibble: 1 × 6
  # mean median    sd   min   max count
  # <dbl>  <dbl> <dbl> <dbl> <dbl> <int>
  #   1  71.0     43  88.5     2   800  2029
  
d2 <- d2%>%
  dplyr::rename(
    Date = "date_ymd",
    Odour_TON = "Odour",
    NO3_mg.L = "Nitrate",
    NH3_mg.L = "Ammonia N",
    SRP_ug.L = "Phosphate (ortho)",
    TP_ug.L = "Phosphate (total)",
    Temp_C = "Temperature",
    Org_N_mg.L = "Organic N"
    )

# select parameters of interest

df <- d2%>%
  select(Date, year, NO3_mg.L, NH3_mg.L, SRP_ug.L, TP_ug.L, Temp_C,  Org_N_mg.L, Odour_TON)

df$Date<- as.Date(df$Date, "%m/%d/%Y")

df <- mutate(df,
             nMonth = as.numeric(format(Date,'%m')))

df <- df[rowSums(is.na(df)) <= 4, ]  

# df %>% filter(is.na(NH3_mg.L)) #871 NAs
# df %>% filter(is.na(NO3_mg.L)) #788 NAs
# df %>% filter(is.na(SRP_ug.L)) #868 NAs
# df %>% filter(is.na(Odour_TON)) #72 NAs

## Add flow d## Add flow d## Add flow data

flow <- read_csv("Data/daily flow.csv")

flow<- flow %>%
  select(date, combined_05JG004.cms, SK05JG006.cms, RC_IC_cms)%>%
  dplyr::rename("Date" = date)


full_df<- merge(df, flow, by="Date", all=TRUE)


full_df <- full_df[rowSums(is.na(full_df)) <= 7, ]       # Apply rowSums & is.na


full_df <- mutate(full_df,
                  DOY = as.numeric(format(Date,'%j')),
                  nMonth = as.numeric(format(Date,'%m')))



# add anthony imputations from 2003-2004 for SRP and ammonia

#impute from nearest neighbours

# # ammonia
# full_df %>% filter(year == 2001) # c(0.12, 0.07) for month 3 I get (0.12, 0.07)
# full_df %>% filter(year == 2003) # c(0.04, 0) for month 1,2
# full_df %>% filter(year == 2003) # c(0.1, 0.02) for month 11,12
# full_df %>% filter(year == 2004) # c(0.1, 0.02) for month 1,2,3


full_df <- full_df %>% 
  mutate(NH3_mg.L = case_when(
    is.na(NH3_mg.L) & year == c(2002) & nMonth %in% c(3)~ mean(c(0.12, 0.07)),
    is.na(NH3_mg.L) & year == c(2003) & nMonth %in% c(1, 2,3) ~ mean(c(0.04, 0)),
    is.na(NH3_mg.L) & year == c(2003) & nMonth %in% c(11, 12) ~ mean(c(0.1, 0.02)),
    is.na(NH3_mg.L) & year == c(2004) & nMonth %in% c(1, 2, 3)~ mean(c(0.1, 0.02)),
    is.na(NH3_mg.L) & year == c(2022) & nMonth %in% c(5,6,7,8,9,10)~ mean(c(0.0, 0.00)),
    TRUE ~ as.numeric(NH3_mg.L)
  ))


## SRP

# full_df %>% filter(year %in% 2004) # c(10, 5.82)


full_df <- full_df %>% 
  mutate(SRP_ug.L = case_when(
    is.na(SRP_ug.L) & year %in% c(2003) & nMonth %in% c(12) ~ mean(c(10, 5.82)),
    is.na(SRP_ug.L) & year %in% c(2004) & nMonth %in% c(1,2,3,4,5,6,7) ~ mean(c(10, 5.82)),
    TRUE ~ as.numeric(SRP_ug.L)
  ))


# OrgN
# full_df %>% filter(year == 2002) # c(0.012, 0.07) for month 3
# full_df %>% filter(year == 2003) # c(0.73, 0.66) for month 1,2,3
# full_df %>% filter(year == 2003) # c(0.75, 0.053) for month 11,12
# full_df %>% filter(year == 2004) # c(0.75, 0.053) for month 1,2,3

full_df <- full_df %>% 
  mutate(Org_N_mg.L = case_when(
    is.na(Org_N_mg.L) & year == c(2002) & nMonth %in% c(3)~ mean(c(0.012, 0.07)),
    is.na(Org_N_mg.L) & year == c(2003) & nMonth %in% c(1, 2) ~ mean(c(0.73, 0.66)),
    is.na(Org_N_mg.L) & year == c(2003) & nMonth %in% c(11, 12) ~ mean(c(0.75, 0.053)),
    is.na(Org_N_mg.L) & year == c(2004) & nMonth %in% c(1, 2, 3)~ mean(c(0.75, 0.053)),
    TRUE ~ as.numeric(Org_N_mg.L)
  ))


## TP

full_df %>% filter(year %in% 2004) # c(45, 94.29)



full_df <- full_df %>% 
  mutate(TP_ug.L = case_when(
    is.na(TP_ug.L) & year %in% c(2003:2004) ~ mean(c(45, 94.29)),
    TRUE ~ as.numeric(TP_ug.L)
  ))



sum(full_df$NH3_mg.L == 0, na.rm = TRUE) #289
sum(full_df$NO3_mg.L == 0, na.rm = TRUE) #532
sum(full_df$SRP_ug.L == 0, na.rm = TRUE) #177

## change 0s for N and P to limits of detection

full_df$NO3_mg.L[full_df$NO3_mg.L == 0] <- (0.057/2)
full_df$NH3_mg.L[full_df$NH3_mg.L == 0] <- (0.086/2)
full_df$SRP_ug.L[full_df$SRP_ug.L == 0] <- (3/2)


# parameter_summary <- full_df %>%
#   dplyr::summarise(
#     mean = mean(Odour_TON, na.rm = TRUE),
#     median = median(Odour_TON, na.rm = TRUE),
#     sd = sd(Odour_TON, na.rm = TRUE),
#     min = min(Odour_TON, na.rm = TRUE),
#     max = max(Odour_TON, na.rm = TRUE),
#     count = n()
#   )
# 
# print(parameter_summary)
# mean median       sd min max count
# 1 71.00393     43 88.50783   2 800  2029


write.csv(full_df, file="data/bpgamdataCLEAN_TON.csv", row.names=F)


