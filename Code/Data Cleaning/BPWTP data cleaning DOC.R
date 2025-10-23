###############################################################################
####### cH.1 GAMs -- data processing ##########
####### Danielle Spence ##########
####### Created 3/30/2021 ########
###############################################################################
### Clear memory

rm(list = ls())

library(pacman)
p_load(tidyverse, reshape2, dplyr, tidyr, readxl, janitor,
       install = TRUE)
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
         year %in% c(1984:2022))


df<- df%>%
  filter(parameter %in% c("Nitrate","Ammonia N","Phosphate (ortho)","Phosphate (total)","Temperature","Organic N", "DOC",  "TOC"))


df$date_ymd<- as.Date(df$date_ymd, "%m/%d/%Y")

df<- df%>% dplyr::select(-station)

 parameter_summary <- df %>%
       filter(parameter == "DOC") %>%
       dplyr::summarise(
             mean = mean(result, na.rm = TRUE),
             median = median(result, na.rm = TRUE),
             sd = sd(result, na.rm = TRUE),
             min = min(result, na.rm = TRUE),
             max = max(result, na.rm = TRUE),
             count = n()
         )

   print(parameter_summary)
   
   # # A tibble: 1 × 6
   # mean median    sd   min   max count
   # <dbl>  <dbl> <dbl> <dbl> <dbl> <int>
   # 6.63    6.2  1.82  1.96    15  2985

# add missing 2001 data


missing_2001 <- read_excel("Data/Raw data/ROUTINE LAB DATA 2001(v2).xlsx",
                           sheet = "Weekly Data", col_names = TRUE, range = "A8:BD67") %>% 
  filter(Parameters %in% c("Nitrate","Ammonia N","Phosphate (ortho)","Phosphate (total)","Temperature","Organic N", "Raw TOC", "Raw DOC (GF diss)")) %>% 
  select(-c("...3":"...4")) %>% 
  pivot_longer(cols = -c(Parameters, Units), names_to = "date_ymd", values_to = "result") %>% 
  mutate_at(c("date_ymd", "result"), as.numeric) %>% 
  mutate(date_ymd = excel_numeric_to_date(date_ymd),
         year = year(date_ymd)) %>% 
  select(parameter = Parameters, date_ymd, year, result) 

missing_2001 <- missing_2001 %>%
  mutate(parameter = case_when(
    parameter == "Raw DOC (GF diss)" ~ "DOC", 
    parameter == "Raw TOC" ~ "TOC",
    TRUE ~ parameter))


bp_longterm_sans2001 <- df %>% filter(!year == 2001) 



bp_longterm_infill <- bind_rows(bp_longterm_sans2001, missing_2001) %>% 
  arrange(date_ymd)



## Pivot wider to make parameters into columns


d2<- bp_longterm_infill %>%
  tidyr::pivot_wider(names_from = "parameter", values_from = "result",
                     values_fn = mean)


### infill missing values from the 90s by modelling off TOC

# 38 NAs not including 1991—1993
bp_longterm_infill_NA_sans9193 <- d2 %>% 
  filter(is.na(DOC), !year %in% c(1991:1993)) %>% 
  mutate(rownum = row_number()) %>% 
  select(rownum, everything())

df %>% 
  filter(parameter == "TOC") %>% 
  ggplot(aes(yday(date_ymd), result)) +
  facet_wrap(~ year) + 
  geom_point() +
  geom_point(data = bp_longterm_infill, aes(yday(date_ymd), result), col = "green")

bp_doc <- bp_longterm_infill %>% 
  filter(parameter == "DOC") %>% 
  rename(DOC = result)

bp_toc <- df %>% 
  filter(parameter == "TOC") %>% 
  mutate(date_ymd = as.character(date_ymd)) %>% 
  add_row(date_ymd = "1996-12-30", result = NA) %>%
  add_row(date_ymd = "2001-12-31", result = NA) %>% 
  add_row(date_ymd = "2002-12-30", result = NA) %>% 
  add_row(date_ymd = "2018-01-01", result = NA) %>% 
  add_row(date_ymd = "2019-01-07", result = NA) %>% 
  rename(TOC_mg.L = result) %>% 
  mutate(date_ymd = ymd(date_ymd),
         year = year(date_ymd),
         month = month(date_ymd),
         week = week(date_ymd)) %>% 
  arrange(date_ymd) %>% 
  select(date_ymd, TOC_mg.L)

bp_doc_toc <- merge(bp_doc, bp_toc, by='date_ymd', all = TRUE) 
bp_doc_toc_cc <- bp_doc_toc %>% filter(!is.na(DOC) & !is.na(TOC_mg.L))

m1 <- lm(bp_doc_toc_cc$DOC ~ bp_doc_toc_cc$TOC_mg.L)
 summary(m1)
 
 # Residuals:
 #   Min      1Q  Median      3Q     Max 
 # -5.0226 -0.2867  0.0031  0.3585  1.9894 
 # 
 # Coefficients:
 #   Estimate Std. Error t value Pr(>|t|)    
 # (Intercept)             0.61572    0.05737   10.73   <2e-16 ***
 #   bp_doc_toc_cc$TOC_mg.L  0.81624    0.00766  106.56   <2e-16 ***
 #   ---
 #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 # 
 # Residual standard error: 0.6743 on 2107 degrees of freedom
 # Multiple R-squared:  0.8435,	Adjusted R-squared:  0.8434 
 # F-statistic: 1.135e+04 on 1 and 2107 DF,  p-value: < 2.2e-16
 # 
 # par(mfrow = c(2, 2)); plot(m1)

bp_doc_toc %>% 
  ggplot(aes(TOC_mg.L, DOC)) + 
  theme_bw(base_size = 12) +
  geom_point(col = "steelblue", alpha = 1/5, size = 3) + 
  geom_abline(slope = 1, intercept = c(0, 0), col = "grey", lty = 2) +
  geom_abline(intercept = m1$coefficients[1], 
              slope = m1$coefficients[2],
              lwd = 1, colour = "forestgreen") +
  labs(x = "TOC", y = "DOC")

bp_doc_toc_filled <- bp_doc_toc %>%
  mutate(note = ifelse(is.na(DOC) & is.na(TOC_mg.L), "Missing both DOC and TOC", 
                       ifelse(is.na(DOC) & !is.na(TOC_mg.L), "Modeled with TOC", "Measured")),
         DOC = ifelse(is.na(DOC), m1$coefficients[2]*TOC_mg.L + m1$coefficients[1], DOC)) 


bp_doc_toc_filled <- mutate(bp_doc_toc_filled,
             year = as.numeric(format(date_ymd,'%Y')))

bp_doc_toc_filled %>% 
  mutate(DOC_mg.L = ifelse(is.na(DOC), 0, DOC)) %>% 
  ggplot(aes(yday(date_ymd), DOC_mg.L, col = note)) + 
  facet_wrap(~ year) + 
  geom_point()

bp_doc_toc_filled_90s<- bp_doc_toc_filled %>%
  filter(year %in% c(1991:1993))%>%
  select(date_ymd, DOC)

# 1991-02-11 until 1993-06-07 modeled with TOC.

## add the DOC back into the dataset

# Perform a left join to merge missing DOC data
d2 <- d2 %>%
  left_join(bp_doc_toc_filled_90s, by = "date_ymd") %>%
  # If needed, update the DOC column with new values from missing_DOC_data
  mutate(DOC = coalesce(DOC.y, DOC.x)) %>%
  select(-DOC.x, -DOC.y)  # Remove extra columns created by the join

# see summary stats, have they changed sig?

parameter_summary <- d2 %>%
  dplyr::summarise(
    mean = mean(DOC, na.rm = TRUE),
    median = median(DOC, na.rm = TRUE),
    sd = sd(DOC, na.rm = TRUE),
    min = min(DOC, na.rm = TRUE),
    max = max(DOC, na.rm = TRUE),
    count = n()
  )

print(parameter_summary)


# # A tibble: 1 × 6
# mean median    sd   min   max count
# <dbl>  <dbl> <dbl> <dbl> <dbl> <int>
#  6.51    6.1  1.85  1.96    15  2609

d2 <- d2[rowSums(is.na(d2)) <= 6, ]  

d2 <- d2%>%
  dplyr::rename(
    Date = "date_ymd",
    NO3_mg.L = "Nitrate",
    NH3_mg.L = "Ammonia N",
    SRP_ug.L = "Phosphate (ortho)",
    TP_ug.L = "Phosphate (total)",
    Temp_C = "Temperature",
    Org_N_mg.L = "Organic N"
    )


## Add flow d## Add flow d## Add flow data

flow <- read_csv("Data/daily flow.csv")

flow<- flow %>%
  select(date, combined_05JG004.cms, SK05JG006.cms, RC_IC_cms)%>%
  dplyr::rename("Date" = date)


full_df<- merge(d2, flow, by="Date", all=TRUE)


full_df <- full_df[rowSums(is.na(full_df)) <= 7, ]       # Apply rowSums & is.na


full_df <- mutate(full_df,
                  DOY = as.numeric(format(Date,'%j')),
                  nMonth = as.numeric(format(Date,'%m')))

full_df <- full_df[rowSums(is.na(full_df)) <= 3, ] 
# add anthony imputations from 2003-2004 for SRP and ammonia

#impute from nearest neighbours

# # ammonia
# full_df %>% filter(year == 2002) # c(0.012, 0.07) for month 3
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


# # OrgN
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

#full_df %>% filter(year %in% 2004) # c(45, 94.29)



full_df <- full_df %>% 
  mutate(TP_ug.L = case_when(
    is.na(TP_ug.L) & year %in% c(2003:2004) ~ mean(c(45, 94.29)),
    TRUE ~ as.numeric(TP_ug.L)
  ))




## change 0s for N and P to limits of detection (need to find out for TP, TN)

sum(full_df$NH3_mg.L == 0, na.rm = TRUE) #289
sum(full_df$NO3_mg.L == 0, na.rm = TRUE) #512
sum(full_df$SRP_ug.L == 0, na.rm = TRUE) #177

full_df$NO3_mg.L[full_df$NO3_mg.L == 0] <- (0.057/2)
full_df$NH3_mg.L[full_df$NH3_mg.L == 0] <- (0.086/2)
full_df$SRP_ug.L[full_df$SRP_ug.L == 0] <- (3/2)



parameter_summary <- full_df %>%
  dplyr::summarise(
    mean = mean(DOC, na.rm = TRUE),
    median = median(DOC, na.rm = TRUE),
    sd = sd(DOC, na.rm = TRUE),
    min = min(DOC, na.rm = TRUE),
    max = max(DOC, na.rm = TRUE),
    count = n()
  )

print(parameter_summary)
# mean      median   sd  min max count
# 6.51      6.1   1.97  1.96  15  1233
# median, mean a bit less, SD a bit higher

write.csv(full_df, file="data/bpgamdataCLEAN_DOC.csv", row.names=F)



full_df %>% filter(is.na(NH3_mg.L)) #51 NAs
full_df %>% filter(is.na(NO3_mg.L)) #41 NAs
full_df %>% filter(is.na(SRP_ug.L)) #57 NAs
full_df %>% filter(is.na(DOC)) #114 NAs

