###############################################################################
####### SOI and PDO Data  ##########
####### Danielle Spence ##########
####### Created 6/5/2024 ########
###############################################################################

rm(list = ls())

library(pacman)
p_load(tidyverse,  readr,  dplyr,  rsoi, install = TRUE)


setwd("C:/Users/danis/OneDrive/R/DW Indicators")

## Load in PDO, downloaded from: https://www.ncei.noaa.gov/access/monitoring/pdo/

pdo <- read.csv("data/ersst.v5.pdo.dat.txt", sep="")

pdo<- pdo%>%
  pivot_longer(c(Jan:Dec), names_to = "month", values_to = "PDO")

pdo$month<- mapvalues(pdo$month,from = c("Jan","Feb","Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), 
                      to = c(1,2,3,4,5,6,7,8,9,10,11,12))

pdo$Date = as.Date(paste0(pdo$Year, "-", pdo$month, "-01"))

# calculate 3-month average for PDO

pdo$PDO_3MON_AVG <- stats::filter(pdo$PDO, rep(1 / 3, 3), sides = 2)

#pdo$date_lag<- pdo$Date %m-% months(6)


pdo1<- pdo%>%
  select(Date, PDO, PDO_3MON_AVG)


pdo2<- do.call("rbind", lapply(1:nrow(pdo), function(i) 
  data.frame(date = seq(pdo1$Date[i], 
                        (seq(pdo1$Date[i],length=2,by="months") - 1)[2], by = "1 days"), 
             PDO_3MON_AVG = pdo1$PDO_3MON_AVG[i],
             pdo = pdo1$PDO[i])))



##----- add in ONI data (source is also NOAA), 

oni <- read.csv("data/Raw Data/oni.ascii.txt", sep="")

oni <- oni%>%
dplyr::rename(year = YR, 
         season = SEAS)

unique(oni$season)

#soi<- soi%>%
# select(Date, Year, SOI, SOI_3MON_AVG)

# Helper function to get months from 3-letter season string
season_month_map <- tibble::tibble(
  season = c("DJF", "JFM", "FMA", "MAM", "AMJ", "MJJ",
             "JJA", "JAS", "ASO", "SON", "OND", "NDJ"),
  months = list(c(12, 1, 2),   # DJF
                c(1, 2, 3),    # JFM
                c(2, 3, 4),    # FMA
                c(3, 4, 5),    # MAM
                c(4, 5, 6),    # AMJ
                c(5, 6, 7),    # MJJ
                c(6, 7, 8),    # JJA
                c(7, 8, 9),    # JAS
                c(8, 9, 10),   # ASO
                c(9, 10, 11),  # SON
                c(10, 11, 12), # OND
                c(11, 12, 1))  # NDJ
)

# Expand ONI values
oni_expanded <- oni %>%
  left_join(season_month_map, by = "season") %>%
  unnest(months) %>%
  dplyr::rename(month = months)

oni_expanded <- oni_expanded %>%
  mutate(
    adj_year = case_when(
      month == 12 & season %in% c("DJF") ~ year - 1,
      month == 1  & season == "NDJ" ~ year + 1,
      TRUE ~ year
    )
  ) %>%
  select(year = adj_year, month, ANOM) %>%
  arrange(year, month)

# expand to have all days to match with full dataset

oni_expanded$Date = as.Date(paste0(oni_expanded$year, "-", oni_expanded$month, "-01"))

oni1<- do.call("rbind", lapply(1:nrow(oni_expanded), function(i) 
  data.frame(date = seq(oni$Date[i], 
                        (seq(oni$Date[i],length=2,by="months") - 1)[2], by = "1 days"), 
             ONI_3MON_AVG = oni$ANOM[i])))


clim<- merge( pdo1,soi1, by="date" )

write.csv(clim, file="data/SOI_PDO data no lag.csv", row.names=F)



##----- add in SOI (source is also NOAA), lag 6 months

soi <- read.csv("data/soi.long.data.txt", sep="")

#soi<- soi%>%
# select(Date, Year, SOI, SOI_3MON_AVG)

soi<- soi%>%
  pivot_longer(c(Jan:Dec), names_to = "month", values_to = "SOI")

soi$month<- mapvalues(soi$month,from = c("Jan","Feb","Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), 
                      to = c(1,2,3,4,5,6,7,8,9,10,11,12))

soi$Date = as.Date(paste0(soi$Year, "-", soi$month, "-01"))

# calculate 3-month average for PDO

soi$SOI_3MON_AVG <- stats::filter(soi$SOI, rep(1 / 3, 3), sides = 2)


# repeat for each date for merging with full dataset

soi1<- do.call("rbind", lapply(1:nrow(soi), function(i) 
  data.frame(date = seq(soi$Date[i], 
                        (seq(soi$Date[i],length=2,by="months") - 1)[2], by = "1 days"), 
             SOI_3MON_AVG = soi$SOI_3MON_AVG[i],
             soi = soi$SOI[i])))


clim<- merge( pdo1,soi1, by="date" )

write.csv(clim, file="data/SOI_PDO data no lag.csv", row.names=F)


## comparing to climate trends....

clim1 <- mutate(clim,
             year = as.numeric(format(date,'%Y')))%>% 
  filter(year %in% c(1984:2022))

annP <-aggregate(clim1$PDO_3MON_AVG,  
               by=list(clim1$year),  
                FUN=mean) 
annP<- annP %>% dplyr::rename(PDO = x)

annS <-aggregate(clim1$SOI_3MON_AVG,  
                 by=list(clim1$year),  
                 FUN=mean) 
annS<- annS %>% dplyr::rename(SOI = x)

clim2<- merge(annP, annS, by="Group.1")
clim2<- clim2%>%
  dplyr::rename(Year = 'Group.1')


clim2<- clim2%>% 
  mutate(PDO_desc = ifelse(clim2$PDO < 0, "Wet", "Dry"),
         SOI_desc = ifelse(clim2$SOI > 0, "Cool_Wet", "Warm_Dry"))

clim2<- clim2%>% 
  mutate(Desc = case_when(
    PDO < 0 & SOI > 0 ~ "Cool_wet",  # Condition 1
    PDO > 0 & SOI < 0 ~ "Warm_dry",  # Condition 2
    TRUE ~ "Other"                         # Default case
  ))

write.csv(clim2, file="data/annual SOI PDO no lag.csv", row.names=F)
