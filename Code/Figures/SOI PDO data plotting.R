###############################################################################
####### BPWTP Data -- GAMs ##########
####### Danielle Spence ##########
####### Created 5/4/2023 ########
###############################################################################
### Clear memory
rm(list = ls())

library(pacman)
p_load(tidyverse, ggplot2, readr, GGally, tidyr,  lubridate, 
       cowplot, tibble, viridis,  zoo, install = TRUE)

setwd("C:/Users/danis/OneDrive/R/DW indicators") 


##----------------------------------------------------------------------------##
## 1. Read in data 
##----------------------------------------------------------------------------##


df <- read_csv("data/SOI_PDO data no lag.csv")

df <- mutate(df,
             year = as.numeric(format(date,'%Y')))%>% 
  filter(year %in% c(1990:2022))

df_annual <- aggregate(df$SOI_3MON_AVG,  
                       by=list(df$year),  
                       FUN=mean,
                       na.rm=TRUE) 


df_annual<- df_annual %>% dplyr::rename(soi = x)

# Compute moving average
df_annual$moving_avg <- rollmean(df_annual$soi, k = 5, fill = NA, align = "right")

p1<- ggplot(df_annual, aes(x = Group.1, y = soi)) +
  geom_col(aes(fill = soi > 0), show.legend = FALSE) +  
  scale_y_continuous(limits=c(-2, 1.5))+
  geom_line(aes(y = moving_avg), color = "black", size = 1.2) +  # Moving average trend
  scale_fill_manual(values = c("red", "blue")) +  
  theme_bw() +
  labs(x = "Year", y = "Average SOI")
p1

# PDO
df_annualP <- aggregate(df$PDO_3MON_AVG,  
                       by=list(df$year),  
                       FUN=mean,
                       na.rm=TRUE) 


df_annualP<- df_annualP %>% dplyr::rename(pdo = x)

# Compute moving average
df_annualP$moving_avg <- rollmean(df_annualP$pdo, k = 5, fill = NA, align = "right")

p2<- ggplot(df_annualP, aes(x = Group.1, y = pdo)) +
  geom_col(aes(fill = pdo > 0), show.legend = FALSE) +  
  geom_line(aes(y = moving_avg), color = "black", size = 1.2) +  # Moving average trend
  scale_y_continuous(limits=c(-2.5, 1.5))+
  scale_fill_manual(values = c("blue", "red")) +  
  theme_bw() +
  labs(x = "Year", y = "Average PDO")
p2



p_all<- plot_grid(p1, p2)
p_all


ggsave('output/SOI PDO averages.png', p_all, height = 6, width  = 8)


