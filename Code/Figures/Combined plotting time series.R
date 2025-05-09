###############################################################################
####### BPWTP Data -- GAMs ##########
####### Danielle Spence ##########
####### Created 5/4/2023 ########
###############################################################################
### Clear memory
rm(list = ls())

library(pacman)
p_load(tidyverse, ggplot2, mgcv, gratia, readr, GGally, dplyr,  mgcViz, lubridate,
       cowplot, tibble, patchwork, install = TRUE)

setwd("C:/Users/danis/OneDrive/R/DW indicators")

papertheme <- theme_bw(base_size=14, base_family = 'Arial') 


## Load in plots ###


##### -------------------------DOY--------------------------- ###########


# DOC

DOC_DOY<- readRDS("DOY_DOC_AC.rds")  # load the ggplot object
print(DOC_DOY)


# Odour

Odour_DOY <- readRDS("DOY_Odour_AC.rds")  # load the ggplot object
print(Odour_DOY)


# Turb

TURB_DOY <- readRDS("DOY_TURB_AC.rds")  # load the ggplot object
print(TURB_DOY)

# TDS

TDS_DOY <- readRDS("DOY_Predicted_TDS_AC.rds")  # load the ggplot object
print(TDS_DOY)

# 
# # cond
# 
# Conductivity_DOY <- readRDS("DOY_Conductivity_notemp.rds")  # load the ggplot object
# #print(Conductivity_DOY)

## DOY not significant for TDS 




# p_DOY<- plot_grid(TURB_DOY, Conductivity_DOY, Odour_DOY, DOC_DOY )
# p_DOY


#ggsave('output/DOY ALL.png', p_DOY, height = 8, width  = 10)




#####-------------------------year-------------------------- ###########

# DOC

DOC_year <- readRDS("year_DOC_AC.rds")  # load the ggplot object
print(DOC_year)


# Odour

Odour_year <- readRDS("year_Odour_AC.rds")  # load the ggplot object
print(Odour_year)


# TURB

TURB_year <- readRDS("year_TURB_AC.rds")  # load the ggplot object
print(TURB_year)


# TDS
TDS_year <- readRDS("year_Predicted_TDS_AC.rds")  # load the ggplot object
print(TDS_year)

# 
# p_year<- plot_grid(TURB_year, Conductivity_year, Odour_year, DOC_year, align = "hv")
# p_year


#ggsave('output/year ALL.png', p_year, height = 8, width  = 10)


#####------------------------- time--------------------------------- ###########

# DOC

DOC_time <- readRDS("time_DOC.rds")  # load the ggplot object
print(DOC_time)


# TDS

TDS_time <- readRDS("time_Predicted_TDS.rds")  # load the ggplot object
print(TDS_time)



# TURB

TURB_time <- readRDS("time_TURB.rds")  # load the ggplot object
print(TURB_time)


Odour_time <- readRDS("time_Odour.rds")  # load the ggplot object
print(Odour_time)


# p_time<- plot_grid(TURB_time, Conductivity_time, DOC_time, align = "hv")
# p_time
# 

#ggsave('output/time ALL.png', p_time, height = 8, width  = 10)

p_all<-plot_grid(DOC_DOY,  DOC_year,Odour_DOY , Odour_year, TDS_DOY,TDS_year,TURB_DOY ,TURB_year,
          ncol = 2, align = "hv")

p_all

ggsave('output/DW time series_AC.png', p_all, height = 12, width  = 10, dpi = 300)

### with conductivity 

p_all<- plot_grid(DOC_DOY, DOC_year, DOC_time, DOC_temp,  
                  Odour_DOY,Odour_year,Odour_time,Odour_temp, 
                  Conductivity_DOY, Conductivity_year,Conductivity_time,Conductivity_temp, 
                  TURB_DOY,TURB_year,TURB_time,TURB_temp,
                  align = "hv", nrow= 4, ncol = 4)
p_all

ggsave('output/DW ALL_indicators Conductivity_notemp.png', p_all, height = 8, width  = 12)

##### --------------------------SOI_PDO ----------------------------###########

# DOC

DOC_SOIPDO <- readRDS("SOIPDO_DOC_notemp.rds")  # load the ggplot object
print(DOC_SOIPDO)

# TDS

TDS_SOIPDO <- readRDS("SOIPDO_TDS_notemp.rds")  # load the ggplot object
print(TDS_SOIPDO)


# Conductivity

#Conductivity_SOIPDO <- readRDS("SOIPDO_Conductivity_notemp.rds")  # load the ggplot object


# Odour

Odour_SOIPDO <- readRDS("SOIPDO_Odour_notemp.rds")  # load the ggplot object
print(Odour_SOIPDO)

# TURB

TURB_SOIPDO <- readRDS("SOIPDO_TURB_notemp.rds")  # load the ggplot object
print(TURB_SOIPDO)



p_SOIPDO <- (TURB_SOIPDO| TDS_SOIPDO ) /
  (Odour_SOIPDO| DOC_SOIPDO)
p_SOIPDO

ggsave('output/SOIPDO ALL predicted TDS_notemp.png', p_SOIPDO, height = 10, width  = 14)



##### ---------------------------TIME ---------------------------###########


# DOC

DOC_time <- readRDS("time_DOC_temp.rds")  # load the ggplot object
#print(DOC_time)

# TDS

TDS_time <- readRDS("time_TDS_notemp.rds")  # load the ggplot object
#print(TDS_time)


# Odour

Odour_time <- readRDS("time_Odour_notemp.rds")  # load the ggplot object
##print(Odour_time)


# TURB

TURB_time <- readRDS("time_TURB.rds")  # load the ggplot object
#print(TURB_time)


# conductivity
#Conductivity_time <- readRDS("time_Conductivity_notemp.rds")  # load the ggplot object
#print(Conductivity_time) 



p_all<- (TURB_time | TDS_time)/
  (Odour_time | DOC_time)

p_all


ggsave('output/TIME ALL TDS predicted_notemp.png', p_all, height = 10, width  = 14)


