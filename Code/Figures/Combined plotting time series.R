###############################################################################
####### BPWTP Data -- GAMs ##########
####### Danielle Spence ##########
####### Created 5/4/2023 ########
###############################################################################
### Clear memory
rm(list = ls())

library(pacman)
p_load(tidyverse, ggplot2, mgcv, gratia, readr, GGally, dplyr,  mgcViz, lubridate,
       cowplot, tibble, patchwork, lehuynh,install = TRUE)

setwd("C:/Users/danis/OneDrive/R/DW indicators")

papertheme <- theme_bw(base_size=14, base_family = 'Arial') 


## Load in plots ###


##### -------------------------DOY--------------------------- ###########

# chla

chla_DOY<- readRDS("DOY_chla.rds")  # load the ggplot object
print(chla_DOY)

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


# p_DOY<- plot_grid(TURB_DOY, Conductivity_DOY, Odour_DOY, DOC_DOY )
# p_DOY


#ggsave('output/DOY ALL.png', p_DOY, height = 8, width  = 10)


#####-------------------------year-------------------------- ###########

# chla

chla_year <- readRDS("year_chla.rds")  # load the ggplot object
print(chla_year)

# DOC

DOC_year <- readRDS("year_DOC_AC.rds")  # load the ggplot object
print(DOC_year)

ggsave_elsevier("output/DOC year.jpeg", plot = DOC_year, width = "one_column", height =  90)

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



p_all<-plot_grid(DOC_DOY,  DOC_year,Odour_DOY , Odour_year, TDS_DOY,TDS_year,TURB_DOY ,TURB_year,
          ncol = 3, align = "hv")

p_all


ggsave_elsevier("output/Figure_2.jpeg", plot = p_all, width = "full_page", height =  240)



##------------------------------------ Model Residuals--------------------------#

TDS_res <- readRDS("TDS_AC_residuals.rds")  # load the ggplot object
print(TDS_res)

DOC_res <- readRDS("DOC_AC_residuals.rds")  # load the ggplot object
print(DOC_res)

turb_res <- readRDS("turb_AC_residuals.rds")  # load the ggplot object
print(turb_res)

Odour_res <- readRDS("Odour_AC_residuals.rds")  # load the ggplot object
print(Odour_res)



p_all_res<-plot_grid(DOC_res,  Odour_res,TDS_res,turb_res,
                 ncol = 2, align = "hv", labels = c("A", "B", "C", "D"), label_size = 16, label_x = 0.01, label_y = 1)

p_all_res

ggsave('output/residuals AC timeseries.jpeg', p_all_res, height = 12, width  = 15, dpi = 300)

