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


##### -------------------------DIN--------------------------- ###########


# DOC

DOC_DIN<- readRDS("DIN_DOC_notemp.rds")  # load the ggplot object
#print(DOC_DIN)


# Odour

Odour_DIN <- readRDS("DIN_Odour_notemp.rds")  # load the ggplot object
#print(Odour_DIN)


# Turb

TURB_DIN <- readRDS("DIN_TURB_notemp.rds")  # load the ggplot object
#print(TURB_DIN)

# TDS

TDS_DIN <- readRDS("DIN_TDS_notemp.rds")  # load the ggplot object
#print(TDS_DIN)

# 
# # cond
# 
# Conductivity_DIN <- readRDS("DIN_Conductivity_notemp.rds")  # load the ggplot object
# #print(Conductivity_DIN)

## DIN not significant for TDS 




# p_DIN<- plot_grid(TURB_DIN, Conductivity_DIN, Odour_DIN, DOC_DIN )
# p_DIN


#ggsave('output/DIN ALL.png', p_DIN, height = 8, width  = 10)




#####-------------------------SRP-------------------------- ###########

# DOC

DOC_SRP <- readRDS("SRP_DOC_notemp.rds")  # load the ggplot object
#print(DOC_SRP)


# Odour

Odour_SRP <- readRDS("SRP_Odour_notemp.rds")  # load the ggplot object
#print(Odour_SRP)


# TURB

TURB_SRP <- readRDS("SRP_TURB_notemp.rds")  # load the ggplot object
#print(TURB_SRP)


# # Conductivity
# 
# Conductivity_SRP <- readRDS("SRP_Conductivity_notemp.rds")  # load the ggplot object
# #print(Conductivity_SRP)

# 
# p_SRP<- plot_grid(TURB_SRP, Conductivity_SRP, Odour_SRP, DOC_SRP, align = "hv")
# p_SRP


#ggsave('output/SRP ALL.png', p_SRP, height = 8, width  = 10)


# 
# #####----------------------------- TEMP-----------------------------###########
# 
# # DOC
# 
# DOC_temp <- readRDS("temp_DOC_notemp.rds")  # load the ggplot object
# #print(DOC_temp)
# 
# # TDS
# 
# TDS_temp <- readRDS("temp_TDS_notemp.rds")  # load the ggplot object
# #print(TDS_temp)
# 
# # Odour
# 
# Odour_temp <- readRDS("temp_Odour_notemp.rds")  # load the ggplot object
# #print(Odour_temp)
# 
# 
# # TURB
# 
# TURB_temp <- readRDS("temp_TURB_notemp.rds")  # load the ggplot object
# #print(TURB_temp)
# 
# # Conductivity
# 
# Conductivity_temp <- readRDS("temp_Conductivity_notemp.rds")  # load the ggplot object
# #print(Conductivity_temp)
# 
# 
# p_TEMP<- plot_grid(TURB_temp, Conductivity_temp, Odour_temp, DOC_temp, align = "hv" )
# p_TEMP
# 
# 
# #ggsave('output/TEMP ALL.png', p_TEMP, height = 8, width  = 10)
# 



#####------------------------- QLD--------------------------------- ###########

# DOC

DOC_QLD <- readRDS("QLD_DOC_notemp.rds")  # load the ggplot object
#print(DOC_QLD)


# TDS

TDS_QLD <- readRDS("QLD_TDS_notemp.rds")  # load the ggplot object
#print(TDS_QLD)


# QLD NOT SIGNIFICANT FOR Odour


# TURB

TURB_QLD <- readRDS("QLD_TURB_notemp.rds")  # load the ggplot object
#print(TURB_QLD)

# TURB

#Conductivity_QLD <- readRDS("QLD_Conductivity_notemp.rds")  # load the ggplot object
#print(Conductivity_QLD)

# p_QLD<- plot_grid(TURB_QLD, Conductivity_QLD, DOC_QLD, align = "hv")
# p_QLD
# 

#ggsave('output/QLD ALL.png', p_QLD, height = 8, width  = 10)



##--- load in insig ----# 

# TDS_SRP

TDS_SRP <- readRDS("SRP_TDS_notemp.rds")  # load the ggplot object
#print(TDS_SRP)

# Odour_QLD

Odour_QLD <- readRDS("QLD_Odour_notemp.rds")  # load the ggplot object
#print(TDS_QLD)




p_all<- plot_grid(DOC_QLD, DOC_SRP, DOC_DIN,   
                  Odour_QLD, Odour_SRP,Odour_DIN,
                  TDS_QLD,TDS_SRP, TDS_DIN, 
                  TURB_QLD,TURB_SRP,TURB_DIN,
                  align = "hv", nrow= 4, ncol = 3)
p_all

ggsave('output/DW ALL_indicators preditced TDS same scale.png', p_all, height = 8, width  = 12, dpi = 300)

### with conductivity 

p_all<- plot_grid(DOC_DIN, DOC_SRP, DOC_QLD, DOC_temp,  
                  Odour_DIN,Odour_SRP,Odour_QLD,Odour_temp, 
                  Conductivity_DIN, Conductivity_SRP,Conductivity_QLD,Conductivity_temp, 
                  TURB_DIN,TURB_SRP,TURB_QLD,TURB_temp,
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

Odour_SOIPDO <- readRDS("SOIPDO_Odour.rds")  # load the ggplot object
print(Odour_SOIPDO)

# TURB

TURB_SOIPDO <- readRDS("SOIPDO_TURB_notemp.rds")  # load the ggplot object
print(TURB_SOIPDO)



p_SOIPDO <- (TURB_SOIPDO| TDS_SOIPDO ) /
  (Odour_SOIPDO| DOC_SOIPDO)
p_SOIPDO

ggsave('output/SOIPDO ALL predicted TDS_notemp_colours.png', p_SOIPDO, height = 11, width  = 18)



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


