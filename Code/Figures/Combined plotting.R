###############################################################################
####### BPWTP Data -- GAMs ##########
####### Danielle Spence ##########
####### Created 5/4/2023 ########
###############################################################################
### Clear memory
rm(list = ls())

library(pacman)
p_load(tidyverse, ggplot2, mgcv, gratia, readr, GGally, dplyr,  mgcViz, lubridate,
       cowplot, tibble, patchwork, lehuynh, install = TRUE)

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

# TDS_SRP

TDS_SRP <- readRDS("SRP_TDS_notemp.rds")  # load the ggplot object
#print(TDS_SRP)


# p_SRP<- plot_grid(TURB_SRP, Conductivity_SRP, Odour_SRP, DOC_SRP, align = "hv")
# p_SRP


#ggsave('output/SRP ALL.png', p_SRP, height = 8, width  = 10)



#####------------------------- QLD--------------------------------- ###########

# DOC

DOC_QLD <- readRDS("QLD_DOC_notemp.rds")  # load the ggplot object
#print(DOC_QLD)


# TDS

TDS_QLD <- readRDS("QLD_TDS_notemp.rds")  # load the ggplot object
#print(TDS_QLD)

# Odour_QLD

Odour_QLD <- readRDS("QLD_Odour_notemp.rds")  # load the ggplot object
#print(TDS_QLD)

# TURB

TURB_QLD <- readRDS("QLD_TURB_notemp.rds")  # load the ggplot object
#print(TURB_QLD)


# p_QLD<- plot_grid(TURB_QLD, Conductivity_QLD, DOC_QLD, align = "hv")
# p_QLD
# 

#ggsave('output/QLD ALL.png', p_QLD, height = 8, width  = 10)



##--- plot it all together now ----# 

p_all<- plot_grid(TURB_QLD,TURB_SRP,TURB_DIN,
                  TDS_QLD,TDS_SRP, TDS_DIN,
                  Odour_QLD, Odour_SRP,Odour_DIN,
                  DOC_QLD, DOC_SRP, DOC_DIN,
                  align = "hv", nrow= 4, ncol = 3)
p_all


ggsave_elsevier("output/Figure_3.jpeg", plot = p_all, width = "full_page", height =  190)


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

#ggsave('output/Figure_4.jpeg', p_SOIPDO, height = 8, width  = 12, dpi = 1000)


ggsave_elsevier("output/Figure_4.jpeg", plot = p_SOIPDO, width = "full_page", height =  190)


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



p_all<- (TURB_time | TDS_time)/
  (Odour_time | DOC_time)

p_all


ggsave('output/TIME ALL TDS predicted_notemp.jpeg', p_all, height = 10, width  = 14)


##### ---------------------------Model diagnostics ---------------------------###########

## model diagnostics
TDS_res <- readRDS("TDS_residuals.rds")  # load the ggplot object
print(TDS_res)

DOC_res <- readRDS("DOC_residuals.rds")  # load the ggplot object
print(DOC_res)

turb_res <- readRDS("turb_residuals.rds")  # load the ggplot object
print(turb_res)

Odour_res <- readRDS("Odour_residuals.rds")  # load the ggplot object
print(Odour_res)



p_all_res<-plot_grid(DOC_res,  Odour_res,TDS_res,turb_res,
                     ncol = 2, align = "hv", labels = c("A", "B", "C", "D"), label_size = 16, label_x = 0.01, label_y = 1)

p_all_res

ggsave('output/residuals full timeseries.jpeg', p_all_res, height = 12, width  = 15, dpi = 300)
