 ###############################################################################
####### BPWTP Data -- GAMs ##########
####### Danielle Spence ##########
####### Created 5/4/2023 ########
###############################################################################
### Clear memory
rm(list = ls())

library(pacman)
p_load(tidyverse, ggplot2, mgcv, gratia, readr, GGally, dplyr, lubridate, dplyr,
        cowplot, tibble,viridis, install = TRUE)

setwd("C:/Users/danis/OneDrive/R/DW Indicators")

papertheme <- theme_bw(base_size=12, base_family = 'Arial') 

##----------------------------------------------------------------------------##
## 1. Read in data 
##----------------------------------------------------------------------------##

df <- read_csv("data/bpgamdataCLEAN_DOC.csv")
# 1353 obs

df$Date<- as.Date(df$Date, "%m/%d/%Y")

# add in Year and nMonth for numeric month and a proper Date class

df <- mutate(df,
             year = as.numeric(format(Date,'%Y')),
             DOY = as.numeric(format(Date,'%j')),
             nMonth = as.numeric(format(Date,'%m')),
             week = as.numeric(format(Date, '%W')))%>% 
  filter(year %in% c(1990:2022))

# 1670 obs

## order the data
df <- df[with(df, order(Date)),]


# rename for ease of use
df <- df%>%
  dplyr::rename(
    date = "Date",
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


#Remove rows without DOC data
df <- df[complete.cases(df$DOC),]
# 913 obs

df$DIN <- rowSums(df[,c("NO3", "NH3")], na.rm=TRUE)

## add in PDO and SOI 

clim <- read_csv("data/SOI_PDO data no lag.csv")

df1<- merge(df, clim, by="date")




##----------------------------------------------------------------------------##
## 3. Estimate annual and seasonal variability
##----------------------------------------------------------------------------##
# 
df_clean <- df1 %>%
  dplyr::ungroup() %>%                    # Remove any previous grouping
  select(year, DOC)


library(dplyr)

df_clean %>%
  group_by(year) %>%
  summarise(
    mean_year = mean(DOC, na.rm = TRUE),
    seasonal_range = max(DOC, na.rm = TRUE) - min(DOC, na.rm = TRUE)
  ) %>%
  summarise(
    avg_seasonal_range = mean(seasonal_range, na.rm = TRUE),
    avg_seasonal_range_pct = mean(seasonal_range / mean_year * 100, na.rm = TRUE),
    between_year_range = max(mean_year) - min(mean_year),
    between_year_range_pct = (max(mean_year) - min(mean_year)) / mean(mean_year) * 100
  )




##----------------------------------------------------------------------------##
## 5. Specify the best fit model (see Fitting script)
##----------------------------------------------------------------------------##


m2 <- gam(DOC ~ s(SRP) +  
            s(DIN)+ 
            s(QLD)+
            te(PDO_3MON_AVG,SOI_3MON_AVG)+
            te(year, DOY, bs = c("cr", "cc")),
          knots=list(DOY=c(0, 366.5)),
          select = TRUE,
          data = df1, method = "REML", family = tw(link="log"))


saveRDS(m2, file = "modelDOC.rds")

# save model outputs as text
#  sink("model_summaryDOC.txt")
#  summary(m2)
#  sink()  


# check model fit...

summary(m2)


k.check(m2)

p1<- draw(m,  residuals = TRUE)& theme_bw() 
p1

#ggsave('output/DOC GAM.png', p1, height = 10, width  = 10)

##----------------------------------------------------------------------------##
## 6. Assess model fit & autocorrelation
##----------------------------------------------------------------------------##

## i. Appraise model fit
p2<-appraise(m2, point_col = 'steelblue', point_alpha = 0.5, n_bins = 'fd') & 
  theme(plot.tag = element_text(face = 'bold')) & theme_bw()
p2

saveRDS(p2, "DOC_residuals.rds")

#ggsave('output/DOC GAM residuals no temp.png', p2, height = 10, width  = 10)


# ii. Check for autocorrelation
layout(matrix(1:2, ncol = 2))
acf(resid(m2), lag.max = 36, main = "ACF") 
pacf(resid(m2), lag.max = 36, main = "pACF")
layout(1)


#----------------------------------------------------------------------------##
## 7. Plotting on the response scale/fitted values
##----------------------------------------------------------------------------##

## ---------- DIN -------------------##

new_data_DIN<- with(df1, expand.grid(DIN = seq(min(DIN), max(DIN),length = 200),
                                           SRP = median(SRP, na.rm=TRUE),
                                           QLD = median(QLD),
                                           SOI_3MON_AVG = median(SOI_3MON_AVG, na.rm=TRUE),
                                           PDO_3MON_AVG = median(PDO_3MON_AVG, na.rm=TRUE),
                                           year= median(year),
                                           DOY = median(DOY)))

DIN.pred <- predict(m2, newdata = new_data_DIN, type = "terms", se.fit = TRUE)

whichCols <- grep("DIN", colnames(DIN.pred$fit))
whichColsSE <- grep("DIN", colnames(DIN.pred$se.fit))
new_data_DIN <- cbind(new_data_DIN, Fitted = DIN.pred$fit[, whichCols], 
                           se.Fitted = DIN.pred$se.fit[, whichColsSE])
limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

## make into original limits
new_data_DIN <- with(new_data_DIN, transform(new_data_DIN, Fittedplus = Fitted + se.Fitted))
new_data_DIN <- with(new_data_DIN, transform(new_data_DIN, Fittedminus = Fitted - se.Fitted))

shiftDIN <- attr(predict(m2, newdata = new_data_DIN, type = "iterms"), "constant")
DIN.pdatnorm <- new_data_DIN
DIN.pdatnorm <- with(DIN.pdatnorm, transform(DIN.pdatnorm, Fitted = Fitted + shiftDIN, 
                                                       Fittedplus = Fittedplus + shiftDIN, 
                                                       Fittedminus = Fittedminus + shiftDIN))


DIN.pdatnorm$Fitted<-exp(DIN.pdatnorm$Fitted)
DIN.pdatnorm$Fittedplus<-exp(DIN.pdatnorm$Fittedplus)
DIN.pdatnorm$Fittedminus<-exp(DIN.pdatnorm$Fittedminus)


DINquants <- quantile(df1$DIN, c(.05,.95), na.rm = TRUE)

DINplot <- ggplot(DIN.pdatnorm, aes(x = DIN, y = Fitted)) +
  papertheme +
  scale_y_continuous(limits=c(3, 15))+
  annotate("rect", xmin=DINquants[1], xmax=DINquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
  annotate("text", 
           x = -Inf, y = Inf,                  # Top-left corner
           hjust = -0.1, vjust = 1.5,          # Adjust alignment to keep text inside plot
           label = expression(italic("p")*" = 0.0797"))+
  geom_line() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25, fill = '#165459B2') +  
  geom_rug(aes(x=DIN), data = df1, stat = "identity", position = "identity", 
           sides = "b", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
  xlab(expression(paste("DIN"~"(mg/L)"))) + ylab(expression(paste("DOC"~"(mg/L)")))

DINplot


saveRDS(DINplot, "DIN_DOC_notemp.rds")  # save the ggplot object


# range of control
max(DIN.pdatnorm$Fitted)-min(DIN.pdatnorm$Fitted)
# = NONSIG

### ---------- SRP -------------------##

new_data_SRP <- with(df1, expand.grid(SRP = seq(min(SRP, na.rm = TRUE), max(SRP, na.rm = TRUE), length = 200),
                                     DIN = median(DIN),
                                     QLD = median(QLD),
                                     SOI_3MON_AVG = median(W_temp, na.rm=TRUE), 
                                     PDO_3MON_AVG = median(PDO_3MON_AVG, na.rm=TRUE),
                                     year= median(year),
                                     DOY = median(DOY)))


SRP.pred <- predict(m2, newdata = new_data_SRP, type = "terms", se.fit = TRUE)



whichCols <- grep("SRP", colnames(SRP.pred$fit))
whichColsSE <- grep("SRP", colnames(SRP.pred$se.fit))
new_data_SRP <- cbind(new_data_SRP, Fitted = SRP.pred$fit[, whichCols], 
                      se.Fitted = SRP.pred$se.fit[, whichColsSE])


limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

## make into original limits
new_data_SRP <- with(new_data_SRP, transform(new_data_SRP, Fittedplus = Fitted + se.Fitted))
new_data_SRP <- with(new_data_SRP, transform(new_data_SRP, Fittedminus = Fitted - se.Fitted))

shiftSRP <- attr(predict(m2, newdata = new_data_SRP, type = "iterms"), "constant")
SRP.pdatnorm <- new_data_SRP
SRP.pdatnorm <- with(SRP.pdatnorm, transform(SRP.pdatnorm, Fitted = Fitted + shiftSRP, 
                                             Fittedplus = Fittedplus + shiftSRP, 
                                             Fittedminus = Fittedminus + shiftSRP))

SRP.pdatnorm$Fitted<-exp(SRP.pdatnorm$Fitted)
SRP.pdatnorm$Fittedplus<-exp(SRP.pdatnorm$Fittedplus)
SRP.pdatnorm$Fittedminus<-exp(SRP.pdatnorm$Fittedminus)


SRPquants <- quantile(df1$SRP, c(.05,.95), na.rm = TRUE)


SRPplot <- ggplot(SRP.pdatnorm, aes(x = SRP, y = Fitted)) +
  papertheme+
  scale_y_continuous(limits=c(3, 15))+
  annotate("rect", xmin=SRPquants[1], xmax=SRPquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
  annotate("text", 
           x = -Inf, y = Inf,                  # Top-left corner
           hjust = -0.1, vjust = 1.5,          # Adjust alignment to keep text inside plot
           label = expression(italic("p")*" < 0.0001***"))+
  geom_line() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25, fill = '#165459B2') +  
  geom_rug(aes(x=SRP), data = df1, stat = "identity", position = "identity", 
           sides = "b", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
  xlab(expression(paste("SRP ","(", mu, "g/L)"))) + ylab(expression(paste("DOC"~"(mg/L)")))

SRPplot



saveRDS(SRPplot, "SRP_DOC_notemp.rds")  # save the ggplot objec


# range of effect

max_idx <- which.max(SRP.pdatnorm$Fitted)
min_idx <- which.min(SRP.pdatnorm$Fitted)
range_effect_SRP <- SRP.pdatnorm$Fitted[max_idx] - SRP.pdatnorm$Fitted[min_idx]
# 73.3

#lower range of CI
range_lwr <- SRP.pdatnorm$Fittedminus[max_idx] - SRP.pdatnorm$Fittedplus[min_idx]
range_upr <- SRP.pdatnorm$Fittedplus[max_idx] - SRP.pdatnorm$Fittedminus[min_idx]


# compare to seasonal and interannual variability

# Intra-annual variation (average seasonal swing)
within_year_var <- df1 %>%
  group_by(year) %>%
  summarise(range_within = max(DOC) - min(DOC), .groups = "drop") %>%
  summarise(avg_within = mean(range_within))

# Interannual variation (range of annual means)
between_year_var <- df1 %>%
  group_by(year) %>%
  summarise(mean_year = mean(DOC), .groups = "drop") %>%
  summarise(range_between = max(mean_year) - min(mean_year))



# relative_to_within = 
range_effect_SRP / within_year_var$avg_within

#1.24

#relative_to_between = 
range_effect_SRP / between_year_var$range_between
#0.0.64


## ---------- QLD -------------------##

new_data_QLD <- with(df1, expand.grid(QLD = seq(min(QLD), max(QLD), length = 200),
                                      DIN = median(DIN),
                                      SRP = median(SRP, na.rm=TRUE),
                                      SOI_3MON_AVG = median(W_temp, na.rm=TRUE), 
                                      PDO_3MON_AVG = median(PDO_3MON_AVG, na.rm=TRUE),
                                      year= median(year),
                                      DOY = median(DOY)))
QLD.pred <- predict(m2, newdata = new_data_QLD, type = "terms", se.fit = TRUE)

whichCols <- grep("QLD", colnames(QLD.pred$fit))
whichColsSE <- grep("QLD", colnames(QLD.pred$se.fit))
new_data_QLD <- cbind(new_data_QLD, Fitted = QLD.pred$fit[, whichCols], 
                       se.Fitted = QLD.pred$se.fit[, whichColsSE])
limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

## make into original limits
new_data_QLD <- with(new_data_QLD, transform(new_data_QLD, Fittedplus = Fitted + se.Fitted))
new_data_QLD <- with(new_data_QLD, transform(new_data_QLD, Fittedminus = Fitted - se.Fitted))

shiftQLD <- attr(predict(m2, newdata = new_data_QLD, type = "iterms"), "constant")
QLD.pdatnorm <- new_data_QLD
QLD.pdatnorm <- with(QLD.pdatnorm, transform(QLD.pdatnorm, Fitted = Fitted + shiftQLD, 
                                               Fittedplus = Fittedplus + shiftQLD, 
                                               Fittedminus = Fittedminus + shiftQLD))



QLD.pdatnorm$Fitted<-exp(QLD.pdatnorm$Fitted)
QLD.pdatnorm$Fittedplus<-exp(QLD.pdatnorm$Fittedplus)
QLD.pdatnorm$Fittedminus<-exp(QLD.pdatnorm$Fittedminus)


QLDquants <- quantile(df1$QLD, c(.05,.95), na.rm = TRUE)

QLDplot <- ggplot(QLD.pdatnorm, aes(x = QLD, y = Fitted)) +
  papertheme +
  scale_y_continuous(limits=c(3, 15))+
  annotate("rect", xmin=QLDquants[1], xmax=QLDquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
  annotate("text",
                    x = -Inf, y = Inf,                 
                    hjust = -0.1, vjust = 1.5,          
             label = expression(italic("p")*" = 0.000322***"))+
  geom_line() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25, fill = '#165459B2') +  
 geom_rug(aes(x=QLD), data = df1, stat = "identity", position = "identity", 
           sides = "b", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
  xlab(expression(paste("QLD (m/s)"))) + ylab(expression(paste("DOC"~"(mg/L)")))

QLDplot


saveRDS(QLDplot, "QLD_DOC_notemp.rds")  # save the ggplot objec

# range of effect

max_idx <- which.max(QLD.pdatnorm$Fitted)
min_idx <- which.min(QLD.pdatnorm$Fitted)
range_effect_QLD <- QLD.pdatnorm$Fitted[max_idx] - QLD.pdatnorm$Fitted[min_idx]
# 0.76

#lower range of CI
range_lwr <- QLD.pdatnorm$Fittedminus[max_idx] - QLD.pdatnorm$Fittedplus[min_idx]
range_upr <- QLD.pdatnorm$Fittedplus[max_idx] - QLD.pdatnorm$Fittedminus[min_idx]



# relative_to_within = 
range_effect_QLD / within_year_var$avg_within

#0.23

#relative_to_between = 
range_effect_QLD / between_year_var$range_between
#0.118


# PUT IT ALL TOGETHER

p_all<- plot_grid(DINplot, SRPplot, QLDplot, ncol = 2)
p_all


#ggsave('output/DOC GAM_Response scale.png', p_all, height = 10, width  = 10)



## ---------- SOI * PDO -------------------##



new_data_SOI_3MON_AVG <- with(df1, expand.grid(SOI_3MON_AVG = seq(min(SOI_3MON_AVG), max(SOI_3MON_AVG), length = 200),
                                               PDO_3MON_AVG = seq(min(PDO_3MON_AVG), max(PDO_3MON_AVG), length = 200),
                                               SRP = median(SRP, na.rm=TRUE),
                                               DIN = median(DIN, na.rm=TRUE),
                                               QLD = median(QLD, na.rm=TRUE),
                                               year= median(year),
                                               DOY = median(DOY)))

SOI_3MON_AVG.pred <- predict(m2, newdata = new_data_SOI_3MON_AVG, type = "terms")

whichCols <- grep("PDO_3MON_AVG,SOI_3MON_AVG", colnames(SOI_3MON_AVG.pred))

new_data_SOI_3MON_AVG <- cbind(new_data_SOI_3MON_AVG, Fitted = SOI_3MON_AVG.pred[, whichCols])

shiftcomb <- attr(SOI_3MON_AVG.pred, "constant")
SOI_3MON_AVG.pdatnorm <- new_data_SOI_3MON_AVG
SOI_3MON_AVG.pdatnorm <- with(SOI_3MON_AVG.pdatnorm, transform(SOI_3MON_AVG.pdatnorm, Fitted = Fitted + shiftcomb))

toofar <- exclude.too.far(SOI_3MON_AVG.pdatnorm$PDO_3MON_AVG, SOI_3MON_AVG.pdatnorm$SOI_3MON_AVG, df1$PDO_3MON_AVG, df1$SOI_3MON_AVG, dist=0.1)
SOI_3MON_AVG.pdatnorm$DOC <- SOI_3MON_AVG.pdatnorm$Fitted
SOI_3MON_AVG.pdatnorm$DOC[toofar] <- NA


SOI_3MON_AVG.pdatnorm$DOC<-exp(SOI_3MON_AVG.pdatnorm$DOC)

names(new_data_SOI_3MON_AVG)[which(names(new_data_SOI_3MON_AVG)=='SOI_3MON_AVG.pred')] <- 'DOC'

comboplot <- ggplot(SOI_3MON_AVG.pdatnorm, aes(x = PDO_3MON_AVG, y = SOI_3MON_AVG, z=DOC)) + 
  theme_bw(base_family = 'Arial') +
  theme(legend.position='top') +
  geom_raster(aes(fill=DOC)) + 
  scale_fill_distiller(palette = "Spectral", direction = -1, na.value='transparent') +
  geom_point(data=df1, aes(x=PDO_3MON_AVG, y=SOI_3MON_AVG, z=NULL)) +
  geom_contour(colour = "black", binwidth = 0.5) +
  scale_x_continuous(expand = c(0, 0)) +  # Remove padding on x-axis
  scale_y_continuous(expand = c(0, 0)) +  # Remove padding on y-axis
  theme(legend.key.width=unit(1,"cm"))+
  xlab("PDO") + ylab("SOI") +
  labs(fill=expression(paste("DOC"~"(mg/L)")))+
  theme(
    axis.text.x = element_text(color = ifelse(seq(-2.5, 2.5, by = 1) < -0.5, "blue",
                                              ifelse(seq(-2.5, 2.5, by = 1) > 0, "red", "black"))),
    axis.text.y = element_text(color = ifelse(seq(-2.5, 2.5, by = 1) < -0.5, "red",
                                              ifelse(seq(-2.5, 2.5, by = 1) > 0, "blue", "black")))
  )+
  theme(legend.position = "top")

comboplot

#ggsave('output/SOI PDO INTERACTION Response scale DOC.png', comboplot, height = 8, width  = 10)


saveRDS(comboplot, "SOIPDO_DOC_notemp.rds")  # save the ggplot objec

# range of effect

SOI_3MON_AVG.pdatnorm1<- na.omit(SOI_3MON_AVG.pdatnorm)

max_idx <- which.max(SOI_3MON_AVG.pdatnorm1$DOC)
min_idx <- which.min(SOI_3MON_AVG.pdatnorm1$DOC)
range_effect_SOI <- SOI_3MON_AVG.pdatnorm1$DOC[max_idx] - SOI_3MON_AVG.pdatnorm1$DOC[min_idx]
# 119.7

#lower range of CI
range_lwr <- SOI_3MON_AVG.pdatnorm1$Fittedminus[max_idx] - SOI_3MON_AVG.pdatnorm1$Fittedplus[min_idx]
range_upr <- SOI_3MON_AVG.pdatnorm1$Fittedplus[max_idx] - SOI_3MON_AVG.pdatnorm1$Fittedminus[min_idx]


# relative_to_within = 
range_effect_SOI / within_year_var$avg_within

#1.06

#relative_to_between = 
range_effect_SOI / between_year_var$range_between
#0.55


##-----------------------------------temporal change------------------------------- ####


new_data_year <- with(df1, expand.grid(DOY = seq(min(DOY), max(DOY),length = 200),
                                      year = seq(min(year), max(year),length = 200),
                                      SOI_3MON_AVG = median(SOI_3MON_AVG),
                                      PDO_3MON_AVG = median(PDO_3MON_AVG),
                                      SRP = median(SRP, na.rm=TRUE),
                                      DIN = median(DIN, na.rm=TRUE),
                                      QLD = median(QLD, na.rm=TRUE)))

year.pred <- predict(m2, newdata = new_data_year, type = "terms")

whichCols <- grep("year,DOY", colnames(year.pred))

new_data_year <- cbind(new_data_year, Fitted = year.pred[, whichCols])

shiftcomb <- attr(year.pred, "constant")
year.pdatnorm <- new_data_year
year.pdatnorm <- with(year.pdatnorm, transform(year.pdatnorm, Fitted = Fitted + shiftcomb))

toofar <- exclude.too.far(year.pdatnorm$DOY, year.pdatnorm$year, df1$DOY, df1$year, dist=0.1)
year.pdatnorm$DOC <- year.pdatnorm$Fitted
year.pdatnorm$DOC[toofar] <- NA


year.pdatnorm$DOC<-exp(year.pdatnorm$DOC)


names(new_data_year)[which(names(new_data_year)=='year.pred')] <- 'DOC'

comboplot_time <- ggplot(year.pdatnorm, aes(x = year, y = DOY, z=DOC)) + 
  theme_minimal(base_size=14, base_family = 'Arial') +
  theme(legend.position='top') +
  geom_raster(aes(fill=DOC)) + 
  scale_fill_distiller(palette = "Spectral", direction = -1, na.value='transparent') +
  geom_point(data=df1, aes(x=year, y=DOY, z=NULL)) +
  geom_contour(colour = "black", binwidth = 1) +
  theme(legend.key.width=unit(1.5,"cm")) +
  scale_y_continuous(expand = c(0, 0), breaks = c(1, 91,  182, 274, 335), labels = c("Jan", "Apr", "Jul", "Oct", "Dec")) +
  scale_x_continuous(expand = c(0, 0)) +  # Remove padding on x-axis
  xlab("Year") + ylab("DOY") +
  labs(fill=expression(paste("DOC (mg/L)   ")))+
  theme(legend.position = "top")


comboplot_time

#ggsave('output/DOC change over time_nolag.png', comboplot_time, height = 6, width  = 8)


saveRDS(comboplot_time, "time_DOC_temp.rds")  # save the ggplot objec


##--------------------------------------------------##
## References
# Simpson, 2014: https://fromthebottomoftheheap.net/2014/05/09/modelling-seasonal-data-with-gam/
# Hayes et al, 2020: https://doi.org/10.1002/lol2.10164
# Simpson 2020 : https://stats.stackexchange.com/questions/471267/plotting-gams-on-response-scale-with-multiple-smooth-and-linear-terms
# Wilk et al., 2018: https://doi.org/10.1029/2018JG004506 (plotting response values)
