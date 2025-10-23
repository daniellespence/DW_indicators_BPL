 ###############################################################################
####### BPWTP Data -- GAMs ##########
####### Danielle Spence ##########
####### Created 5/4/2023 ########
###############################################################################
### Clear memory
rm(list = ls())

library(pacman)
p_load(tidyverse, ggplot2, mgcv, gratia, readr, GGally, dplyr, lubridate,
        cowplot, tibble,viridis, install = TRUE)

setwd("C:/Users/danis/OneDrive/R/DW indicators")

papertheme <- theme_bw(base_size=12, base_family = 'Arial') 

##----------------------------------------------------------------------------##
## 1. Read in data 
##----------------------------------------------------------------------------##

df <- read_csv("data/bpgamdataCLEAN_Turb.csv")

df$Date<- as.Date(df$Date, "%m/%d/%Y")

# add in Year and nMonth for numeric month and a proper Date class

df <- mutate(df,
             year = as.numeric(format(Date,'%Y')),
             DOY = as.numeric(format(Date,'%j')),
             nMonth = as.numeric(format(Date,'%m')),
             week = as.numeric(format(Date, '%W')))%>% 
  filter(year %in% c(1990:2022))
# 974 obs


## order the data
df <- df[with(df, order(Date)),]


# rename for ease of use
df <- df%>%
  dplyr::rename(
    date = "Date",
   turb = "Turbidity",
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


#Remove rows without turb data
df <- df[complete.cases(df$turb),]
# 974 obs

df$DIN <- rowSums(df[,c("NO3", "NH3")], na.rm=TRUE)


## add in PDO and SOI 

clim <- read_csv("data/SOI_PDO data no lag.csv")

df1<- merge(df, clim, by="date")

##----------------------------------------------------------------------------##
## 3. Estimate annual and seasonal variability
##----------------------------------------------------------------------------##

df_clean <- df1 %>%
  ungroup() %>%                    # Remove any previous grouping
  select(year, turb)


library(dplyr)

df_clean %>%
  group_by(year) %>%
  summarise(
    mean_year = mean(turb, na.rm = TRUE),
    seasonal_range = max(turb, na.rm = TRUE) - min(turb, na.rm = TRUE)
  ) %>%
  summarise(
    avg_seasonal_range = mean(seasonal_range, na.rm = TRUE),
    avg_seasonal_range_pct = mean(seasonal_range / mean_year * 100, na.rm = TRUE),
    between_year_range = max(mean_year) - min(mean_year),
    between_year_range_pct = (max(mean_year) - min(mean_year)) / mean(mean_year) * 100
  )


##----------------------------------------------------------------------------##
## 5. Run the full model (see Fitting script)
##----------------------------------------------------------------------------##


m2 <- gam(turb ~ s(SRP) +  
             s(DIN)+ 
             s(QLD)+
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
            knots=list(DOY=c(0, 366.5)),
           select = TRUE,
           data = df1, method = "REML", family = Gamma(link = "log"))

 saveRDS(m2, file = "modelturb.rds")
 sink("model_summaryturb.txt")
 summary(m2)
 sink()
 
summary(m2)
# r sq = 0.514

k.check(m2)# k-index looks good 

p1<- draw(m2,  residuals = TRUE)& theme_bw() 
p1

#ggsave('output/turb GAM 1995 on TN TP.png', p1, height = 10, width  = 10)

##----------------------------------------------------------------------------##
## 6. Assess model fit & autocorrelation
##----------------------------------------------------------------------------##

## i. Appraise model fit
p2<-appraise(m2, point_col = 'steelblue', point_alpha = 0.5, n_bins = 'fd') & 
  theme(plot.tag = element_text(face = 'bold')) & theme_bw()
p2

saveRDS(p2, "turb_residuals.rds")


ggsave('output/turb GAM residuals no temp.png', p2, height = 10, width  = 10)


# ii. Check for autocorrelation
layout(matrix(1:2, ncol = 2))
acf(resid(m2), lag.max = 36, main = "ACF") 
pacf(resid(m2), lag.max = 36, main = "pACF")
layout(1)


#----------------------------------------------------------------------------##
## 7. Plotting on the response scale/fitted values
##----------------------------------------------------------------------------##

## ---------- DIN -------------------##


new_data_din <- with(df1, expand.grid(DIN = seq(min(DIN), max(DIN),length = 200),
                                           SRP = median(SRP, na.rm=TRUE),
                                           QLD = median(QLD),
                                           SOI_3MON_AVG = median(SOI_3MON_AVG, na.rm=TRUE),
                                           PDO_3MON_AVG = median(PDO_3MON_AVG, na.rm=TRUE),
                                           year = median(year),
                                           DOY = median(DOY)))

din.pred <- predict(m2, newdata = new_data_din, type = "terms", se.fit = TRUE)

whichCols <- grep("DIN", colnames(din.pred$fit))
whichColsSE <- grep("DIN", colnames(din.pred$se.fit))
new_data_din <- cbind(new_data_din, Fitted = din.pred$fit[, whichCols], 
                           se.Fitted = din.pred$se.fit[, whichColsSE])
limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

## make into original limits
new_data_din <- with(new_data_din, transform(new_data_din, Fittedplus = Fitted + se.Fitted))
new_data_din <- with(new_data_din, transform(new_data_din, Fittedminus = Fitted - se.Fitted))

shiftdin <- attr(predict(m2, newdata = new_data_din, type = "iterms"), "constant")
din.pdatnorm <- new_data_din
din.pdatnorm <- with(din.pdatnorm, transform(din.pdatnorm, Fitted = Fitted + shiftdin, 
                                                       Fittedplus = Fittedplus + shiftdin, 
                                                       Fittedminus = Fittedminus + shiftdin))

din.pdatnorm$Fitted<-exp(din.pdatnorm$Fitted)
din.pdatnorm$Fittedplus<-exp(din.pdatnorm$Fittedplus)
din.pdatnorm$Fittedminus<-exp(din.pdatnorm$Fittedminus)

dinquants <- quantile(df1$DIN, c(.05,.95), na.rm = TRUE)

DINplot <- ggplot(din.pdatnorm, aes(x = DIN, y = Fitted)) +
  papertheme +
  scale_y_continuous(limits=c(0, 10))+
  annotate("rect", xmin=dinquants[1], xmax=dinquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
  annotate("text", 
           x = -Inf, y = Inf,                  # Top-left corner
           hjust = -0.1, vjust = 1.5,          # Adjust alignment to keep text inside plot
           label = "italic(p) == 0.0579", parse = TRUE)+
  geom_line() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25, fill = '#165459B2') +  
  geom_rug(aes(x=DIN), data = df1, stat = "identity", position = "identity", 
           sides = "b", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
  xlab(expression(paste("DIN"~"(mg/L)"))) + ylab("Turbidity (NTU)")

DINplot


saveRDS(DINplot, "DIN_TURB_notemp.rds")



### ---------- SRP -------------------##

new_data_srp <- with(df1, expand.grid(SRP = seq(min(SRP, na.rm = TRUE), max(SRP, na.rm = TRUE), length = 200),
                                     DIN = median(DIN),
                                     QLD = median(QLD),
                                     SOI_3MON_AVG = median(W_temp, na.rm=TRUE), 
                                     PDO_3MON_AVG = median(PDO_3MON_AVG, na.rm=TRUE),
                                     year= median(year),
                                     DOY = median(DOY)))


SRP.pred <- predict(m2, newdata = new_data_srp, type = "terms", se.fit = TRUE)




whichCols <- grep("SRP", colnames(SRP.pred$fit))
whichColsSE <- grep("SRP", colnames(SRP.pred$se.fit))
new_data_srp <- cbind(new_data_srp, Fitted = SRP.pred$fit[, whichCols], 
                      se.Fitted = SRP.pred$se.fit[, whichColsSE])


limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

## make into original limits
new_data_srp <- with(new_data_srp, transform(new_data_srp, Fittedplus = Fitted + se.Fitted))
new_data_srp <- with(new_data_srp, transform(new_data_srp, Fittedminus = Fitted - se.Fitted))

shiftSRP <- attr(predict(m2, newdata = new_data_srp, type = "iterms"), "constant")
SRP.pdatnorm <- new_data_srp
SRP.pdatnorm <- with(SRP.pdatnorm, transform(SRP.pdatnorm, Fitted = Fitted + shiftSRP, 
                                             Fittedplus = Fittedplus + shiftSRP, 
                                             Fittedminus = Fittedminus + shiftSRP))

# backtransform fitted data..
SRP.pdatnorm$Fitted<-exp(SRP.pdatnorm$Fitted)
SRP.pdatnorm$Fittedplus<-exp(SRP.pdatnorm$Fittedplus)
SRP.pdatnorm$Fittedminus<-exp(SRP.pdatnorm$Fittedminus)


SRPquants <- quantile(df1$SRP, c(.05,.95), na.rm = TRUE)


SRPplot <- ggplot(SRP.pdatnorm, aes(x = SRP, y = Fitted)) +
  papertheme+
  scale_y_continuous(limits=c(0, 10))+
  annotate("rect", xmin=SRPquants[1], xmax=SRPquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
  annotate("text", 
           x = -Inf, y = Inf,                  # Top-left corner
           hjust = -0.1, vjust = 1.5,          # Adjust alignment to keep text inside plot
           label = "italic(p) == 0.0790", parse = TRUE)+
  geom_line() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25, fill = '#165459B2') +  
   geom_rug(aes(x=SRP), data = df1, stat = "identity", position = "identity", 
           sides = "b", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
  xlab(expression(paste("SRP ","(", mu, "g/L)"))) + ylab("Turbidity (NTU)")

SRPplot


saveRDS(SRPplot, "SRP_TURB_notemp.rds")


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
  scale_y_continuous(limits=c(0, 10))+
  annotate("rect", xmin=QLDquants[1], xmax=QLDquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
  annotate("text", 
           x = -Inf, y = Inf,                  # Top-left corner
           hjust = -0.1, vjust = 1.5,          # Adjust alignment to keep text inside plot
           label = expression(italic("p")*" < 0.0001***"))+
  geom_line() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25, fill = '#165459B2') +  
 geom_rug(aes(x=QLD), data = df1, stat = "identity", position = "identity", 
           sides = "b", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
  xlab(expression(paste("QLD (",m^3, "/", s,")"))) + ylab("Turbidity (NTU)")

QLDplot


saveRDS(QLDplot, "QLD_TURB_notemp.rds")

# range of effect

max_idx <- which.max(QLD.pdatnorm$Fitted)
min_idx <- which.min(QLD.pdatnorm$Fitted)
range_effect_QLD <- QLD.pdatnorm$Fitted[max_idx] - QLD.pdatnorm$Fitted[min_idx]
# 2.62

#lower range of CI
range_lwr <- QLD.pdatnorm$Fittedminus[max_idx] - QLD.pdatnorm$Fittedplus[min_idx]
range_upr <- QLD.pdatnorm$Fittedplus[max_idx] - QLD.pdatnorm$Fittedminus[min_idx]


# compare to seasonal and interannual variability

# Intra-annual variation (average seasonal swing)
within_year_var <- df1 %>%
  group_by(year) %>%
  summarise(range_within = max(turb) - min(turb), .groups = "drop") %>%
  summarise(avg_within = mean(range_within))

# Interannual variation (range of annual means)
between_year_var <- df1 %>%
  group_by(year) %>%
  summarise(mean_year = mean(turb), .groups = "drop") %>%
  summarise(range_between = max(mean_year) - min(mean_year))


# relative_to_within = 
range_effect_QLD / within_year_var$avg_within

#0.178

#relative_to_between = 
range_effect_QLD / between_year_var$range_between
#0.294



p_all<- plot_grid(DINplot,SRPplot, tempplot,QLDplot, ncol = 2)
p_all


#ggsave('output/turb GAM_Response scale .png', p_all, height = 10, width  = 10)


## ---------- SOI * PDO -------------------##



new_data_SOI_3MON_AVG <- with(df1, expand.grid(SOI_3MON_AVG = seq(min(SOI_3MON_AVG), max(SOI_3MON_AVG), length = 200),
                                               PDO_3MON_AVG = seq(min(PDO_3MON_AVG), max(PDO_3MON_AVG), length = 200),
                                               SRP = median(SRP, na.rm=TRUE),
                                               DIN = median(DIN, na.rm=TRUE),
                                               QLD = median(QLD, na.rm=TRUE),
                                               year= median(year),
                                               DOY = median(DOY)))


SOI_3MON_AVG.pred <- predict(m2, newdata = new_data_SOI_3MON_AVG, type = "terms", se.fit = TRUE)

whichCols <- grep("PDO_3MON_AVG,SOI_3MON_AVG", colnames(SOI_3MON_AVG.pred$fit))
whichColsSE <- grep("PDO_3MON_AVG,SOI_3MON_AVG", colnames(SOI_3MON_AVG.pred$se.fit))
new_data_SOI_3MON_AVG <- cbind(new_data_SOI_3MON_AVG, Fitted = SOI_3MON_AVG.pred$fit[, whichCols], 
                               se.Fitted = SOI_3MON_AVG.pred$se.fit[, whichColsSE])

limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

## make into original limits
new_data_SOI_3MON_AVG <- with(new_data_SOI_3MON_AVG, transform(new_data_SOI_3MON_AVG, Fittedplus = Fitted + se.Fitted))
new_data_SOI_3MON_AVG <- with(new_data_SOI_3MON_AVG, transform(new_data_SOI_3MON_AVG, Fittedminus = Fitted - se.Fitted))

shiftcomb <- attr(SOI_3MON_AVG.pred, "constant")
SOI_3MON_AVG.pdatnorm <- new_data_SOI_3MON_AVG
SOI_3MON_AVG.pdatnorm <- with(SOI_3MON_AVG.pdatnorm, transform(SOI_3MON_AVG.pdatnorm, Fitted = Fitted + shiftcomb, 
                                                               Fittedplus = Fittedplus + shiftcomb, 
                                                               Fittedminus = Fittedminus + shiftcomb))

toofar <- exclude.too.far(SOI_3MON_AVG.pdatnorm$PDO_3MON_AVG, SOI_3MON_AVG.pdatnorm$SOI_3MON_AVG, df1$PDO_3MON_AVG, df1$SOI_3MON_AVG, dist=0.1)
SOI_3MON_AVG.pdatnorm$turb <- SOI_3MON_AVG.pdatnorm$Fitted
SOI_3MON_AVG.pdatnorm$turb[toofar] <- NA

SOI_3MON_AVG.pdatnorm$turb<-exp(SOI_3MON_AVG.pdatnorm$turb)
SOI_3MON_AVG.pdatnorm$Fittedplus<-exp(SOI_3MON_AVG.pdatnorm$Fittedplus)
SOI_3MON_AVG.pdatnorm$Fittedminus<-exp(SOI_3MON_AVG.pdatnorm$Fittedminus)

names(new_data_SOI_3MON_AVG)[which(names(new_data_SOI_3MON_AVG)=='SOI_3MON_AVG.pred')] <- 'turb'

comboplot <- ggplot(SOI_3MON_AVG.pdatnorm, aes(x = PDO_3MON_AVG, y = SOI_3MON_AVG, z=turb)) + 
  theme_bw(base_family = 'Arial') +
  theme(legend.position='top') +
  geom_raster(aes(fill=turb)) + 
  scale_fill_distiller(palette = "Spectral", direction = -1, na.value='transparent') +
  geom_point(data=df1, aes(x=PDO_3MON_AVG, y=SOI_3MON_AVG, z=NULL)) +
  geom_contour(colour = "black", binwidth = 0.5) +
  scale_x_continuous(expand = c(0, 0)) +  # Remove padding on x-axis
  scale_y_continuous(expand = c(0, 0)) +  # Remove padding on y-axis
  theme(legend.key.width=unit(1,"cm"))+
  xlab("PDO") + ylab("SOI") +
  labs(fill="Turbidity (NTU)")+
  theme(
    axis.text.x = element_text(color = ifelse(seq(-2.5, 2.5, by = 1) < -0.5, "blue",
                                              ifelse(seq(-2.5, 2.5, by = 1) > 0, "red", "black"))),
    axis.text.y = element_text(color = ifelse(seq(-2.5, 2.5, by = 1) < -0.5, "red",
                                              ifelse(seq(-2.5, 2.5, by = 1) > 0, "blue", "black")))
  )+
  theme(legend.position = "top")

comboplot


#ggsave('output/SOI PDO INTERACTION Response scale turb.png', comboplot, height = 8, width  = 10)


saveRDS(comboplot, "SOIPDO_TURB_notemp.rds")


# range of control

SOI_3MON_AVG.pdatnorm1<- na.omit(SOI_3MON_AVG.pdatnorm)

max_idx <- which.max(SOI_3MON_AVG.pdatnorm1$turb)
min_idx <- which.min(SOI_3MON_AVG.pdatnorm1$turb)
range_effect <- SOI_3MON_AVG.pdatnorm1$turb[max_idx] - SOI_3MON_AVG.pdatnorm1$turb[min_idx]


#lower range of CI
range_lwr <- SOI_3MON_AVG.pdatnorm1$Fittedminus[max_idx] - SOI_3MON_AVG.pdatnorm1$Fittedplus[min_idx]
range_upr <- SOI_3MON_AVG.pdatnorm1$Fittedplus[max_idx] - SOI_3MON_AVG.pdatnorm1$Fittedminus[min_idx]


# relative_to_within = 
range_effect_SOI / within_year_var$avg_within

#1.06

#relative_to_between = 
range_effect_SOI / between_year_var$range_between
#0.55

##---------------------------- temporal change------------------------------------- ####


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
year.pdatnorm$TURB <- year.pdatnorm$Fitted
year.pdatnorm$TURB[toofar] <- NA


year.pdatnorm$TURB<-exp(year.pdatnorm$TURB)


names(new_data_year)[which(names(new_data_year)=='year.pred')] <- 'TURB'

comboplot_time <- ggplot(year.pdatnorm, aes(x = year, y = DOY, z=TURB)) + 
  theme_minimal(base_size=14, base_family = 'Arial') +
  theme(legend.position='top') +
  geom_raster(aes(fill=TURB)) +
  scale_fill_distiller(palette = "Spectral", direction = -1, na.value='transparent') +
  geom_point(data=df1, aes(x=year, y=DOY, z=NULL)) +
  geom_contour(colour = "black", binwidth = 5) +
  scale_y_continuous(expand = c(0, 0), breaks = c(1, 91,  182, 274, 335), labels = c("Jan", "Apr", "Jul", "Oct", "Dec")) +
  scale_x_continuous(expand = c(0, 0)) +  # Remove padding on x-axis
  theme(legend.key.width=unit(1.5,"cm")) +
  xlab("Year") + ylab("DOY") +
  labs(fill=expression(paste("Turbidity (NTU)")))+
  theme(legend.position = "top")


comboplot_time

#ggsave('output/TURB change over time.png', comboplot_time, height = 6, width  = 8)


saveRDS(comboplot_time, "time_TURB.rds")  # save the ggplot object



##--------------------------------------------------##
## References
# Simpson, 2014: https://fromthebottomoftheheap.net/2014/05/09/modelling-seasonal-data-with-gam/
# Hayes et al, 2020: https://doi.org/10.1002/lol2.10164
# Simpson 2020 : https://stats.stackexchange.com/questions/471267/plotting-gams-on-response-scale-with-multiple-smooth-and-linear-terms
# Wilk et al., 2018: https://doi.org/10.1029/2018JG004506 (plotting response values)
