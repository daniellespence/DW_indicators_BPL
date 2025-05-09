 ###############################################################################
####### BPWTP Data -- GAMs ##########
####### Danielle Spence ##########
####### Created 5/4/2023 ########
###############################################################################
### Clear memory
rm(list = ls())

library(pacman)
p_load(tidyverse, ggplot2, mgcv, gratia, readr, GGally, dplyr,  mgcViz, lubridate,
        cowplot, tibble, install = TRUE)

setwd("C:/Users/danis/OneDrive/R/DW indicators")

papertheme <- theme_bw(base_size=12, base_family = 'Arial') 

##----------------------------------------------------------------------------##
## 1. Read in data 
##----------------------------------------------------------------------------##

df <- read_csv("data/bpgamdataCLEAN_TON.csv")

df$Date<- as.Date(df$Date, "%m/%d/%Y")

# add in Year and nMonth for numeric month and a proper Date class

df <- mutate(df,
             year = as.numeric(format(Date,'%Y')),
             DOY = as.numeric(format(Date,'%j')),
             nMonth = as.numeric(format(Date,'%m')),
             week = as.numeric(format(Date, '%W')))%>% 
  filter(year %in% c(1990:2022))
# 1435 obs

## order the data
df <- df[with(df, order(Date)),]


# rename for ease of use
df <- df%>%
  dplyr::rename(
    date = "Date",
   Odour = "Odour_TON",
   NO3 = "NO3_mg.L",
    NH3 = "NH3_mg.L",
    SRP = "SRP_ug.L",
    TP = "TP_ug.L",
    W_temp = "Temp_C",
   QBP = "combined_05JG004.cms",
    QLD = "SK05JG006.cms",
    QWS = "RC_IC_cms",
    orgN = "Org_N_mg.L")


#Remove rows without Odour data
df <- df[complete.cases(df$Odour),]
# 1398 obs

df$DIN <- rowSums(df[,c("NO3", "NH3")], na.rm=TRUE)

## add in PDO and SOI 

clim <- read_csv("data/SOI_PDO data no lag.csv")

df1<- merge(df, clim, by="date")

##----------------------------------------------------------------------------##
## 3. Check annual and seasonal variability
##----------------------------------------------------------------------------##


df_clean <- df1 %>%
  ungroup() %>%                    # Remove any previous grouping
  select(year, Odour)    


dx<- df_clean %>%
  group_by(year) %>%
  summarize(
    mean_conc = mean(Odour, na.rm = TRUE),
    range_percent = ((max(Odour, na.rm = TRUE) - min(Odour, na.rm = TRUE)) / mean_conc) * 100
  )
dx

mean(dx$range_percent)

# Step 1: Calculate the mean concentration for each year
yearly_means <- df1 %>%
  group_by(year) %>%
  summarize(mean_conc = mean(Odour, na.rm = TRUE), .groups = "drop")

# Step 2: Compute the max, min, and average of yearly means
summary_stats <- yearly_means %>%
  summarize(
    max_mean = max(mean_conc),
    min_mean = min(mean_conc),
    overall_mean = mean(mean_conc),
    percent_variation = ((max_mean - min_mean) / overall_mean) * 100
  )

print(summary_stats)


##----------------------------------------------------------------------------##
## 5. Run the full model (see Fitting script)
##----------------------------------------------------------------------------##

m2 <- gam(Odour ~ s(SRP) +  
             s(DIN)+ 
             #s(W_temp)+
             s(QLD)+
            #ti(DOY, QLD, bs = c("cc", "tp"))+ 
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs=c('cr', 'cc')),
             #s(DOY, bs = 'cc')+
             #s(Year, bs = 're'),
            knots=list(DOY=c(0, 366.5)),
           select = TRUE,
           data = df1, method = "REML", family = Gamma(link = "log"))



saveRDS(m2, file = "modelOdour.rds")

 # sink("model_summaryOdour_notemp.txt")
 # summary(m2)
 # sink()


summary(m2) # Deviance explained =75%, REML = 4010, r2 = 0.493, n=861
# all sig

k.check(m2)# k-index looks good 

p1<- draw(m2,  residuals = TRUE)& theme_bw() 
p1

#ggsave('output/Odour GAM 1990.png', p1, height = 10, width  = 10)

##----------------------------------------------------------------------------##
## 6. Assess model fit & autocorrelation
##----------------------------------------------------------------------------##

## i. Appraise model fit
p2<-appraise(m2, point_col = 'steelblue', point_alpha = 0.5, n_bins = 'fd') & 
  theme(plot.tag = element_text(face = 'bold')) & theme_bw()
p2

ggsave('output/Odour GAM residuals.png', p2, height = 10, width  = 10)


# ii. Check for autocorrelation
layout(matrix(1:2, ncol = 2))
acf(resid(m), lag.max = 36, main = "ACF") 
pacf(resid(m), lag.max = 36, main = "pACF")
layout(1)


#----------------------------------------------------------------------------##
## 7. Plotting on the response scale/fitted values
##----------------------------------------------------------------------------##

# how much does TON vary on an annual basis?

# Group by year and calculate variation metrics
annual_var <- aggregate(Odour ~ year, data = df1, function(x) {
  if (all(is.na(x))) return(NA)
  max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
})

mean(annual_var$Odour)
#230

min(df1$Odour)
max(df1$Odour)
#The average range in DOC concentration each year is 230 TON
# the range is 2 to 800 TON

## ---------- DIN -------------------##


new_data_DIN <- with(df1, expand.grid(DIN = seq(min(DIN), max(DIN),length = 200),
                                           SRP = median(SRP, na.rm=TRUE),
                                           #W_temp = median(W_temp, na.rm=TRUE),
                                           QLD = median(QLD),
                                           SOI_3MON_AVG = median(SOI_3MON_AVG, na.rm=TRUE),
                                           PDO_3MON_AVG = median(PDO_3MON_AVG, na.rm=TRUE),
                                           year = median(year),
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
  scale_y_continuous(limits=c(15, 150))+
  annotate("rect", xmin=DINquants[1], xmax=DINquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
  annotate("text", 
           x = -Inf, y = Inf,                  # Top-left corner
           hjust = -0.1, vjust = 1.5,          # Adjust alignment to keep text inside plot
           label = expression(italic("p")*" = 0.0344*"))+
  geom_line() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25, fill = '#165459B2') +  
  geom_rug(aes(x=DIN), data = df1, stat = "identity", position = "identity", 
           sides = "b", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
  xlab(expression(paste("DIN"~"("*"mg"~L^{-1}*")"))) + ylab("Odour (T.O.N)")

#DINplot


saveRDS(DINplot, "DIN_Odour_notemp.rds")  # save the ggplot objec

# range of control
max(DIN.pdatnorm$Fitted)-min(DIN.pdatnorm$Fitted)
# = 15.3

### ---------- SRP -------------------##

new_data_SRP <- with(df1, expand.grid(SRP = seq(min(SRP, na.rm = TRUE), max(SRP, na.rm = TRUE), length = 200),
                                     DIN = median(DIN),
                                     #W_temp = median(W_temp, na.rm=TRUE),
                                     QLD = median(QLD),
                                     SOI_3MON_AVG = median(W_temp, na.rm=TRUE), 
                                     PDO_3MON_AVG = median(PDO_3MON_AVG, na.rm=TRUE),
                                     year = median(year),
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

# backtransform fitted data...
SRP.pdatnorm$Fitted<-exp(SRP.pdatnorm$Fitted)
SRP.pdatnorm$Fittedplus<-exp(SRP.pdatnorm$Fittedplus)
SRP.pdatnorm$Fittedminus<-exp(SRP.pdatnorm$Fittedminus)


SRPquants <- quantile(df1$SRP, c(.05,.95), na.rm = TRUE)


SRPplot <- ggplot(SRP.pdatnorm, aes(x = SRP, y = Fitted)) +
  papertheme+
  scale_y_continuous(limits=c(15, 150))+
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
  xlab(expression(paste("SRP ","(", mu, "g"~L^{-1}*")"))) + ylab("Odour (T.O.N)")

#SRPplot



saveRDS(SRPplot, "SRP_Odour_notemp.rds")  # save the ggplot objec


# range of control
max(SRP.pdatnorm$Fitted)-min(SRP.pdatnorm$Fitted)
# = 73.3


## ---------- QLD -------------------##

new_data_QLD <- with(df1, expand.grid(QLD = seq(min(QLD), max(QLD), length = 200),
                                      DIN = median(DIN),
                                      SRP = median(SRP, na.rm=TRUE),
                                      #W_temp = median(W_temp, na.rm=TRUE),
                                      SOI_3MON_AVG = median(W_temp, na.rm=TRUE), 
                                      PDO_3MON_AVG = median(PDO_3MON_AVG, na.rm=TRUE),
                                      year = median(year),
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
  scale_y_continuous(limits=c(15, 150))+
  annotate("rect", xmin=QLDquants[1], xmax=QLDquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
  annotate("text", 
           x = -Inf, y = Inf,                  # Top-left corner
           hjust = -0.1, vjust = 1.5,          # Adjust alignment to keep text inside plot
           label = expression(italic("p")*" = 0.431"))+
  geom_line() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25, fill = '#165459B2') +  
 geom_rug(aes(x=QLD), data = df1, stat = "identity", position = "identity", 
           sides = "b", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
  xlab(expression(paste("QLD (",m^3, "/", s,")"))) + ylab("Odour (T.O.N)")

#QLDplot


saveRDS(QLDplot, "QLD_Odour_notemp.rds")  # save the ggplot objec

# range of control
max(QLD.pdatnorm$Fitted)-min(QLD.pdatnorm$Fitted)
# = 12.5

# p_all<- plot_grid(SRPplot,DINplot, tempplot,QLDplot, ncol = 2)
# p_all


#ggsave('output/Odour GAM_Response scale.png', p_all, height = 10, width  = 10)

## ---------- SOI * PDO -------------------##


new_data_SOI_3MON_AVG <- with(df1, expand.grid(SOI_3MON_AVG = seq(min(SOI_3MON_AVG), max(SOI_3MON_AVG), length = 200),
                                               PDO_3MON_AVG = seq(min(PDO_3MON_AVG), max(PDO_3MON_AVG), length = 200),
                                               SRP = median(SRP, na.rm =TRUE),
                                               DIN = median(DIN, na.rm=TRUE),
                                               #W_temp = median(W_temp, na.rm=TRUE),
                                               QLD = median(QLD, na.rm=TRUE),
                                               year = median(year),
                                               DOY = median(DOY)))

SOI_3MON_AVG.pred <- predict(m2, newdata = new_data_SOI_3MON_AVG, type = "terms")

whichCols <- grep("PDO_3MON_AVG,SOI_3MON_AVG", colnames(SOI_3MON_AVG.pred))

new_data_SOI_3MON_AVG <- cbind(new_data_SOI_3MON_AVG, Fitted = SOI_3MON_AVG.pred[, whichCols])

shiftcomb <- attr(SOI_3MON_AVG.pred, "constant")
SOI_3MON_AVG.pdatnorm <- new_data_SOI_3MON_AVG
SOI_3MON_AVG.pdatnorm <- with(SOI_3MON_AVG.pdatnorm, transform(SOI_3MON_AVG.pdatnorm, Fitted = Fitted + shiftcomb))

toofar <- exclude.too.far(SOI_3MON_AVG.pdatnorm$PDO_3MON_AVG, SOI_3MON_AVG.pdatnorm$SOI_3MON_AVG, df1$PDO_3MON_AVG, df1$SOI_3MON_AVG, dist=0.1)
SOI_3MON_AVG.pdatnorm$Odour <- SOI_3MON_AVG.pdatnorm$Fitted
SOI_3MON_AVG.pdatnorm$Odour[toofar] <- NA


SOI_3MON_AVG.pdatnorm$Odour<-exp(SOI_3MON_AVG.pdatnorm$Odour)

names(new_data_SOI_3MON_AVG)[which(names(new_data_SOI_3MON_AVG)=='SOI_3MON_AVG.pred')] <- 'Odour'

comboplot <- ggplot(SOI_3MON_AVG.pdatnorm, aes(x = PDO_3MON_AVG, y = SOI_3MON_AVG, z=Odour)) + 
  theme_bw(base_size=14) +
  theme(legend.position='top') +
  geom_raster(aes(fill=Odour)) + 
  scale_fill_distiller(palette = "Spectral", direction = -1, na.value='transparent') +
  geom_point(data=df1, aes(x=PDO_3MON_AVG, y=SOI_3MON_AVG, z=NULL)) +
  scale_x_continuous(expand = c(0, 0)) +  # Remove padding on x-axis
  scale_y_continuous(expand = c(0, 0)) +  # Remove padding on y-axis
  geom_contour(colour = "black", binwidth = 10) +
  theme(legend.key.width=unit(2,"cm"))+
  xlab("PDO") + ylab("SOI") +
  labs(fill="Odour (T.O.N)")+
  theme(
    axis.text.x = element_text(color = ifelse(seq(-2.5, 2.5, by = 1) < -0.5, "blue",
                                              ifelse(seq(-2.5, 2.5, by = 1) > 0, "red", "black"))),
    axis.text.y = element_text(color = ifelse(seq(-2.5, 2.5, by = 1) < -0.5, "red",
                                              ifelse(seq(-2.5, 2.5, by = 1) > 0, "blue", "black")))
  )+
  theme(legend.position = "top")

comboplot


#ggsave('output/SOI PDO INTERACTION Response scale Odour.png', comboplot, height = 8, width  = 10)


saveRDS(comboplot, "SOIPDO_Odour.rds")  # save the ggplot objec

# range of control
max(SOI_3MON_AVG.pdatnorm$Odour, na.rm = TRUE)- min(SOI_3MON_AVG.pdatnorm$Odour, na.rm = TRUE)
# = 119.5
# = 14.8 min and 134 max


### -------------------------temporal change ------------------------------ ####


new_data_year <- with(df1, expand.grid(DOY = seq(min(DOY), max(DOY),length = 200),
                                       year = seq(min(year), max(year),length = 200),
                                       SOI_3MON_AVG = median(SOI_3MON_AVG),
                                       PDO_3MON_AVG = median(PDO_3MON_AVG),
                                       SRP = median(SRP, na.rm=TRUE),
                                       DIN = median(DIN, na.rm=TRUE),
                                       #W_temp = median(W_temp, na.rm=TRUE),
                                       QLD = median(QLD, na.rm=TRUE)))

year.pred <- predict(m2, newdata = new_data_year, type = "terms")

whichCols <- grep("year,DOY", colnames(year.pred))

new_data_year <- cbind(new_data_year, Fitted = year.pred[, whichCols])

shiftcomb <- attr(year.pred, "constant")
year.pdatnorm <- new_data_year
year.pdatnorm <- with(year.pdatnorm, transform(year.pdatnorm, Fitted = Fitted + shiftcomb))

toofar <- exclude.too.far(year.pdatnorm$DOY, year.pdatnorm$year, df1$DOY, df1$year, dist=0.1)
year.pdatnorm$Odour <- year.pdatnorm$Fitted
year.pdatnorm$Odour[toofar] <- NA


year.pdatnorm$Odour<-exp(year.pdatnorm$Odour)


names(new_data_year)[which(names(new_data_year)=='year.pred')] <- 'Odour'



comboplot_time <- ggplot(year.pdatnorm, aes(x = year, y = DOY, z=Odour)) + 
  theme_minimal(base_size=14, base_family = 'Arial') +
  theme(legend.position='top') +
  geom_raster(aes(fill=Odour)) + 
  scale_fill_distiller(palette = "Spectral", direction = -1, na.value='transparent') +
  geom_point(data=df1, aes(x=year, y=DOY, z=NULL)) +
  geom_contour(colour = "black", binwidth = 25) +
  theme(legend.key.width=unit(1.5,"cm")) +
   xlab("Year") + ylab("DOY") +
  scale_y_continuous(expand = c(0, 0), breaks = c(1, 91,  182, 274, 335), labels = c("Jan", "Apr", "Jul", "Oct", "Dec")) +
  scale_x_continuous(expand = c(0, 0)) +  # Remove padding on x-axis
  labs(fill=expression(paste("Odour (T.O.N)")))+
  theme(legend.position = "top")


#ggtitle(expression(paste("Chlorophyll ", italic("a"), " (", mu, "g L"^-1*")")))


comboplot_time


#ggsave('output/Odour change over time.png', comboplot_time, height = 6, width  = 8)


saveRDS(comboplot_time, "time_Odour_notemp.rds")  # save the ggplot objecT



### ------------SOI PDO ----------------#


N <- 200
simSOI <- c(-1, 0.3, 1.1)
SOIgroup <- factor(rep(c('-1', '0.3', '1.1'), each=200))  # 7*300
reptimes <- length(simSOI) #3
preddf <- data.frame(PDO_3MON_AVG = rep(seq(min(df1$PDO_3MON_AVG, na.rm=TRUE),
                                            max(df1$PDO_3MON_AVG, na.rm=TRUE),
                                            length = N), times=reptimes),
                     SOI_3MON_AVG = rep(simSOI, each = N)) #300
superN <- nrow(preddf)


SOI.pdat <- with(df1,
                 data.frame(PDO_3MON_AVG = rep(preddf$PDO_3MON_AVG), # 300*7
                            SOI_3MON_AVG = rep(preddf$SOI_3MON_AVG), # 300*7
                            SRP = median(SRP, na.rm=TRUE),
                            DIN = median(DIN, na.rm=TRUE),
                            #W_temp = median(W_temp, na.rm=TRUE),
                            QLD = median(QLD, na.rm=TRUE),
                            year = rep(2003),
                            DOY = median(DOY)))

SOI.pred <- predict(m2, newdata = SOI.pdat, type = "link")
SOI.pdat <- cbind(SOI.pdat, SOI.pred)
SOI.pdat$SOIgroup <- SOIgroup

names(SOI.pdat)[which(names(SOI.pdat)=='SOI.pred')] <- 'Odour'

# backtransform from tweedie distribution
SOI.pdat$Odour<- exp(SOI.pdat$Odour)

# need to take away places where PDO*SOI combo has not occurred in the data!!
toofar <- exclude.too.far(SOI.pdat$PDO_3MON_AVG, SOI.pdat$SOI_3MON_AVG, df1$PDO_3MON_AVG, df1$SOI_3MON_AVG, dist=0.1)
SOI.pdat$Odour[toofar] <- NA



#reorder factors
SOI.pdat$SOIgroup <- factor(SOI.pdat$SOIgroup, levels = c('-1', '0.3', '1.1'))

SOIplot <- ggplot(SOI.pdat, aes(x = PDO_3MON_AVG, y = Odour, group= SOI_3MON_AVG, col=SOIgroup, 
                                lty=SOIgroup)) +
  theme_bw(base_size=14) +
  theme(legend.position='top') +
  geom_line() +
  scale_y_continuous(limits=c(50, 300))+
  scale_color_manual(name=expression(paste(bold('')~'      SOI')), values = c("red", "black", "blue"))+
  scale_linetype_manual(name=expression(paste(bold('')~'      SOI')), values = c("solid", "solid","longdash")) +
  xlab('PDO') + ylab('Odour (T.O.N)')+
  theme(
    axis.text.x = element_text(color = ifelse(seq(-2.5, 2.5, by = 1) < 0, "blue",
                                              ifelse(seq(-2.5, 2.5, by = 1) > 0.5, "red", "black")))
  )

## '#e66101','#fdb863','#b2abd2','#5e3c99' (order red:light:dark)

SOIplot


saveRDS(SOIplot, "soi_Odour.rds")  # save the ggplot object


## now slices for PDO
N <- 200
simPDO <- c(-1.5, 0, 1.6)
PDOgroup <- factor(rep(c('-1.5', '0', '1.6'), each=200)) 

reptimes <- length(simPDO) #3
preddf <- data.frame(SOI_3MON_AVG = rep(seq(min(df1$SOI_3MON_AVG, na.rm=TRUE),
                                            max(df1$SOI_3MON_AVG, na.rm=TRUE),
                                            length = N), times=reptimes),
                     PDO_3MON_AVG = rep(simPDO, each = N)) #300
superN <- nrow(preddf)


PDO.pdat <- with(df1,
                 data.frame(PDO_3MON_AVG = rep(preddf$PDO_3MON_AVG), # 300*7
                            SOI_3MON_AVG = rep(preddf$SOI_3MON_AVG), # 300*7
                            SRP = median(SRP, na.rm=TRUE),
                            DIN = median(DIN, na.rm=TRUE),
                            #W_temp = median(W_temp, na.rm=TRUE),
                            QLD = median(QLD, na.rm=TRUE),
                            year = rep(2003),
                            DOY = median(DOY)))

PDO.pred <- predict(m2, newdata = PDO.pdat, type = "link")
PDO.pdat <- cbind(PDO.pdat, PDO.pred)
PDO.pdat$PDOgroup <- PDOgroup

names(PDO.pdat)[which(names(PDO.pdat)=='PDO.pred')] <- 'Odour'

# backtransform from tweedie distribution
PDO.pdat$Odour<- exp(PDO.pdat$Odour)


# need to take away places where PDO*SOI combo has not occurred in the data!!
toofar <- exclude.too.far(PDO.pdat$PDO_3MON_AVG, PDO.pdat$SOI_3MON_AVG, df1$PDO_3MON_AVG, df1$SOI_3MON_AVG, dist=0.1)
PDO.pdat$Odour[toofar] <- NA


#reorder factors
PDO.pdat$PDOgroup <- factor(PDO.pdat$PDOgroup, levels = c('-1.5', '0', '1.6'))

PDOplot <- ggplot(PDO.pdat, aes(x = SOI_3MON_AVG, y = Odour, group= PDO_3MON_AVG, col=PDOgroup, 
                                lty=PDOgroup)) +
  theme_bw(base_size=14) +
  theme(legend.position='top') +
  geom_line() +
  scale_y_continuous(limits=c(50, 300))+
  scale_color_manual(name=expression(paste(bold('')~'      PDO')), values = c("blue", "black", "red"))+
  scale_linetype_manual(name=expression(paste(bold('')~'      PDO')), values = c("solid", "solid","longdash")) +
  xlab('SOI') + ylab('Odour (T.O.N)')+
  theme(
    axis.text.x = element_text(color = ifelse(seq(-2.5, 2.5, by = 1) <  -0.5, "blue",
                                              ifelse(seq(-2.5, 2.5, by = 1) > 0, "red", "black")))
  )


PDOplot


saveRDS(PDOplot, "pdo_Odour.rds")  # save the ggplot object


pdo_all<- plot_grid(PDOplot, SOIplot)
pdo_all

ggsave('output/PDO SOI separate infl TW.png', pdo_all, height = 6, width  = 8)

## ---------- W_temp -------------------##
# 
# new_data_temp <- with(df1, expand.grid(W_temp = seq(min(W_temp, na.rm = TRUE), max(W_temp, na.rm = TRUE), length = 200),
#                                      DIN = median(DIN),
#                                      SRP = median(SRP, na.rm=TRUE),
#                                      QLD = median(QLD),
#                                      SOI_3MON_AVG = median(W_temp, na.rm=TRUE), 
#                                      PDO_3MON_AVG = median(PDO_3MON_AVG, na.rm=TRUE),
#                                      year= median(year),
#                                      DOY = median(DOY)))
# 
# temp.pred <- predict(m2, newdata = new_data_temp, type = "terms", se.fit = TRUE)
# 
# whichCols <- grep("W_temp", colnames(temp.pred$fit))
# whichColsSE <- grep("W_temp", colnames(temp.pred$se.fit))
# new_data_temp <- cbind(new_data_temp, Fitted = temp.pred$fit[, whichCols], 
#                        se.Fitted = temp.pred$se.fit[, whichColsSE])
# limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)
# 
# ## make into original limits
# new_data_temp <- with(new_data_temp, transform(new_data_temp, Fittedplus = Fitted + se.Fitted))
# new_data_temp <- with(new_data_temp, transform(new_data_temp, Fittedminus = Fitted - se.Fitted))
# 
# shifttemp <- attr(predict(m2, newdata = new_data_temp, type = "iterms"), "constant")
# temp.pdatnorm <- new_data_temp
# temp.pdatnorm <- with(temp.pdatnorm, transform(temp.pdatnorm, Fitted = Fitted + shifttemp, 
#                                                Fittedplus = Fittedplus + shifttemp, 
#                                                Fittedminus = Fittedminus + shifttemp))
# 
# temp.pdatnorm$Fitted<-exp(temp.pdatnorm$Fitted)
# temp.pdatnorm$Fittedplus<-exp(temp.pdatnorm$Fittedplus)
# temp.pdatnorm$Fittedminus<-exp(temp.pdatnorm$Fittedminus)
# 
# 
# tempquants <- quantile(df1$W_temp, c(.05,.95), na.rm = TRUE)
# 
# tempplot <- ggplot(temp.pdatnorm, aes(x = W_temp, y = Fitted)) +
#   papertheme +
#   #scale_y_continuous(limits=c(30, 125))+
#   annotate("rect", xmin=tempquants[1], xmax=tempquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
#   geom_line() +
#   geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
#               alpha = 0.25, fill = '#165459B2') +  
#   geom_rug(aes(x=W_temp), data = df1, stat = "identity", position = "identity", 
#            sides = "b", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
#   xlab(expression(paste("Water temperature (" , degree*C,")"))) + ylab("Odour (T.O.N)")
# 
# tempplot
# 
# 
# saveRDS(tempplot, "temp_Odour_notemp.rds")  # save the ggplot objec
# 
# 
# # range of control
# max(temp.pdatnorm$Fitted)-min(temp.pdatnorm$Fitted)
# # = 29.8

## ==================================================================================
## Testing effect (NOT RESPONSE) by time plots SOURCE CODE (Wilk et al. 2018)
## ==================================================================================


df1<- df1%>%
  select(date, year, nMonth, DOY, Odour, SRP, DIN, QLD, W_temp, SOI_3MON_AVG, PDO_3MON_AVG)
df1<- na.exclude(df1)

df1$date <- as.Date(df1$date)

df1$seasons <- time2season(df1$date,                # Convert dates to seasons
                          out.fmt = "seasons")

testing1 <- predict(m2, type = 'terms')
testing <- as.data.frame(testing1)

tosum <- grep("TN", colnames(testing))
TNeffect <- rowSums(testing[tosum], na.rm = TRUE)
testing <- testing[,-tosum]
testing$TN <- TNeffect


names(testing) <- c("SRP", "W_temp", "QLD", "SOI_PDO", "time", "TN")
testing$Date <- df1$date 
testing$year <- df1$year
testing$DOY <- df1$DOY
testing$month <- df1$nMonth
testing$seasons <- df1$seasons

testing$month<- as.factor(testing$month)
testing$year<- as.factor(testing$year)

peD <- ggplot(testing, aes(x = month, y = TN)) +
  annotate("rect", ymin = -0.05, ymax = 0.05, 
           xmin = -Inf, xmax = Inf, alpha = 0.3, fill = '#165459B2') +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=55, hjust=1, vjust=1, face = "bold")) +
  #facet_wrap('seasons') +  
  ylab("TN effect") + xlab("Month")
peD


##--------SRP--------------##

peS <- ggplot(testing, aes(x = month, y = SRP)) +
  annotate("rect", ymin = -0.05, ymax = 0.05, 
           xmin = -Inf, xmax = Inf, alpha = 0.3, fill = '#165459B2') +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=55, hjust=1, vjust=1, face = "bold")) +
  #facet_wrap('month') +  
 ylab("SRP Effect") + xlab("Month")
peS


##--------TEMP--------------##

peT <- ggplot(testing, aes(x = month, y = W_temp)) +
  annotate("rect", ymin = -0.05, ymax = 0.05, 
        xmin = -Inf, xmax = Inf, alpha = 0.3, fill = '#165459B2') +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=55, hjust=1, vjust=1, face = "bold")) +
  #facet_wrap('month') +  
  ylab("Temp Effect") + xlab("Month")
peT


##--------QLD--------------##
peQ <- ggplot(testing, aes(x = month, y = QLD)) +
  annotate("rect", ymin = -0.05, ymax = 0.05, 
             xmin = -Inf, xmax = Inf, alpha = 0.3, fill = '#165459B2') +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=55, hjust=1, vjust=1, face = "bold")) +
  #facet_wrap('month') +  
  #geom_text(data = labdatGPP, aes(label = label, x = x, y = y, size = 5), 
  #         show.legend = FALSE) +
  ylab("QLD Effect") + xlab("Month")

peQ




##--------SOI*PDO--------------##
## soi 
pesoi <- ggplot(testing, aes(x = month, y = SOI_PDO)) +
  annotate("rect", ymin = -0.05, ymax = 0.05, 
           xmin = -Inf, xmax = Inf, alpha = 0.3, fill = '#165459B2') +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=55, hjust=1, vjust=1, face = "bold")) +
  #facet_wrap('month') +  
  #geom_text(data = labdatGPP, aes(label = label, x = x, y = y, size = 5), 
  #         show.legend = FALSE) +
  ylab("SOI*PDO Effect") + xlab("Month")

pesoi




p_allEFF<- plot_grid(  peD, peS, peT, peQ, pesoi, ncol = 2, align = "hv")
p_allEFF

ggsave('output/sqrt Odour GAM_partial effects seasonal ALL.png', p_allEFF, height = 12, width  = 12)




## make predictors in the same range as before?

## start with TN

minmax <- function(df, colnames) {
  allmin <- as.data.frame(do.call(cbind, lapply(df[,colnames], min, na.rm = TRUE)))
  names(allmin) <- sapply(names(allmin), function(x) paste("min",x, sep = ""))
  allmax <- as.data.frame(do.call(cbind, lapply(df[,colnames], max, na.rm = TRUE)))
  names(allmax) <- sapply(names(allmax), function(x) paste("max",x, sep = ""))
  summ <- as.data.frame(cbind(allmin, allmax))
    summ
}

df1<- df%>% select(SRP, TN, QLD, W_temp)

minmaxes <- do.call(rbind, lapply(df1, minmax, colnames = c("TN", "SRP", "QLD", 
                                                                 "W_temp")))
rownames(minmaxes) <- NULL



##--------------------------------------------------##
## References
# Simpson, 2014: https://fromthebottomoftheheap.net/2014/05/09/modelling-seasonal-data-with-gam/
# Hayes et al, 2020: https://doi.org/10.1002/lol2.10164
# Simpson 2020 : https://stats.stackexchange.com/questions/471267/plotting-gams-on-response-scale-with-multiple-smooth-and-linear-terms
# Wilk et al., 2018: https://doi.org/10.1029/2018JG004506 (plotting response values)
