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
# 1655 obs

df$DIN <- rowSums(df[,c("NO3", "NH3")], na.rm=TRUE)
df$TN <- rowSums(df[,c("NO3", "NH3", "orgN")], na.rm=TRUE)

## add in PDO and SOI 

clim <- read_csv("data/SOI_PDO data no lag.csv")

df1<- merge(df, clim, by="date")

summary(df1$SOI_3MON_AVG)


##----------------------------------------------------------------------------##
## 3. Check trends, distributions, correlations
##----------------------------------------------------------------------------##

 # basic plot

# df1 %>%
#   dplyr::select(DOC, TP,  W_temp,  orgN, QLD) %>%
#   gather() %>%
#   ggplot(aes(value)) +
#   facet_wrap(~ key, scales = "free") +
#   geom_histogram()

# # ##
#  dfx<-df1%>%
#    select(DOC, SRP, DIN, QLD, W_temp, PDO_3MON_AVG, SOI_3MON_AVG)

# # ## visualize correlations
#  p<- ggpairs(dfx[])  # can see that DOC and turbidity, temperature, Org_N and TP are correlated. Ammonia and conductivity are not. (NH3, NO3, and TP correlated; TP and temp). Org N and TP are similarly correlated with DOC (0.640 and 0.667, respectively)
#  p

# ggsave('output/DOC Correlations.png', p, height = 8, width  = 10)

df_clean <- df1 %>%
  dplyr::ungroup() %>%                    # Remove any previous grouping
  select(year, DOC)    


dx<- df_clean %>%
  dplyr::group_by(year) %>%
  summarize(
    mean_conc = mean(DOC, na.rm = TRUE),
    range_percent = ((max(DOC, na.rm = TRUE) - min(DOC, na.rm = TRUE)) / mean_conc) * 100
  )
dx

mean(dx$range_percent)

# Step 1: Calculate the mean concentration for each year
yearly_means <- df1 %>%
  group_by(year) %>%
  summarize(mean_conc = mean(DOC, na.rm = TRUE), .groups = "drop")

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
## 5. Specify the best fit model (see Fitting script)
##----------------------------------------------------------------------------##


m2 <- gam(DOC ~ s(SRP) +  
            s(DIN)+ 
           # s(W_temp)+
            s(QLD)+
           #ti(DOY, QLD, bs = c("cc", "tp"))+ 
            te(PDO_3MON_AVG,SOI_3MON_AVG)+
            te(year, DOY, bs = c("cr", "cc")),
          knots=list(DOY=c(0, 366.5)),
          select = TRUE,
          data = df1, method = "REML", family = tw(link="log"))


saveRDS(m2, file = "modelDOC_notemp.rds")
# 
 sink("model_summaryDOC.txt")
 summary(m2)
 sink()  


# check model fit...

summary(m2) # Deviance explained = 68, REML = 1338, r2 = 0.64
# all predictors except QLD significant, n=861

k.check(m2)# k-index looks good 

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

#ggsave('output/DOC GAM residuals no temp.png', p2, height = 10, width  = 10)


# ii. Check for autocorrelation
layout(matrix(1:2, ncol = 2))
acf(resid(m2), lag.max = 36, main = "ACF") 
pacf(resid(m2), lag.max = 36, main = "pACF")
layout(1)


#----------------------------------------------------------------------------##
## 7. Plotting on the response scale/fitted values
##----------------------------------------------------------------------------##

# how much does DOC vary on an annual basis?

# Group by year and calculate variation metrics
annual_var <- aggregate(DOC ~ year, data = df1, function(x) {
  if (all(is.na(x))) return(NA)
  max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
})

mean(annual_var$DOC)
#3.363852

min(df1$DOC)
#The average range in DOC concentration each year is 3.36 mg/L
# the range in concentration is 3.11 to 15
## ---------- DIN -------------------##

new_data_DIN<- with(df1, expand.grid(DIN = seq(min(DIN), max(DIN),length = 200),
                                           SRP = median(SRP, na.rm=TRUE),
                                           #W_temp = median(W_temp, na.rm=TRUE),
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
  xlab(expression(paste("DIN"~"("*"mg"~L^{-1}*")"))) + ylab(expression(paste("DOC"~"("*"mg"~L^{-1}*")")))

#DINplot


saveRDS(DINplot, "DIN_DOC_notemp.rds")  # save the ggplot object


# range of control
max(DIN.pdatnorm$Fitted)-min(DIN.pdatnorm$Fitted)
# = NONSIG

### ---------- SRP -------------------##

new_data_SRP <- with(df1, expand.grid(SRP = seq(min(SRP, na.rm = TRUE), max(SRP, na.rm = TRUE), length = 200),
                                     DIN = median(DIN),
                                     #W_temp = median(W_temp, na.rm=TRUE),
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
  xlab(expression(paste("SRP ","(", mu, "g"~L^{-1}*")"))) + ylab(expression(paste("DOC"~"("*"mg"~L^{-1}*")")))

#SRPplot



saveRDS(SRPplot, "SRP_DOC_notemp.rds")  # save the ggplot objec


# range of control
max(SRP.pdatnorm$Fitted)-min(SRP.pdatnorm$Fitted)
# = 4.17
min(SRP.pdatnorm$Fitted) #=6.5
max(SRP.pdatnorm$Fitted) #10.6

## ---------- QLD -------------------##

new_data_QLD <- with(df1, expand.grid(QLD = seq(min(QLD), max(QLD), length = 200),
                                      DIN = median(DIN),
                                      SRP = median(SRP, na.rm=TRUE),
                                     # W_temp = median(W_temp, na.rm=TRUE),
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
  xlab(expression(paste("QLD (",m^3, "/", s,")"))) + ylab(expression(paste("DOC"~"("*"mg"~L^{-1}*")")))

#QLDplot


saveRDS(QLDplot, "QLD_DOC_notemp.rds")  # save the ggplot objec

# range of control
max(QLD.pdatnorm$Fitted)-min(QLD.pdatnorm$Fitted)
# = 0.762


# PUT IT ALL TOGETHER

p_all<- plot_grid(DINplot, SRPplot, QLDplot, ncol = 2)
p_all


#ggsave('output/DOC GAM_Response scale.png', p_all, height = 10, width  = 10)



## ---------- SOI * PDO -------------------##



new_data_SOI_3MON_AVG <- with(df1, expand.grid(SOI_3MON_AVG = seq(min(SOI_3MON_AVG), max(SOI_3MON_AVG), length = 200),
                                               PDO_3MON_AVG = seq(min(PDO_3MON_AVG), max(PDO_3MON_AVG), length = 200),
                                               SRP = median(SRP, na.rm=TRUE),
                                               DIN = median(DIN, na.rm=TRUE),
                                               #W_temp = median(W_temp, na.rm=TRUE),
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
  theme_bw(base_size=14) +
  theme(legend.position='top') +
  geom_raster(aes(fill=DOC)) + # change to turn grey background into nothing
  scale_fill_distiller(palette = "Spectral", direction = -1, na.value='transparent') +
  geom_point(data=df1, aes(x=PDO_3MON_AVG, y=SOI_3MON_AVG, z=NULL)) +
  geom_contour(colour = "black", binwidth = 0.5) +
  scale_x_continuous(expand = c(0, 0)) +  # Remove padding on x-axis
  scale_y_continuous(expand = c(0, 0)) +  # Remove padding on y-axis
  theme(legend.key.width=unit(2,"cm"))+
  xlab("PDO") + ylab("SOI") +
  labs(fill=expression(paste("DOC"~"("*"mg"~L^{-1}*")")))+
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

# range of control
max(SOI_3MON_AVG.pdatnorm$DOC, na.rm = TRUE)-min(SOI_3MON_AVG.pdatnorm$DOC, na.rm = TRUE)
# = 3.6



##-----------------------------------temporal change------------------------------- ####


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
#whichColsSE <- grep("year", colnames(year.pred$se.fit))

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
  geom_raster(aes(fill=DOC)) + # change to turn grey background into nothing
  scale_fill_distiller(palette = "Spectral", direction = -1, na.value='transparent') +
  geom_point(data=df1, aes(x=year, y=DOY, z=NULL)) +
  geom_contour(colour = "black", binwidth = 1) +
  theme(legend.key.width=unit(1.5,"cm")) +
  scale_y_continuous(expand = c(0, 0), breaks = c(1, 91,  182, 274, 335), labels = c("Jan", "Apr", "Jul", "Oct", "Dec")) +
  scale_x_continuous(expand = c(0, 0)) +  # Remove padding on x-axis
  xlab("Year") + ylab("DOY") +
  labs(fill=expression(paste("DOC (mg L"^-1*")   ")))+
  theme(legend.position = "top")


comboplot_time

#ggsave('output/DOC change over time_nolag.png', comboplot_time, height = 6, width  = 8)


saveRDS(comboplot_time, "time_DOC_temp.rds")  # save the ggplot objec



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

names(SOI.pdat)[which(names(SOI.pdat)=='SOI.pred')] <- 'DOC'

# backtransform from tweedie distribution
SOI.pdat$DOC<- exp(SOI.pdat$DOC)

# need to take away places where PDO*SOI combo has not occurred in the data!!
toofar <- exclude.too.far(SOI.pdat$PDO_3MON_AVG, SOI.pdat$SOI_3MON_AVG, df1$PDO_3MON_AVG, df1$SOI_3MON_AVG, dist=0.1)
SOI.pdat$DOC[toofar] <- NA



#reorder factors
SOI.pdat$SOIgroup <- factor(SOI.pdat$SOIgroup, levels = c('-1', '0.3', '1.1'))

SOIplot <- ggplot(SOI.pdat, aes(x = PDO_3MON_AVG, y = DOC, group= SOI_3MON_AVG, col=SOIgroup, 
                                lty=SOIgroup)) +
  theme_bw(base_size=14) +
  theme(legend.position='top') +
  geom_line() +
  scale_y_continuous(limits=c(4.5,7.5))+
  scale_color_manual(name=expression(paste(bold('')~'      SOI')), values = c("red", "black", "blue"))+
  scale_linetype_manual(name=expression(paste(bold('')~'      SOI')), values = c("solid", "solid","longdash")) +
  xlab('PDO') + ylab(expression(paste("DOC (mg L"^-1*")   ")))+
theme(
  axis.text.x = element_text(color = ifelse(seq(-2.5, 2.5, by = 1) < 0, "blue",
                                            ifelse(seq(-2.5, 2.5, by = 1) > 0.5, "red", "black")))
)
  
## '#e66101','#fdb863','#b2abd2','#5e3c99' (order red:light:dark)

SOIplot


saveRDS(SOIplot, "soi_DOC.rds")  # save the ggplot object


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

names(PDO.pdat)[which(names(PDO.pdat)=='PDO.pred')] <- 'DOC'

# backtransform from tweedie distribution
PDO.pdat$DOC<- exp(PDO.pdat$DOC)


# need to take away places where PDO*SOI combo has not occurred in the data!!
toofar <- exclude.too.far(PDO.pdat$PDO_3MON_AVG, PDO.pdat$SOI_3MON_AVG, df1$PDO_3MON_AVG, df1$SOI_3MON_AVG, dist=0.1)
PDO.pdat$DOC[toofar] <- NA


#reorder factors
PDO.pdat$PDOgroup <- factor(PDO.pdat$PDOgroup, levels = c('-1.5', '0', '1.6'))

PDOplot <- ggplot(PDO.pdat, aes(x = SOI_3MON_AVG, y = DOC, group= PDO_3MON_AVG, col=PDOgroup, 
                                lty=PDOgroup)) +
  theme_bw(base_size=14) +
  theme(legend.position='top') +
  geom_line() +
  scale_y_continuous(limits=c(4.5,7.5))+
  scale_color_manual(name=expression(paste(bold('')~'      PDO')), values = c("blue", "black", "red"))+
  scale_linetype_manual(name=expression(paste(bold('')~'      PDO')), values = c("solid", "solid","longdash")) +
  xlab('SOI') + ylab(expression(paste("DOC (mg L"^-1*")   ")))+
  theme(
    axis.text.x = element_text(color = ifelse(seq(-2.5, 2.5, by = 1) <  -0.5, "blue",
                                              ifelse(seq(-2.5, 2.5, by = 1) > 0, "red", "black")))
  )

PDOplot


saveRDS(PDOplot, "pdo_DOC.rds")  # save the ggplot object


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
#   annotate("rect", xmin=tempquants[1], xmax=tempquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
#   geom_line() +
#   geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
#               alpha = 0.25, fill = '#165459B2') +  
#    geom_rug(aes(x=W_temp), data = df1, stat = "identity", position = "identity", 
#            sides = "b", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
#   xlab(expression(paste("Water temperature (" , degree*C,")"))) +ylab(expression(paste("DOC"~"("*"mg"~L^{-1}*")")))
# 
# #tempplot
# 
# 
# 
# saveRDS(tempplot, "temp_DOC_nolag.rds")  # save the ggplot objec
# 
# # range of control
# max(temp.pdatnorm$Fitted)-min(temp.pdatnorm$Fitted)
# # = 2.1
# 

##--------------------------------------------------##
## References
# Simpson, 2014: https://fromthebottomoftheheap.net/2014/05/09/modelling-seasonal-data-with-gam/
# Hayes et al, 2020: https://doi.org/10.1002/lol2.10164
# Simpson 2020 : https://stats.stackexchange.com/questions/471267/plotting-gams-on-response-scale-with-multiple-smooth-and-linear-terms
# Wilk et al., 2018: https://doi.org/10.1029/2018JG004506 (plotting response values)
