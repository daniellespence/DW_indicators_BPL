 ###############################################################################
####### BPWTP Data -- GAMs ##########
####### Danielle Spence ##########
####### Created 5/4/2023 ########
###############################################################################
### Clear memory
rm(list = ls())

library(pacman)
p_load(tidyverse, ggplot2, mgcv, gratia, readr, GGally, dplyr,  mgcViz, hydroTSM, ggridges, lubridate, plyr,
        cowplot, rsoi, rpdo, tibble,viridis, install = TRUE)

setwd("C:/Users/danis/OneDrive/R/DW Indicators")

papertheme <- theme_bw(base_size=14, base_family = 'Arial') 

##----------------------------------------------------------------------------##
## 1. Read in data 
##----------------------------------------------------------------------------##

df <- read_csv("data/bpgamdataCLEAN_Conductivity.csv")

df$Date<- as.Date(df$Date, "%m/%d/%Y")

# add in Year and nMonth for numeric month and a proper Date class

df <- mutate(df,
             year = as.numeric(format(Date,'%Y')),
             DOY = as.numeric(format(Date,'%j')),
             nMonth = as.numeric(format(Date,'%m')),
             week = as.numeric(format(Date, '%W')))%>% 
  filter(year %in% c(1990:2022))


## order the data
df <- df[with(df, order(Date)),]


# rename for ease of use
df <- df%>%
  dplyr::rename(
    date = "Date",
    Conductivity = "Conductivity",
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



df$DIN <- rowSums(df[,c("NO3", "NH3")], na.rm=TRUE)

## add in PDO and SOI 

clim <- read_csv("data/SOI_PDO data.csv")

df1<- merge(df, clim, by="date")

## remove 0s/missing values from analysis....

#Remove rows without Chla data
df1 <- df1[complete.cases(df$Conductivity),]

# 1330 obs


##----------------------------------------------------------------------------##
## 3. Check trends, distributions, correlations
##----------------------------------------------------------------------------##

# basic plot

df1 %>%
  ggplot(aes(x = year, y = Conductivity)) +
  geom_line()+
  ggtitle("Conductivity")

df1 %>%
  dplyr::select(Conductivity, SRP,  W_temp,  DIN, QLD) %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()



## correlaitons
dfx<-df1%>%
  select(Conductivity, SRP, DIN, QLD,  W_temp, PDO_3MON_AVG, SOI_3MON_AVG)
## visualize correlations
p<- ggpairs(dfx[])  # can see that Conductivity and turbidity, temperature, Org_N and TP are correlated. Ammonia and conductivity are not. (NH3, NO3, and TP correlated; TP and temp). Org N and TP are similarly correlated with Conductivity (0.640 and 0.667, respectively)
p 


ggsave('output/Conductivity Correlations.png', p, height = 8, width  = 10)


##----------------------------------------------------------------------------##
## 5. Run the full model
##----------------------------------------------------------------------------##

m <- gam(Conductivity ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
           select = TRUE,
           data = df1, method = "REML", family = gaussian, 
           na.action = na.exclude)


m1 <- gam(Conductivity ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
           select = TRUE,
           data = df1, method = "REML",family = tw(), 
           na.action = na.exclude)

m2 <- gam(Conductivity ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
           select = TRUE,
           data = df1, method = "REML",  family = Gamma(link = "log"),
           na.action = na.exclude)

m3 <- gam(Conductivity ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
           select = TRUE,
           data = df1, method = "REML",family = scat(), 
           na.action = na.exclude)

m4 <- gam(Conductivity ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
           select = TRUE,
           data = df1, method = "REML", family = nb(), 
           na.action = na.exclude)

m5 <- gam(Conductivity ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
           select = TRUE,
           data = df1, method = "REML", family = poisson,
           na.action = na.exclude)

AIC(m, m1, m2, m3, m4, m5) # m1 (tw) or m2 (gamma family) is best fit



##--- check which variables should stay/be removed

nullvars <- c('Conductivity', "SRP", "DIN", "QLD", "W_temp", "PDO_3MON_AVG", "SOI_3MON_AVG", "year", "DOY")
mumdat <- df1[complete.cases(df1[,nullvars]),] # want to AIC compare, and NA diff
#   between vars
mummod <- update(m2, data=mumdat, na.action="na.fail")
mummoddin <- update(mummod, .~.-s(DIN))
mummodsrp <- update(mummod, .~.-s(SRP))
mummodwtemp <- update(mummod, .~.-s(W_temp))
mummodQLD <- update(mummod, .~.-s(QLD))
mummodTE <- update(mummod, .~.-te(PDO_3MON_AVG,SOI_3MON_AVG))

AIC(mummod, mummoddin,mummodsrp, mummodwtemp, mummodQLD, mummodTE) # 

# dropping any predictor doesnt improve fit. keep as is.


summary(m2) # Deviance explained = 69%, REML = 6455, r2 = 0.652, n=1138
# SRP,  water temp, PDO*SOI, and DOY*YEar all sig
#QLD not sig

k.check(m2)# k-index looks good 

p1<- draw(m2, residuals = TRUE)& theme_bw() 
p1

#ggsave('output/Conductivity GAM_SQRT SOI PDO interact.png', p1, height = 10, width  = 10)


##----------------------------------------------------------------------------##
## 6. Assess model fit & autocorrelation
##----------------------------------------------------------------------------##


## i. Appraise model fit
p2<-appraise(m2, point_col = 'steelblue', point_alpha = 0.5, n_bins = 'fd') & 
  theme(plot.tag = element_text(face = 'bold')) & theme_bw()
p2


ggsave('output/Conductivity GAM residuals.png', p2, height = 10, width  = 10)


# ii. Check for autocorrelation
layout(matrix(1:2, ncol = 2))
acf(resid(m2), lag.max = 36, main = "ACF") 
pacf(resid(m2), lag.max = 36, main = "pACF")
layout(1)

##----------------------------------------------------------------------------##
## 6. Control for autocorrelation
##----------------------------------------------------------------------------##


ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")

m<- gamm(Conductivity ~ s(SRP) +  
           s(DIN)+ 
           s(W_temp)+
           s(QLD)+
           te(SOI_3MON_AVG, PDO_3MON_AVG)+
           te(year, DOY, bs = c("cr", "cc")),
         select = TRUE,
         data = df1,method = "REML")

## AR(1)
m1 <- gamm(Conductivity ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(SOI_3MON_AVG, PDO_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
           select = TRUE,
           data = df1,correlation = corARMA(form = ~ 1|year, p = 1),
           control = ctrl, method = "REML")
## AR(2)
m2 <- gamm(Conductivity ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(SOI_3MON_AVG, PDO_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
           select = TRUE,
           data = df1, correlation = corARMA(form = ~ 1|year, p = 2),
           control = ctrl, method = "REML")

## AR(3)
m3 <- gamm(Conductivity ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(SOI_3MON_AVG, PDO_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
           select = TRUE,
           data = df1, correlation = corARMA(form = ~ 1|year, p = 3),
           control = ctrl, method = "REML")

anova(m$lme, m1$lme, m2$lme) ## shows that m1 is best


layout(matrix(1:2, ncol = 2))
res <- resid(m1$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(2) errors")
acf(res, lag.max = 36, main = "pACF- AR(2) errors")
layout(1) # this is better


summary(m1$gam) 

k.check(m1$gam) 

## i. Appraise model fit
appraise(m1$gam, point_col = 'steelblue', point_alpha = 0.5, n_bins = 'fd') & 
  theme(plot.tag = element_text(face = 'bold')) & theme_bw()

p1<- draw(m1$gam, residuals = TRUE,  select = smooths(m1$gam)[1:5])& theme_bw()
p1

#ggsave('output/Conductivity GAM_ autocorrelation no interaction SOI.png', p1, height = 10, width  = 10)



#----------------------------------------------------------------------------##
## 7. Plotting on the response scale/fitted values
##----------------------------------------------------------------------------##

## ---------- DIN -------------------##

new_data_DIN <- with(df1, expand.grid(DIN = seq(min(DIN), max(DIN),length = 200),
                                           SRP = median(SRP, na.rm=TRUE),
                                           W_temp = median(W_temp, na.rm=TRUE),
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

#DIN.pdatnorm <- merge(DIN.pdatnorm, minmaxes)

#overs <- with(DIN.pdatnorm, which(DIN < min(DIN) | DIN > max(DIN)))
#DIN.pdatnorm <- DIN.pdatnorm[-overs,]

#labdatoxy <- data.frame(x = c(2.5, 11), y = c(meanpH + 0.03, 8.4), label = c("mean pH", 'mean oxygen'))

DINquants <- quantile(df1$DIN, c(.05,.95), na.rm = TRUE)

DINplot <- ggplot(DIN.pdatnorm, aes(x = DIN, y = Fitted)) +
  papertheme +
  annotate("rect", xmin=DINquants[1], xmax=DINquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
  geom_line() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25, fill = 'gray60') +  
  #geom2ext(data = labdatoxy, aes(label = label, x = x, y = y, size = 5), 
  #         show.legend = FALSE) +
  #geom_abline(slope = 0, intercept = meanConductivity, linetype="dotted") +
  # geom_vline(xintercept = meanDIN, linetype="dotted") +
  geom_rug(aes(x=DIN), data = df1, stat = "identity", position = "identity", 
           sides = "b", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
  xlab(expression(paste("DIN"~"("*"mg"~L^{-1}*")"))) + ylab(expression(paste("Conductivity"~"("*"mg"~L^{-1}*")")))

DINplot


### ---------- SRP -------------------##

new_data_SRP <- with(df1, expand.grid(SRP = seq(min(SRP, na.rm = TRUE), max(SRP, na.rm = TRUE), length = 200),
                                     DIN = median(DIN),
                                     W_temp = median(W_temp, na.rm=TRUE),
                                     QLD = median(QLD),
                                     SOI_3MON_AVG = median(W_temp, na.rm=TRUE), 
                                     PDO_3MON_AVG = median(PDO_3MON_AVG, na.rm=TRUE),
                                     year= median(year),
                                     DOY = median(DOY)))


SRP.pred <- predict(m2, newdata = new_data_SRP, type = "terms", se.fit = TRUE)


# backtransform fitted data...

#SRP.pred$fit<-sqrt1pexp(SRP.pred$fit)
#SRP.pred$se.fit<-sqrt1pexp(SRP.pred$fit)


whichCols <- grep("SRP", colnames(SRP.pred$fit))
whichColsSE <- grep("SRP", colnames(SRP.pred$se.fit))
new_data_SRP <- cbind(new_data_SRP, Fitted = SRP.pred$fit[, whichCols], 
                      se.Fitted = SRP.pred$se.fit[, whichColsSE])

#new_data_SRP$srp<-sqrt1pexp(new_data_SRP$SRP)

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

#SRP.pdatnorm <- merge(SRP.pdatnorm, minmaxes)

#overs <- with(SRP.pdatnorm, which(SRP < min(SRP) | SRP > max(SRP)))
#SRP.pdatnorm <- SRP.pdatnorm[-overs,]

#labdatoxy <- data.frame(x = c(2.5, 11), y = c(meanpH + 0.03, 8.4), label = c("mean pH", 'mean oxygen'))

SRPquants <- quantile(df1$SRP, c(.05,.95), na.rm = TRUE)


SRPplot <- ggplot(SRP.pdatnorm, aes(x = SRP, y = Fitted)) +
  papertheme+
 annotate("rect", xmin=SRPquants[1], xmax=SRPquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
  geom_line() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25, fill = '#165459B2') +  
  #geom2ext(data = labdatoxy, aes(label = label, x = x, y = y, size = 5), 
  #         show.legend = FALSE) +
  #geom_abline(slope = 0, intercept = meanConductivity, linetype="dotted") +
  #geom_vline(xintercept = meanSRP, linetype="dotted") +
  geom_rug(aes(x=SRP), data = df1, stat = "identity", position = "identity", 
           sides = "b", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
  xlab(expression(paste("SRP ","(", mu, "g"~L^{-1}*")"))) + ylab(expression(paste("Conductivity"~"("*"mg"~L^{-1}*")")))

SRPplot

## ---------- W_temp -------------------##

new_data_temp <- with(df1, expand.grid(W_temp = seq(min(W_temp, na.rm = TRUE), max(W_temp, na.rm = TRUE), length = 200),
                                     DIN = median(DIN),
                                     SRP = median(SRP, na.rm=TRUE),
                                     QLD = median(QLD),
                                     SOI_3MON_AVG = median(W_temp, na.rm=TRUE), 
                                     PDO_3MON_AVG = median(PDO_3MON_AVG, na.rm=TRUE),
                                     year= median(year),
                                     DOY = median(DOY)))

temp.pred <- predict(m2, newdata = new_data_temp, type = "terms", se.fit = TRUE)

whichCols <- grep("W_temp", colnames(temp.pred$fit))
whichColsSE <- grep("W_temp", colnames(temp.pred$se.fit))
new_data_temp <- cbind(new_data_temp, Fitted = temp.pred$fit[, whichCols], 
                       se.Fitted = temp.pred$se.fit[, whichColsSE])
limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

## make into original limits
new_data_temp <- with(new_data_temp, transform(new_data_temp, Fittedplus = Fitted + se.Fitted))
new_data_temp <- with(new_data_temp, transform(new_data_temp, Fittedminus = Fitted - se.Fitted))

shifttemp <- attr(predict(m2, newdata = new_data_temp, type = "iterms"), "constant")
temp.pdatnorm <- new_data_temp
temp.pdatnorm <- with(temp.pdatnorm, transform(temp.pdatnorm, Fitted = Fitted + shifttemp, 
                                               Fittedplus = Fittedplus + shifttemp, 
                                               Fittedminus = Fittedminus + shifttemp))

temp.pdatnorm$Fitted<-exp(temp.pdatnorm$Fitted)
temp.pdatnorm$Fittedplus<-exp(temp.pdatnorm$Fittedplus)
temp.pdatnorm$Fittedminus<-exp(temp.pdatnorm$Fittedminus)

#temp.pdatnorm <- merge(temp.pdatnorm, minmaxes)

#overs <- with(temp.pdatnorm, which(temp < min(temp) | temp > max(temp)))
#temp.pdatnorm <- temp.pdatnorm[-overs,]

#labdatoxy <- data.frame(x = c(2.5, 11), y = c(meanpH + 0.03, 8.4), label = c("mean pH", 'mean oxygen'))

tempquants <- quantile(df1$W_temp, c(.05,.95), na.rm = TRUE)

tempplot <- ggplot(temp.pdatnorm, aes(x = W_temp, y = Fitted)) +
  papertheme +
  annotate("rect", xmin=tempquants[1], xmax=tempquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
  geom_line() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25, fill = '#165459B2') +  
  #geom2ext(data = labdatoxy, aes(label = label, x = x, y = y, size = 5), 
  #         show.legend = FALSE) +
  #geom_abline(slope = 0, intercept = meanConductivity, linetype="dotted") +
  #geom_vline(xintercept = meanW_temp, linetype="dotted") +
  geom_rug(aes(x=W_temp), data = df1, stat = "identity", position = "identity", 
           sides = "b", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
  xlab(expression(paste("Water temperature (" , degree*C,")"))) +ylab(expression(paste("Conductivity"~"("*"mg"~L^{-1}*")")))

tempplot

## ---------- QLD -------------------##

new_data_QLD <- with(df1, expand.grid(QLD = seq(min(QLD), max(QLD), length = 200),
                                      DIN = median(DIN),
                                      SRP = median(SRP, na.rm=TRUE),
                                      W_temp = median(W_temp, na.rm=TRUE),
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
#QLD.pdatnorm <- merge(QLD.pdatnorm, minmaxes)

#overs <- with(QLD.pdatnorm, which(QLD < min(QLD) | QLD > max(QLD)))
#QLD.pdatnorm <- QLD.pdatnorm[-overs,]

#labdatoxy <- data.frame(x = c(2.5, 11), y = c(meanpH + 0.03, 8.4), label = c("mean pH", 'mean oxygen'))

QLDquants <- quantile(df1$QLD, c(.05,.95), na.rm = TRUE)

QLDplot <- ggplot(QLD.pdatnorm, aes(x = QLD, y = Fitted)) +
  papertheme +
  annotate("rect", xmin=QLDquants[1], xmax=QLDquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
  geom_line() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25, fill = '#165459B2') +  
 geom_rug(aes(x=QLD), data = df1, stat = "identity", position = "identity", 
           sides = "b", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
  xlab(expression(paste("QLD (",m^3, "/", s,")"))) + ylab(expression(paste("Conductivity"~"("*"mg"~L^{-1}*")")))

QLDplot


## ---------- SOI * PDO -------------------##



new_data_SOI_3MON_AVG <- with(df1, expand.grid(SOI_3MON_AVG = seq(min(SOI_3MON_AVG), max(SOI_3MON_AVG), length = 200),
                                               PDO_3MON_AVG = seq(min(PDO_3MON_AVG), max(PDO_3MON_AVG), length = 200),
                                               SRP = median(SRP),
                                               DIN = median(DIN, na.rm=TRUE),
                                               W_temp = median(W_temp, na.rm=TRUE),
                                               QLD = median(QLD, na.rm=TRUE),
                                               
                                               year= median(year),
                                               DOY = median(DOY)))

SOI_3MON_AVG.pred <- predict(m2, newdata = new_data_SOI_3MON_AVG, type = "terms")

whichCols <- grep("PDO_3MON_AVG,SOI_3MON_AVG", colnames(SOI_3MON_AVG.pred))
#whichColsSE <- grep("SOI_3MON_AVG", colnames(SOI_3MON_AVG.pred$se.fit))

new_data_SOI_3MON_AVG <- cbind(new_data_SOI_3MON_AVG, Fitted = SOI_3MON_AVG.pred[, whichCols])

shiftcomb <- attr(SOI_3MON_AVG.pred, "constant")
SOI_3MON_AVG.pdatnorm <- new_data_SOI_3MON_AVG
SOI_3MON_AVG.pdatnorm <- with(SOI_3MON_AVG.pdatnorm, transform(SOI_3MON_AVG.pdatnorm, Fitted = Fitted + shiftcomb))

toofar <- exclude.too.far(SOI_3MON_AVG.pdatnorm$PDO_3MON_AVG, SOI_3MON_AVG.pdatnorm$SOI_3MON_AVG, df1$PDO_3MON_AVG, df1$SOI_3MON_AVG, dist=0.1)
SOI_3MON_AVG.pdatnorm$Conductivity <- SOI_3MON_AVG.pdatnorm$Fitted
SOI_3MON_AVG.pdatnorm$Conductivity[toofar] <- NA


SOI_3MON_AVG.pdatnorm$Conductivity<-exp(SOI_3MON_AVG.pdatnorm$Conductivity)

names(new_data_SOI_3MON_AVG)[which(names(new_data_SOI_3MON_AVG)=='SOI_3MON_AVG.pred')] <- 'Conductivity'

comboplot <- ggplot(SOI_3MON_AVG.pdatnorm, aes(x = PDO_3MON_AVG, y = SOI_3MON_AVG, z=Conductivity)) + 
  theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top') +
  geom_raster(aes(fill=Conductivity)) + # change to turn grey background into nothing
  scale_fill_distiller(palette = "Spectral", direction = -1, na.value='transparent') +
  # scale_fill_viridis(na.value='transparent') +
  geom_point(data=df1, aes(x=PDO_3MON_AVG, y=SOI_3MON_AVG, z=NULL)) +
  geom_contour(colour = "black", binwidth = 25) +
  theme(legend.key.width=unit(2,"cm"))+
  xlab("PDO") + ylab("SOI") +
  labs(fill=expression(paste("Conductivity"~"("*"mg"~L^{-1}*")")))+
  theme(legend.position = "top")


comboplot


ggsave('output/SOI PDO INTERACTION Response scale Conductivity.png', comboplot, height = 8, width  = 10)



p_all<- plot_grid(SRPplot, DINplot, tempplot, QLDplot, ncol = 2)
p_all


ggsave('output/Conductivity GAM_Response scale.png', p_all, height = 10, width  = 10)



## ==================================================================================
## Testing effect (NOT RESPONSE) by time plots SOURCE CODE (Wilk et al. 2018)
## ==================================================================================


df1<- df1%>%
  select(date, year, nMonth, DOY, Conductivity, SRP, DIN, QLD, W_temp, SOI_3MON_AVG, PDO_3MON_AVG)
df1<- na.exclude(df1)

testing1 <- predict(m2, type = 'terms')
testing <- as.data.frame(testing1)

tosum <- grep("DIN", colnames(testing))
DINeffect <- rowSums(testing[tosum], na.rm = TRUE)
testing <- testing[,-tosum]
testing$DIN <- DINeffect


names(testing) <- c("SRP", "W_temp", "QLD", "SOI_PDO", "time", "DIN")
testing$Date <- df1$date 
testing$Year <- df1$year
testing$DOY <- df1$DOY
testing$Month <- df1$nMonth

testing$Month<- as.factor(testing$Month)
testing$Year<- as.factor(testing$Year)

peD <- ggplot(testing, aes(x = Month, y = DIN)) +
  annotate("rect", ymin = -0.05, ymax = 0.05, 
           xmin = -Inf, xmax = Inf, alpha = 0.3, fill = '#165459B2') +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=55, hjust=1, vjust=1, face = "bold")) +
  #facet_wrap('seasons') +  
  ylab("DIN effect") + xlab("Month")
peD


##--------SRP--------------##

peS <- ggplot(testing, aes(x = Month, y = SRP)) +
  annotate("rect", ymin = -0.05, ymax = 0.05, 
           xmin = -Inf, xmax = Inf, alpha = 0.3, fill = '#165459B2') +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=55, hjust=1, vjust=1, face = "bold")) +
  #facet_wrap('Month') +  
 ylab("SRP Effect") + xlab("Month")
peS


##--------TEMP--------------##

peT <- ggplot(testing, aes(x = Month, y = W_temp)) +
  annotate("rect", ymin = -0.05, ymax = 0.05, 
        xmin = -Inf, xmax = Inf, alpha = 0.3, fill = '#165459B2') +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=55, hjust=1, vjust=1, face = "bold")) +
  #facet_wrap('Month') +  
  ylab("Temp Effect") + xlab("Month")
peT


##--------QLD--------------##
peQ <- ggplot(testing, aes(x = Month, y = QLD)) +
  annotate("rect", ymin = -0.05, ymax = 0.05, 
             xmin = -Inf, xmax = Inf, alpha = 0.3, fill = '#165459B2') +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=55, hjust=1, vjust=1, face = "bold")) +
  #facet_wrap('Month') +  
  #geom2ext(data = labdatGPP, aes(label = label, x = x, y = y, size = 5), 
  #         show.legend = FALSE) +
  ylab("QLD Effect") + xlab("Month")

peQ




##--------SOI*PDO--------------##
## soi 
pesoi <- ggplot(testing, aes(x = Month, y = SOI_PDO)) +
  annotate("rect", ymin = -0.05, ymax = 0.05, 
           xmin = -Inf, xmax = Inf, alpha = 0.3, fill = '#165459B2') +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=55, hjust=1, vjust=1, face = "bold")) +
  #facet_wrap('Month') +  
  #geom2ext(data = labdatGPP, aes(label = label, x = x, y = y, size = 5), 
  #         show.legend = FALSE) +
  ylab("SOI*PDO Effect") + xlab("Month")

pesoi




p_allEFF<- plot_grid(peD, peS, peT, peQ, pesoi, ncol = 2, align = "hv")
p_allEFF

ggsave('output/Conductivity GAM_partial effects seasonal.png', p_allEFF, height = 12, width  = 12)




## make predictors in the same range as before?

## start with DIN

minmax <- function(df, colnames) {
  allmin <- as.data.frame(do.call(cbind, lapply(df[,colnames], min, na.rm = TRUE)))
  names(allmin) <- sapply(names(allmin), function(x) paste("min",x, sep = ""))
  allmax <- as.data.frame(do.call(cbind, lapply(df[,colnames], max, na.rm = TRUE)))
  names(allmax) <- sapply(names(allmax), function(x) paste("max",x, sep = ""))
  summ <- as.data.frame(cbind(allmin, allmax))
    summ
}

df1<- df%>% select(SRP, DIN, QLD, W_temp)

minmaxes <- do.call(rbind, lapply(df1, minmax, colnames = c("DIN", "SRP", "QLD", 
                                                                 "W_temp")))
rownames(minmaxes) <- NULL





##--------------------------------------------------##
## References
# Simpson, 2014: https://fromthebottomoftheheap.net/2014/05/09/modelling-seasonal-data-with-gam/
# Hayes et al, 2020: https://doi.org/10.1002/lol2.10164
# Simpson 2020 : https://stats.stackexchange.com/questions/471267/plotting-gams-on-response-scale-with-multiple-smooth-and-linear-terms
# Wilk et al., 2018: https://doi.org/10.1029/2018JG004506 (plotting response values)
