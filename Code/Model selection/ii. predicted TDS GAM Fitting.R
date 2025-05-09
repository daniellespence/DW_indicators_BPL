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
# 1039 obs


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


#Remove rows without Conductivity data
df <- df[complete.cases(df$Conductivity),]
# 1022 obs


df$DIN <- rowSums(df[,c("NO3", "NH3")], na.rm=TRUE)


## add in PDO and SOI 

clim <- read_csv("data/SOI_PDO data no lag.csv")

df1<- merge(df, clim, by="date")


## add TDS data

tds <- read_csv("data/bpgamdataCLEAN_TDS.csv")

# 305 obs

tds$Date<- as.Date(tds$Date, "%m/%d/%Y")

# add in Year and nMonth for numeric month and a proper Date class


## order the data
tds <- tds[with(tds, order(Date)),]


# rename for ease of use
tds <- tds%>%
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
    orgN = "Org_N_mg.L",
  )

#Remove rows without Chla data
tds <- tds[complete.cases(tds$TDS),]

tds<- tds %>% select(date, TDS)
# 284 obs


df1<- left_join(df1, tds, by="date")

# Fit linear model
lm_model <- lm(TDS ~ Conductivity, data = df1)
summary(lm_model)


# Predict TDS where it's missing
df1 <- df1 %>%
  mutate(Predicted_TDS = ifelse(is.na(TDS), predict(lm_model, newdata = df1), TDS))



##----------------------------------------------------------------------------##
## 3. Check trends, distributions, correlations
##----------------------------------------------------------------------------##

# basic plot

df1 %>%
  ggplot(aes(x = year, y = TDS)) +
  geom_line()+
  ggtitle("TDS")

df1 %>%
  dplyr::select(TDS, SRP,  W_temp,  DIN, pH) %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()


## transform positively skewed data (TDS, SRP, DIN)

df1$srp_sqrt<- sqrt(df1$SRP)
df1$din_sqrt<- sqrt(df1$DIN)
df1$tds_sqrt<- sqrt(df1$TDS)

df1 %>%
  dplyr::select(tds_sqrt, srp_sqrt,  W_temp,  din_sqrt) %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()


## add ph??? re: Zepernick et al. 2024
dfx<-df1%>%
  select(TDS, SRP, DIN, QLD, QWS, QBP, W_temp, PDO_3MON_AVG, SOI_3MON_AVG, pH)
## visualize correlations
p<- ggpairs(dfx[])  # can see that TDS and turbidity, temperature, Org_N and TP are correlated. Ammonia and conductivity are not. (NH3, NO3, and TP correlated; TP and temp). Org N and TP are similarly correlated with TDS (0.640 and 0.667, respectively)
p 


ggsave('output/TDS Correlations.png', p, height = 8, width  = 10)


##----------------------------------------------------------------------------##
## 5. Run the full model
##----------------------------------------------------------------------------##

m <- gam(Predicted_TDS ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
           select = TRUE,
           data = df1, method = "REML", family = gaussian, 
           na.action = na.exclude)


m1 <- gam(Predicted_TDS ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
           select = TRUE,
           data = df1, method = "REML",family = tw(), 
           na.action = na.exclude)

m2 <- gam(Predicted_TDS ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
           select = TRUE,
           data = df1, method = "REML",  family = Gamma(link = "log"),
           na.action = na.exclude)

m3 <- gam(Predicted_TDS ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
           select = TRUE,
           data = df1, method = "REML",family = scat(), 
           na.action = na.exclude)

m4 <- gam(Predicted_TDS ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
           select = TRUE,
           data = df1, method = "REML", family = nb(), 
           na.action = na.exclude)

m5 <- gam(Predicted_TDS ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
           select = TRUE,
           data = df1, method = "REML", family = poisson,
           na.action = na.exclude)

AIC(m, m1, m2, m3, m4) #  m2 (gamma family) is best fit



##--- check which variables should stay/be removed

nullvars <- c('TDS', "SRP", "DIN", "QLD", "W_temp", "PDO_3MON_AVG", "SOI_3MON_AVG", "year", "DOY")
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


summary(m2) # Deviance explained = 80%, REML = 1107, r2 = 0.757, n=198
# SRP, QLD, water temp, PDO*SOI, and DOY*YEar all sig

k.check(m2)# k-index looks good 

p1<- draw(m2, residuals = TRUE)& theme_bw() 
p1

#ggsave('output/TDS GAM_SQRT SOI PDO interact.png', p1, height = 10, width  = 10)


##----------------------------------------------------------------------------##
## 6. Assess model fit & autocorrelation
##----------------------------------------------------------------------------##


## i. Appraise model fit
p2<-appraise(m2, point_col = 'steelblue', point_alpha = 0.5, n_bins = 'fd') & 
  theme(plot.tag = element_text(face = 'bold')) & theme_bw()
p2


ggsave('output/TDS GAM residuals.png', p2, height = 10, width  = 10)


# ii. Check for autocorrelation
layout(matrix(1:2, ncol = 2))
acf(resid(m2), lag.max = 36, main = "ACF") 
pacf(resid(m2), lag.max = 36, main = "pACF")
layout(1)

##----------------------------------------------------------------------------##
## 6. Control for autocorrelation
##----------------------------------------------------------------------------##


ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")

m<- gamm(TDS ~ s(SRP) +  
           s(DIN)+ 
           s(W_temp)+
           s(QLD)+
           te(SOI_3MON_AVG, PDO_3MON_AVG)+
           te(year, DOY, bs = c("cr", "cc")),
         select = TRUE,
         data = df1,method = "REML")

## AR(1)
m1 <- gamm(TDS ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(SOI_3MON_AVG, PDO_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
           select = TRUE,
           data = df1,correlation = corARMA(form = ~ 1|year, p = 1),
           control = ctrl, method = "REML")
## AR(2)
m2 <- gamm(TDS ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(SOI_3MON_AVG, PDO_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
           select = TRUE,
           data = df1, correlation = corARMA(form = ~ 1|year, p = 2),
           control = ctrl, method = "REML")

## AR(3)
m3 <- gamm(TDS ~ s(SRP) +  
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

#ggsave('output/TDS GAM_ autocorrelation no interaction SOI.png', p1, height = 10, width  = 10)


##----------------------------------------------------------------------------##
## 5. Check which fit is best for time
##----------------------------------------------------------------------------##


# year not added

m <- gam(Predicted_TDS ~ s(SRP) +  
           s(DIN)+ 
           #s(W_temp)+
           s(QLD)+
           te(PDO_3MON_AVG,SOI_3MON_AVG),
         #te(year, DOY, bs = c("cr", "cc")),
         select = TRUE,
         data = df1, method = "REML", Gamma(link = "log"))

m1 <- gam(Predicted_TDS ~ s(SRP) +  
            s(DIN)+ 
            #s(W_temp)+
            s(QLD)+
            te(PDO_3MON_AVG,SOI_3MON_AVG)+
            s(year,  bs = "cr")+
            s(DOY, bs = "cc"),
          knots=list(DOY=c(0, 366.5)),
          select = TRUE,
          data = df1, method = "REML", family = Gamma(link = "log"))

m2 <- gam(Predicted_TDS ~ s(SRP) +  
            s(DIN)+ 
            #s(W_temp)+
            s(QLD)+
            te(PDO_3MON_AVG,SOI_3MON_AVG)+
            te(year, DOY, bs = c("cr", "cc"), k = c(10,12)),
          knots=list(DOY=c(0, 366.5)),
          select = TRUE,
          data = df1, method = "REML", family = tw(link = "log"))

AIC(m, m1, m2) #m2 
BIC(m, m1, m2)

summary(m2) # Deviance explained =73%, REML = 3504, r2 = 0.463

k.check(m1)# k-index looks good 

p1<- draw(m2,  residuals = TRUE)& theme_bw() 
p1


# model 1 is a better fit in terms of AIC and BIC.

##--------------------------------------------------##
## References
# Simpson, 2014: https://fromthebottomoftheheap.net/2014/05/09/modelling-seasonal-data-with-gam/
# Hayes et al, 2020: https://doi.org/10.1002/lol2.10164
# Simpson 2020 : https://stats.stackexchange.com/questions/471267/plotting-gams-on-response-scale-with-multiple-smooth-and-linear-terms
# Wilk et al., 2018: https://doi.org/10.1029/2018JG004506 (plotting response values)
