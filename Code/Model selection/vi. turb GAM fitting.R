 ###############################################################################
####### BPWTP Data -- GAMs ##########
####### Danielle Spence ##########
####### Created 5/4/2023 ########
###############################################################################
### Clear memory
rm(list = ls())

library(pacman)
p_load(tidyverse, ggplot2, mgcv, gratia, readr, GGally, dplyr,  mgcViz,  lubridate, plyr,
        cowplot, tibble,viridis,  install = TRUE)

setwd("C:/Users/danis/OneDrive/R/DW Indicators") 

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
df$TN <- rowSums(df[,c("NO3", "NH3", "orgN")], na.rm=TRUE)

## add in PDO and SOI 

clim <- read_csv("data/SOI_PDO data no lag.csv")

df1<- merge(df, clim, by="date")

##----------------------------------------------------------------------------##
## 3. Check trends, distributions, correlations
##----------------------------------------------------------------------------##

# basic plot

df %>%
  ggplot(aes(x = year, y = turb)) +
  geom_line()+
  ggtitle("turb")

df %>%
  dplyr::select(turb, SRP,  DIN) %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()


## transform positively skewed data (turb, SRP, DIN)

df1$turb_sqrt<- sqrt(df1$turb)
df1$srp_sqrt<- sqrt(df1$SRP)
df1$din_sqrt<- sqrt(df1$DIN)

df1 %>%
  dplyr::select(turb_sqrt, srp_sqrt,  W_temp,  din_sqrt) %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()


## add ph??? re: Zepernick et al. 2024
dfx<-df1%>%
  select(turb, SRP, DIN, QLD, W_temp, PDO_3MON_AVG, SOI_3MON_AVG)
## visualize correlations
p<- ggpairs(dfx[])  # can see that turb and turbidity, temperature, Org_N and TP are correlated. Ammonia and conductivity are not. (NH3, NO3, and TP correlated; TP and temp). Org N and TP are similarly correlated with turb (0.640 and 0.667, respectively)
p 


ggsave('output/Correlations.png', p, height = 6, width  = 8)


##----------------------------------------------------------------------------##
## 1. Check for best fit family
##----------------------------------------------------------------------------##


ma <- gam(turb ~ s(SRP) +  
            s(DIN)+
            s(W_temp)+
            s(QLD)+
            te(PDO_3MON_AVG,SOI_3MON_AVG)+
            te(year, DOY, bs = c("cr", "cc"), k=c(10,12)),
          select = TRUE,
          data = df1, method = "ML", family = scat())

mb <- gam(turb ~ s(SRP) +  
            s(DIN)+
            s(W_temp)+
            s(QLD)+
            te(PDO_3MON_AVG,SOI_3MON_AVG)+
            te(year, DOY, bs = c("cr", "cc"), k=c(10,12)),
          select = TRUE,
          data = df1, method = "ML", family = gaussian)

mc <- gam(turb ~ s(SRP) +  
            s(DIN)+
            s(W_temp)+
            s(QLD)+
            te(PDO_3MON_AVG,SOI_3MON_AVG)+
            te(year, DOY, bs = c("cr", "cc"), k=c(10,12)),
          select = TRUE,
          data = df1, method = "ML", family = tw())


md <- gam(turb ~ s(SRP) +  
            s(DIN)+
            s(W_temp)+
            s(QLD)+
            te(PDO_3MON_AVG,SOI_3MON_AVG)+
            te(year, DOY, bs = c("cr", "cc"), k=c(10,12)),
          knots=list(DOY=c(0, 366.5)),
          select = TRUE,
          data = df1, method = "ML", family = Gamma(link = "log"))

me <- gam(turb ~ s(SRP) +  
            s(DIN)+
            s(W_temp)+
            s(QLD)+
            te(PDO_3MON_AVG,SOI_3MON_AVG)+
            te(year, DOY, bs = c("cr", "cc"), k=c(10,12)),
          select = TRUE,
          data = df1, method = "ML",  family = inverse.gaussian(link = "log"))

# Compare AIC values
AIC(ma, mb, mc, md, me) # tweedie or gamma best fit
summary(md)
draw(md)

##--- check which variables should stay/be removed

nullvars <- c('turb', "SRP", "DIN", "QLD", "W_temp", "PDO_3MON_AVG", "SOI_3MON_AVG", "year", "DOY")
mumdat <- df1[complete.cases(df1[,nullvars]),] # want to AIC compare, and NA diff
#   between vars
mummod <- update(mc, data=mumdat, na.action="na.fail")
mummoddin <- update(mummod, .~.-s(DIN))
mummodsrp <- update(mummod, .~.-s(SRP))
mummodwtemp <- update(mummod, .~.-s(W_temp))
mummodQLD <- update(mummod, .~.-s(QLD))
mummodTE <- update(mummod, .~.-te(PDO_3MON_AVG,SOI_3MON_AVG))

AIC(mummod, mummoddin, mummodsrp, mummodwtemp, mummodQLD, mummodTE) # removing predictors doesnt improve model fit


##----------------------------------------------------------------------------##
## 5. check whether SOI and PDO should be included as individual or tensor
##----------------------------------------------------------------------------##

mi <- gam(turb ~ s(SRP) +  
            s(DIN)+
            #s(W_temp)+
            s(QLD)+
            s(PDO_3MON_AVG)+
            s(SOI_3MON_AVG)+
            te(year, DOY, bs = c("cr", "cc"), k=c(10,12)),
          #s(fYear, bs="re"),
          select = TRUE,
          data = df1, method = "REML", family = tw())

summary(mi)
gam.check(mi)
draw(mi, all=TRUE)


# year as tensor
mI <- gam(turb ~ s(SRP) +  
            s(DIN)+
            # s(W_temp)+
            s(QLD)+
            ti(PDO_3MON_AVG) +
            ti(SOI_3MON_AVG)+
            ti(PDO_3MON_AVG,SOI_3MON_AVG)+
            te(year, DOY, bs = c("cr", "cc"), k=c(10,12)),
          select = TRUE,
          data = df1, method = "REML", family = tw())

anova(mi, mI, test = "LRT") #tensor interaction is better.
draw(mI)



##----------------------------------------------------------------------------##
## 6. Assess model fit & autocorrelation
##----------------------------------------------------------------------------##

## i. Appraise model fit
r<- appraise(m_t, point_col = 'steelblue', point_alpha = 0.5, n_bins = 'fd') & 
  theme(plot.tag = element_text(face = 'bold')) & theme_bw()
r

ggsave('output/turb GAM RESIDUALS TW SQ RT.png', r, height = 6, width  = 8)

# ii. Check for autocorrelation

layout(matrix(1:2, ncol = 2))
acf(resid(m_t), lag.max = 36, main = "ACF") 
pacf(resid(m_t), lag.max = 36, main = "pACF")
layout(1)


##----------------------------------------------------------------------------##
## 6. Control for autocorrelation
##----------------------------------------------------------------------------##


ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")

m<- gamm(turb ~ s(SRP, k=5) +  
           s(DIN, k=6)+
           s(QLD, k=5)+ 
           te(PDO_3MON_AVG,SOI_3MON_AVG)+
           te(year, DOY, bs = c("cr", "cc"), k=c(10,12)),
        # select = TRUE,
         data = df1, method = "REML", 
         family = Gamma(link = "log"))

## AR(1)
m1 <- gamm(turb ~ s(SRP, k=10) +  
             s(DIN, k=10)+
             s(QLD, k=5)+ 
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc"), k=c(10,12)),
          # select = TRUE,
           data = df1, family = Gamma(link = "log"),
           correlation = corARMA(form = ~ 1|year, p = 1),
           control = ctrl, method = "REML")
## AR(2)
m2 <- gamm(turb ~ s(SRP, k=5) +  
             s(DIN, k=6)+
             s(QLD, k=5)+ 
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc"), k=c(10,12)),
           #select = TRUE,
           family = Gamma(link = "log"),
           data = df1, correlation = corARMA(form = ~ 1|year, p = 2),
           control = ctrl, method = "REML")

## AR(3)
m3 <- gamm(turb ~ s(SRP, k=5) +  
             s(DIN, k=6)+
             s(QLD, k=5)+ 
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc"), k=c(10,12)),
           
           data = df1, family = Gamma(link = "log"),
           correlation = corARMA(form = ~ 1|year, p = 3),
           control = ctrl, method = "REML")

anova(m$lme, m1$lme, m2$lme, m3$lme) ## shows that m1 is best


layout(matrix(1:2, ncol = 2))
res <- resid(m1$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
acf(res, lag.max = 36, main = "pACF- AR(1) errors")
layout(1) # this is better

summary(m1$gam) 

k.check(m1$gam) 

## i. Appraise model fit
r1<- appraise(m1$gam, point_col = 'steelblue', point_alpha = 0.5, n_bins = 'fd') & 
  theme(plot.tag = element_text(face = 'bold')) & theme_bw()
r1

ggsave('output/turb GAM_preds RESIDUALS AUTOCORR GAMMA.png', r1, height = 6, width  = 8)

p1<- draw(m1$gam, residuals = TRUE)& theme_bw()
p1

#ggsave('output/turb GAM_ autocorrelation no interaction SOI.png', p1, height = 10, width  = 10)

##----------------------------------------------------------------------------##
## 5. Check which fit is best for time
##----------------------------------------------------------------------------##


# year not added

m <- gam(turb ~ s(SRP) +  
            s(DIN)+ 
            #s(W_temp)+
            s(QLD)+
            te(PDO_3MON_AVG,SOI_3MON_AVG),
          #te(year, DOY, bs = c("cr", "cc")),
          select = TRUE,
          data = df1, method = "REML", Gamma(link = "log"))

m1 <- gam(turb ~ s(SRP) +  
            s(DIN)+ 
            #s(W_temp)+
            s(QLD)+
            te(PDO_3MON_AVG,SOI_3MON_AVG)+
            s(year,  bs = "cr")+
            s(DOY, bs = "cc"),
          knots=list(DOY=c(0, 366.5)),
          select = TRUE,
          data = df1, method = "REML", family = Gamma(link = "log"))

m2 <- gam(turb ~ s(SRP) +  
            s(DIN)+ 
            s(W_temp)+
            s(QLD)+
            te(PDO_3MON_AVG,SOI_3MON_AVG)+
            te(year, DOY, bs = c("cr", "cc")),
          knots=list(DOY=c(0, 366.5)),
          select = TRUE,
          data = df1, method = "REML", family = Gamma(link = "log"))

AIC(m, m1, m2) #m2 

summary(m2) # Deviance explained =71%, REML = 3504, r2 = 0.463

k.check(m1)# k-index looks good 

p1<- draw(m1,  residuals = TRUE)& theme_bw() 
p1


# model 2 is a better fit in terms of AIC and BIC.

##--------------------------------------------------##
## References
# Simpson, 2014: https://fromthebottomoftheheap.net/2014/05/09/modelling-seasonal-data-with-gam/
# Hayes et al, 2020: https://doi.org/10.1002/lol2.10164
# Simpson 2020 : https://stats.stackexchange.com/questions/471267/plotting-gams-on-response-scale-with-multiple-smooth-and-linear-terms
# Wilk et al., 2018: https://doi.org/10.1029/2018JG004506 (plotting response values)
