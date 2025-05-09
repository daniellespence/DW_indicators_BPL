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

setwd("C:/Users/danis/OneDrive/R/DW indicators")

papertheme <- theme_bw(base_size=14, base_family = 'Arial') 

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
df$TN <- rowSums(df[,c("NO3", "NH3", "orgN")], na.rm=TRUE)

## add in PDO and SOI 

clim <- read_csv("data/SOI_PDO data no lag.csv")

df1<- merge(df, clim, by="date")


##----------------------------------------------------------------------------##
## 3. Check trends, distributions, correlations
##----------------------------------------------------------------------------##

# basic plot

df1 %>%
  ggplot(aes(x = year, y = TON)) +
  geom_line()+
  ggtitle("TON")

df1 %>%
  dplyr::select(TON, SRP,  W_temp,  DIN, pH) %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()

# assess correlations

dfx<-df1%>%
  select(TON, SRP, TP, DIN, QLD,W_temp, PDO_3MON_AVG, SOI_3MON_AVG)
## visualize correlations
p<- ggpairs(dfx[])  # can see that TON and turbidity, temperature, Org_N and TP are correlated. Ammonia and conductivity are not. (NH3, NO3, and TP correlated; TP and temp). Org N and TP are similarly correlated with TON (0.640 and 0.667, respectively)
p 


ggsave('output/TON Correlations.png', p, height = 8, width  = 10)


##----------------------------------------------------------------------------##
## 5. Run the full model
##----------------------------------------------------------------------------##

m <- gam(Odour ~ s(SRP) +  
            s(DIN)+ 
            s(W_temp)+
            s(QLD)+
            te(PDO_3MON_AVG,SOI_3MON_AVG)+
            te(year, DOY, bs = c("cr", "cc")),
            select = TRUE,
            data = df1, method = "REML")


m1 <- gam(Odour ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
           select = TRUE,
           data = df1, method = "REML", family = tw())

m2 <- gam(Odour ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
          knots=list(DOY=c(0, 366.5)),
           select = TRUE,
           data = df1, method = "REML", family = Gamma(link = "log"))

m3 <- gam(Odour ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
           select = TRUE,
           data = df1, method = "REML", family = scat())

m4 <- gam(Odour ~ s(SRP) +  
            s(DIN)+ 
            s(W_temp)+
            s(QLD)+
            te(PDO_3MON_AVG,SOI_3MON_AVG)+
            te(year, DOY, bs = c("cr", "cc")),
          select = TRUE,
          data = df1, method = "REML", family = nb())

m5 <- gam(Odour ~ s(SRP) +  
            s(DIN)+ 
            s(W_temp)+
            s(QLD)+
            te(PDO_3MON_AVG,SOI_3MON_AVG)+
            te(year, DOY, bs = c("cr", "cc")),
          select = TRUE,
          data = df1, method = "REML", family = poisson)

AIC(m, m1, m2, m3, m4, m5) #m2 (gamma family is best fit)

summary(m2) # Deviance explained =71%, REML = 3504, r2 = 0.463

k.check(m2)# k-index looks good 

p1<- draw(m2,  residuals = TRUE)& theme_bw() 
p1

#ggsave('output/TON GAM_SQRT SOI PDO interact.png', p1, height = 10, width  = 10)


##--- check which variables should stay/be removed

nullvars <- c('Odour', "SRP", "DIN", "QLD", "W_temp", "PDO_3MON_AVG", "SOI_3MON_AVG", "year", "DOY")
mumdat <- df1[complete.cases(df1[,nullvars]),] # want to AIC compare, and NA diff
#   between vars
mummod <- update(m2, data=mumdat, na.action="na.fail")
mummoddin <- update(mummod, .~.-s(DIN))
mummodsrp <- update(mummod, .~.-s(SRP))
mummodwtemp <- update(mummod, .~.-s(W_temp))
mummodQLD <- update(mummod, .~.-s(QLD))
mummodTE <- update(mummod, .~.-te(PDO_3MON_AVG,SOI_3MON_AVG))

AIC(mummod, mummoddin,mummodsrp, mummodwtemp, mummodQLD, mummodTE) # 

summary(mummodwtemp)

# AIC does not improve when removing additional predictors.


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
acf(resid(m2), lag.max = 36, main = "ACF") 
pacf(resid(m2), lag.max = 36, main = "pACF")
layout(1)

##----------------------------------------------------------------------------##
## 6. Control for autocorrelation
##----------------------------------------------------------------------------##


ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")

m<- gamm(Odour ~ s(SRP) +  
           s(DIN)+ 
           s(W_temp)+
           s(QLD)+
           te(PDO_3MON_AVG,SOI_3MON_AVG)+
           te(year, DOY, bs = c("cr", "cc")),
         
         data = df1,method = "REML", family = Gamma(link = "log"))

## AR(1)
m1 <- gamm(Odour ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
           #select = TRUE,
           data = df1, correlation = corARMA(form = ~ 1|year, p = 1),
           control = ctrl, method = "REML", family = Gamma(link = "log"))
## AR(2)
m2 <- gamm(Odour ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
           #select = TRUE,
           data = df1, correlation = corARMA(form = ~ 1|year, p = 2),
           control = ctrl, method = "REML", family = Gamma(link = "log"))

## AR(3)
m3 <- gamm(Odour ~ s(SRP) +  
             s(DIN)+ 
             s(W_temp)+
             s(QLD)+
             te(PDO_3MON_AVG,SOI_3MON_AVG)+
             te(year, DOY, bs = c("cr", "cc")),
           #select = TRUE,
           data = df1, correlation = corARMA(form = ~ 1|year, p = 3),
           control = ctrl, method = "REML", family = Gamma(link = "log"))

AIC(m$gam, m1$gam, m2$gam, m3$gam) ## shows that m1 is best


layout(matrix(1:2, ncol = 2))
res <- resid(m1$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(2) errors")
acf(res, lag.max = 36, main = "pACF- AR(2) errors")
layout(1) # this is better


#summary(m2) 

k.check(m1$gam) 

## i. Appraise model fit
appraise(m1$gam, point_col = 'steelblue', point_alpha = 0.5, n_bins = 'fd') & 
  theme(plot.tag = element_text(face = 'bold')) & theme_bw()

p1<- draw(m1$gam, residuals = TRUE,  select = smooths(m1$gam)[1:5])& theme_bw()
p1

#ggsave('output/Odour GAM_ autocorrelation no interaction SOI.png', p1, height = 10, width  = 10)



##----------------------------------------------------------------------------##
## 5. Compare time
##----------------------------------------------------------------------------##

m <- gam(Odour ~ s(SRP) +  
           s(DIN)+ 
           #s(W_temp)+
           s(QLD)+
           te(PDO_3MON_AVG,SOI_3MON_AVG),
           #te(year, DOY, bs = c("cr", "cc")),
         select = TRUE,
         data = df1, method = "REML", family = Gamma(link = "log"))


m1 <- gam(Odour ~ s(SRP) +  
            s(DIN)+ 
            #s(W_temp)+
            s(QLD)+
            te(PDO_3MON_AVG,SOI_3MON_AVG)+
            s(year,  bs = "cr")+
            s(DOY, bs = "cc"),
          knots=list(DOY=c(0, 366.5)),
          select = TRUE,
          data = df1, method = "REML", family = Gamma(link = "log"))

m2 <- gam(Odour ~ s(SRP) +  
            s(DIN)+ 
            s(W_temp)+
            s(QLD)+
            te(PDO_3MON_AVG,SOI_3MON_AVG)+
            te(year, DOY, bs = c("cr", "cc")),
          knots=list(DOY=c(0, 366.5)),
          select = TRUE,
          data = df1, method = "REML", family = Gamma(link = "log"))

AIC(m, m1, m2) #m1 

summary(m1) # Deviance explained =71%, REML = 3504, r2 = 0.463

k.check(m1)# k-index looks good 

p1<- draw(m1,  residuals = TRUE)& theme_bw() 
p1





##--------------------------------------------------##
## References
# Simpson, 2014: https://fromthebottomoftheheap.net/2014/05/09/modelling-seasonal-data-with-gam/
# Hayes et al, 2020: https://doi.org/10.1002/lol2.10164
# Simpson 2020 : https://stats.stackexchange.com/questions/471267/plotting-gams-on-response-scale-with-multiple-smooth-and-linear-terms
# Wilk et al., 2018: https://doi.org/10.1029/2018JG004506 (plotting response values)
