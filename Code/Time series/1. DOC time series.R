###############################################################################
####### BPWTP Data -- data analysis ##########
####### Danielle Spence ##########
####### Created 5/4/2023 ########
###############################################################################
### Clear memory

rm(list = ls())

library(pacman)
p_load(tidyverse, ggplot2, mgcv, gratia, readr, GGally, dplyr, lubridate, plyr, itsadug,
       cowplot, tibble,viridis, install = TRUE)

setwd("C:/Users/danis/OneDrive/R/DW Indicators")

papertheme <- theme_bw(base_size=14, base_family = 'Arial') 

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
    )


#Remove rows without DOC data
df <- df[complete.cases(df$DOC),]
# 1655 obs


##----------------------------------------------------------------------------##
## 2. Check trends, distributions
##----------------------------------------------------------------------------##

df %>%
  ggplot(aes(x = year, y = DOC)) +
  geom_point()+
  geom_line()+
  ggtitle("Taste and Odour Number")

df %>%
  select(DOC) %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram() # DOC has a strong positive skew, suggest using gamma or tweedie dist


##----------------------------------------------------------------------------##
## 3. Modelling DOC over the Years
##----------------------------------------------------------------------------##

m<- gam(DOC ~ s(DOY, bs="cc", k=12)+
           s(year, bs = "cr", k=12), 
           knots=list(DOY=c(0, 366.5)),
           data = df, method = "REML", family = tw(link="log"))

summary(m) # deviance explained = 67%; REML  = 1426
draw(m)

k.check(m) 

##----------------------------------------------------------------------------##
## 4. Checking & controlling for autocorrelation
##----------------------------------------------------------------------------##

layout(matrix(1:2, ncol = 2))
acf(resid(m), lag.max = 36, main = "ACF") 
pacf(resid(m), lag.max = 36, main = "pACF")
layout(1)

## Controlling for autocorrelation

ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")

m<- gamm(DOC ~ s(DOY, bs="cc")+
           s(year, bs = "cr"), 
         knots=list(DOY=c(0, 366.5)), data = df, method = "REML", family = Gamma(link="log"))

## AR(1)
m1 <- gamm(DOC ~ s(DOY, bs="cc")+
             s(year, bs = "cr"), 
           knots=list(DOY=c(0, 366.5)),
            data = df, correlation = corARMA(form = ~ 1|year, p = 1),
            control = ctrl, method = "REML", family = Gamma(link="log"))
## AR(2)
m2 <- gamm(DOC ~ s(DOY, bs="cc")+
             s(year, bs = "cr"), 
           knots=list(DOY=c(0, 366.5)),
           data = df, correlation = corARMA(form = ~ 1|year, p = 2),
           control = ctrl, method = "REML", family = Gamma(link="log"))

## AR(3)
m3 <- gamm(DOC ~ s(DOY, bs="cc")+
             s(year, bs = "cr"), 
           knots=list(DOY=c(0, 366.5)),
            data = df, correlation = corARMA(form = ~ 1|year, p = 3),
            control = ctrl, method = "REML", family = Gamma(link="log"))

anova(m$lme, m1$lme, m2$lme, m3$lme) ## shows that m1 is best


layout(matrix(1:2, ncol = 2))
res <- resid(m2$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(2) errors")
acf(res, lag.max = 36, main = "pACF- AR(2) errors")
layout(1) # this is better

summary(m2$gam) 

k.check(m2$gam) 

appraise(m2$gam, method = "simulate") # patterns in residuals vs. linear and bs. observed...need to decide whether log is the corect transformation? Similar to count data...

p<-draw(m2$gam)
p
#ggsave('output/DOC year.png', p, height = 8, width  = 10)



##----------------------------------------------------------------------------##
## 5. PLot predictions
##----------------------------------------------------------------------------##


new_data_DOY <- with(df, expand.grid(DOY = seq(min(DOY, na.rm = TRUE), max(DOY, na.rm = TRUE), length = 200),
                                      year = median(year)))


DOY.pred <- predict(m2$gam, newdata = new_data_DOY, type = "terms", se.fit = TRUE)



whichCols <- grep("DOY", colnames(DOY.pred$fit))
whichColsSE <- grep("DOY", colnames(DOY.pred$se.fit))
new_data_DOY <- cbind(new_data_DOY, Fitted = DOY.pred$fit[, whichCols], 
                      se.Fitted = DOY.pred$se.fit[, whichColsSE])


limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

## make into original limits
new_data_DOY <- with(new_data_DOY, transform(new_data_DOY, Fittedplus = Fitted + se.Fitted))
new_data_DOY <- with(new_data_DOY, transform(new_data_DOY, Fittedminus = Fitted - se.Fitted))

shiftDOY <- attr(predict(m2$gam, newdata = new_data_DOY, type = "iterms"), "constant")
DOY.pdatnorm <- new_data_DOY
DOY.pdatnorm <- with(DOY.pdatnorm, transform(DOY.pdatnorm, Fitted = Fitted + shiftDOY, 
                                             Fittedplus = Fittedplus + shiftDOY, 
                                             Fittedminus = Fittedminus + shiftDOY))

# backtransform fitted data...
DOY.pdatnorm$Fitted<-exp(DOY.pdatnorm$Fitted)
DOY.pdatnorm$Fittedplus<-exp(DOY.pdatnorm$Fittedplus)
DOY.pdatnorm$Fittedminus<-exp(DOY.pdatnorm$Fittedminus)




DOYplot <- ggplot(DOY.pdatnorm, aes(x = DOY, y = Fitted)) +
  papertheme+
  geom_line() +
  geom_point(data = df, aes(DOY,DOC), alpha = 0.1) +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25, fill = '#165459B2') +  
  scale_x_continuous(expand = c(0, 0), breaks = c(5, 60,  121, 182, 244, 305, 350), labels = c("Jan", "Mar", "May", "Jul", "Sept", "Nov", "Dec")) +
  scale_y_continuous(expand = c(0, 0)) +  # Remove padding on x-axis
  xlab(expression(paste("DOY "))) + ylab(expression(paste("DOC (mg ", L^-1,")")))

DOYplot



saveRDS(DOYplot, "DOY_DOC_AC.rds")


## intra-annual range in concentraiton?

dx<- DOY.pdatnorm %>%
  summarize(
    mean_conc = mean(Fitted, na.rm = TRUE),
    range_percent = ((max(Fitted, na.rm = TRUE) - min(Fitted, na.rm = TRUE)) / mean_conc) * 100
  )
dx

mean(dx$range_percent)



##plot year


new_data_year <- with(df, expand.grid(year = seq(min(year, na.rm = TRUE), max(year, na.rm = TRUE), length = 200),
                                     DOY = median(DOY)))


year.pred <- predict(m2$gam, newdata = new_data_year, type = "terms", se.fit = TRUE)



whichCols <- grep("year", colnames(year.pred$fit))
whichColsSE <- grep("year", colnames(year.pred$se.fit))
new_data_year <- cbind(new_data_year, Fitted = year.pred$fit[, whichCols], 
                      se.Fitted = year.pred$se.fit[, whichColsSE])


limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

## make into original limits
new_data_year <- with(new_data_year, transform(new_data_year, Fittedplus = Fitted + se.Fitted))
new_data_year <- with(new_data_year, transform(new_data_year, Fittedminus = Fitted - se.Fitted))

shiftyear <- attr(predict(m2$gam, newdata = new_data_year, type = "iterms"), "constant")
year.pdatnorm <- new_data_year
year.pdatnorm <- with(year.pdatnorm, transform(year.pdatnorm, Fitted = Fitted + shiftyear, 
                                             Fittedplus = Fittedplus + shiftyear, 
                                             Fittedminus = Fittedminus + shiftyear))

# backtransform fitted data...
year.pdatnorm$Fitted<-exp(year.pdatnorm$Fitted)
year.pdatnorm$Fittedplus<-exp(year.pdatnorm$Fittedplus)
year.pdatnorm$Fittedminus<-exp(year.pdatnorm$Fittedminus)




yearplot <- ggplot(year.pdatnorm, aes(x = year, y = Fitted)) +
  papertheme+
  geom_line() +
  geom_point(data = df, aes(year,DOC), alpha = 0.1) +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25, fill = '#165459B2') +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +  # Remove padding on x-axis
 xlab(expression(paste("Year"))) + ylab(expression(paste("DOC (mg ", L^-1,")")))

yearplot

# inter-annnual range in DOC

dy<- year.pdatnorm %>%
  summarize(
    mean_conc = mean(Fitted, na.rm = TRUE),
    range_percent = ((max(Fitted, na.rm = TRUE) - min(Fitted, na.rm = TRUE)) / mean_conc) * 100
  )
dy

mean(dy$range_percent)


saveRDS(yearplot, "year_DOC_AC.rds")

p_all<- DOYplot + yearplot
p_all

saveRDS(p_all, "time_DOC.rds")


##--------------------------------------------------##
## References
# Simpson, 2014: https://fromthebottomoftheheap.net/2014/05/09/modelling-seasonal-data-with-gam/
# Painter et al., 2023: https://doi.org/10.1002/ecs2.4472
# selection of tensor prods: https://pub.towardsai.net/gams-and-smoothing-splines-part-2-tensor-product-splines-97928f226a2c
# plotting predicted vs. observed: Baron, 2023: https://harvest.usask.ca/handle/10388/14623