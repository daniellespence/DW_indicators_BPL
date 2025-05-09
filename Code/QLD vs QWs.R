###############################################################################
####### BPWTP Data -- GAMs ##########
####### Danielle Spence ##########
####### Created 8/26/2024 ########
###############################################################################
### Clear memory
rm(list = ls())

library(pacman)
p_load(tidyverse, ggplot2, readr, GGally, dplyr, patchwork, install = TRUE)

library(dplyr)

setwd("C:/Users/danis/OneDrive/R/DW indicators")


##----------------------------------------------------------------------------##
## 1. Read in data 
##----------------------------------------------------------------------------##


df <- read_csv("Data/daily flow.csv")

# add in Year and nMonth for numeric month and a proper Date class

df <- df%>% 
  filter(year %in% c(1990:2022))



# rename for ease of use
df <- df %>%
  dplyr::rename(
    QBP = "combined_05JG004.cms",
    QLD = "SK05JG006.cms",
    QWS = "RC_IC_cms",
    )


df<- df%>%
  select(date, year, DOY, QWS, QLD)

#Remove rows without TON data
df <- df[complete.cases(df$QWS),]
df <- df[complete.cases(df$QLD),]

# Fit a linear regression model
m <- lm(QLD ~ QWS, data = df)

# Display the summary of the model
summary(m)

# Add the regression line to a scatter plot
p<-ggplot(df, aes(x=QWS, y=QLD, colour=year))+
  geom_point(size=2)+
  scale_colour_viridis_c(option = "mako", direction = 1, name = "Year") +
  xlab(expression(paste("Flow from watershed sources (",m^3,  s^-1,")")))+ ylab(expression(paste("Flow from Lake Diefenbaker ( ",m^3,  s^-1,")")))+
  theme_bw(base_size = 14)
p
ggsave('output/scatterplot QWS QLD 1990.png', p, height = 6, width  = 8)


## 
dfx<-df%>%
  select(QWS, QLD)
## visualize correlations
p<- ggpairs(dfx[]) + theme_bw()  # can see that DOC and turbidity, temperature, Org_N and TP are correlated. Ammonia and conductivity are not. (NH3, NO3, and TP correlated; TP and temp). Org N and TP are similarly correlated with DOC (0.640 and 0.667, respectively)
p 


#ggsave('output/QWS QLD Correlations.png', p, height = 8, width  = 10)

summary(df$QWS)

# Convert yQWS# Convert y to a factor to indicate zero vs non-zero
df$QWS_group <- factor(ifelse(df$QWS >= 0.2069, "Above/Equal", "Below"))


p <-ggplot(df, aes(x = QLD, fill=QWS_group)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Above/Equal" = "#165460A1", "Below" = "#165430A1"), name="Average QWS")+
  ylab(expression(paste("Conditional density")))+ xlab(expression(paste("Flow from Lake Diefenbaker ( ",m^3, "/", s,")")))+
  theme_bw(base_size = 14)
p

#ggsave('output/conditional density QWS QLD.png', p, height = 8, width  = 10)



##-----------------------------------Flow ratio vs. year------------------------------- ####
df$QWS1<- df$QWS+0.00000000000000000000000001
df$ratioLD <- (df$QLD / df$QWS1)

df$ratioWS <- (df$QWS / df$QLD)

# annual_summary <- df %>%
#   group_by(year) %>%
#   summarize(Average_ratio = mean(ratio, na.rm = TRUE), 
#             Max_ratio = max(ratio, na.rm = TRUE),
#             .groups = 'drop')


max <- ggplot(df, aes(x = year, y = ratioLD))+
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = 'Year', y = expression(paste('QLD:QWS')))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(df$ratioLD) * 1.01))

max


#ggsave('output/QLD_QWS year.png', max, height = 8, width  = 10)


## showing seasonal flows


# first create season

seasons = function(x){
  if(x %in% 3:5) return("Spring")
  if(x %in% 6:8) return("Summer")
  if(x %in% 9:11) return("Fall")
  if(x %in% c(12,1,2)) return("Winter")
  
}


df$Season = sapply(month(df$date), seasons)





max <- ggplot(df, aes(x = year, y = ratioLD))+
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = 'Year', y = expression(paste("QLD:QWS")))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(df$ratioLD) * 1.01))

max


p1<- max + facet_grid(~factor(Season, levels=c('Fall', 'Winter', 'Spring', 'Summer')))
p1


#ggsave('output/QLD_QWS seaspnal.png', p1, height = 8, width  = 10)

##-----------------------------------Flows across years------------------------------- ####

long_flow_data <- df %>%
  pivot_longer(
    cols = c(QWS, QLD),
    names_to = "Source",
    values_to = "Total_Volume"
  )
annual_summary <- long_flow_data %>%
  group_by(year, Source, Season) %>%
  summarize(Average_Volume = mean(Total_Volume, na.rm = TRUE), 
            Max_Volume = max(Total_Volume, na.rm = TRUE),
            .groups = 'drop')


max <- ggplot(annual_summary, aes(x = year, y = Max_Volume, fill = Source))+
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(
    values = c("QWS" = "#165459B2", "QLD" = "gray60")
  ) +
  #facet_wrap(~ Season, ncol = 2) +
  labs(x = 'Year', y = expression(paste("Maximum flow rate (",m^3, s^-1,")")))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(annual_summary$Max_Volume) * 1.01))

max

p1<- max + facet_grid(~factor(Season, levels=c('Fall', 'Winter', 'Spring', 'Summer')))
p1


mean <- ggplot(annual_summary, aes(x = year, y = Average_Volume, fill = Source))+
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(
    values = c("QWS" = "#165459B2", "QLD" = "gray60")
  ) +
  #facet_wrap(~ Season, ncol = 2) +
  labs(x = 'Year', y = expression(paste("Average flow rate (",m^3, s^-1,")")))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(annual_summary$Average_Volume) * 1.01))

mean

p2<- mean + facet_grid(~factor(Season, levels=c('Fall', 'Winter', 'Spring', 'Summer')))
p2


combined_plot <- p2 / p1 +
  plot_layout(guides = "collect") &   theme(legend.position = "bottom")
combined_plot

ggsave('output/annual max and mean seasonal QWS QLD_leg.png', combined_plot, height = 8, width  = 8)




#---- summer averages?
annual_summer <- long_flow_data %>%
  filter(Season == 'Summer',
         Source == 'QLD')

annual_summer<- annual_summer %>%
  group_by(year) %>%
  summarize(Average_Volume = mean(Total_Volume, na.rm = TRUE), 
            Max_Volume = max(Total_Volume, na.rm = TRUE),
            .groups = 'drop')


max <- ggplot(annual_summer, aes(x = year, y = Average_Volume))+
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = 'Year', y = expression(paste("Average Summer QLD")))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(annual_summer$Average_Volume) * 1.01))

max
ggsave('output/annual mean summer QLD.png', max, height = 6, width  = 8)


# QLD
# average summer = 5.5
# average max = 7.9
# 2014 & 2015 average = 1.0 and 4.3
# 2014 & 2015 max = 2.6 and 8.5

mean(annual_summer$Average_Volume)

mean# select only summer months
dfy<- df%>%
  filter(month(date) %in% c(6, 7, 8))

dfy <- dfy %>%
  mutate(Month = month(date, label = TRUE, abbr = TRUE))  # Convert to month factor


p2<- ggplot(dfy, aes(x = Month, y = QLD)) +
  geom_boxplot(fill = "#165459B2", color = "black") +
  labs(
    x = "Year",
    y = expression(paste("Average flow rate ( ",m^3,s^-1,")")),
    fill = "Year"
  ) +
  theme_bw()+
  #theme(legend.position = "none")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(dfy$QLD) * 1.01))
p2

ggsave('output/summer QLD boxplot.png', p2, height = 6, width  = 8)


