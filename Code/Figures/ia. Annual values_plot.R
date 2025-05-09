###############################################################################
####### CH.1 BPWTP -- data processing ##########
####### Danielle Spence ##########
####### Created 3/30/2021 ########
###############################################################################
### Clear memory
rm(list = ls())

library(pacman)
p_load(tidyverse, dplyr, tidyr, geomtextpath, ggsci, cowplot, patchwork,
       install = TRUE)


setwd("C:/Users/danis/OneDrive/R/DW Indicators")

# Example data with warm (red) and cold (blue) phases

#SOI
phases <- data.frame(
  phase = c("Warm", "Cold", "Warm", "Cold", "Warm", "Cold", "Warm", "Cold", "Warm", "Cold"),
  start = c(1990, 1995, 1997, 1998, 2002, 2007, 2014, 2016, 2018, 2020),
  end = c(1995, 1997, 1998, 2002, 2007, 2014, 2016, 2018, 2020, 2022)
)


#PDO
phases_PDO <- data.frame(
  phase = c("Cold","Warm","Cold", "Warm", "Cold", "Warm", "Cold", "Warm", "Cold"),
  start = c(1990, 1992, 1994, 1995, 1998, 2003, 2004, 2014, 2016),
  end = c(1992, 1994, 1995, 1998, 2003, 2004, 2014, 2016, 2022)
)

#In-phase

in_phase <- data.frame(
  phase = c("Warm","Warm","Cold", "Warm", "Cold", "Warm", "Cold", "Cold"),
  start = c(1992, 1997, 1998, 2003, 2007, 2014, 2016, 2020 ),
  end = c(1993, 1998, 2001, 2004, 2014, 2016, 2019, 2022 )
)

##----------------------------------------------------------------------------##
## 1. conductivity
##----------------------------------------------------------------------------##
## Read in conductivity data 
df_c <- read_csv("data/bpgamdataCLEAN_Conductivity.csv")


# add in Year and nMonth for numeric month and a proper Date class


df_c$Date<- as.Date(df_c$Date, "%m/%d/%Y")

df_c <- mutate(df_c,
             year = as.numeric(format(Date,'%Y')),
             DOY = as.numeric(format(Date,'%j')),
             nMonth = as.numeric(format(Date,'%m')))%>% 
  filter(year %in% c(1990:2022))


## order the data
df_c <- df_c[with(df_c, order(Date)),]


# rename for ease of use
df_c <- df_c%>%
  dplyr::rename(
    date = "Date" )

#Remove rows without Conductivity data
df_c <- df_c[complete.cases(df_c$Conductivity),]


ann <-aggregate(df_c$Conductivity,  
                by=list(df_c$year),  
                FUN=mean,
                na.rm=TRUE) 
ann<- ann %>% dplyr::rename(mean = x)


ann1 <-aggregate(df_c$Conductivity,  
                by=list(df_c$year),  
                FUN = max,
                na.rm=TRUE)
ann1 <- ann1 %>% dplyr::rename(max = x)


ann2 <-aggregate(df_c$Conductivity,  
                by=list(df_c$year),  
                FUN = min,
                na.rm=TRUE)
ann2 <- ann2 %>% dplyr::rename(min = x)

anncl<- left_join(ann, ann1,  by= "Group.1")
anncl<- left_join(anncl, ann2, by= "Group.1")

anncl <- mutate(anncl,
               year = as.numeric(Group.1))
               
p_anCL<- anncl %>%
  ggplot(aes(x = year, y = mean)) +
  geom_point()+
  geom_ribbon(aes(ymin = min, ymax = max), fill = '#165459B2', alpha = 1/2) +
  geom_line()+
  theme_classic() +
  labs(fill = "SOI Phase")+
  ggtitle(expression(paste("Conductivity")))+
  xlab("Year") + ylab(expression(paste("Conductivity" ," (", mu, "S cm"^-1*")")))+
  theme_bw(base_size = 14)+
  theme(axis.text.y = element_text(face = "bold"))+
  theme(axis.title.y = element_text(face = "bold")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank()) #+
  
p_anCL

p_anCL <- p_anCL +
  geom_rect(data = phases, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = phase),
            alpha = 0.2,  inherit.aes = FALSE) +                 # Transparency for shading
  scale_fill_manual(values = c("Warm" = "red", "Cold" = "blue")) 
p_anCL 


pCL<- df_c %>%
  ggplot(aes(x = year, y = Conductivity)) +
  geom_line()+
  ggtitle(expression(paste("Conductivity")))+
  xlab("Year") + ylab(expression(paste("Conductivity" ," (", mu, "S cm"^-1*")")))+
  theme_bw(base_size = 14)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank()) 
pCL



##----------------------------------------------------------------------------##
## 1. predicted TDS
##----------------------------------------------------------------------------##
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
    date = "Date"
  )

#Remove rows without Chla data
tds <- tds[complete.cases(tds$TDS),]

tds<- tds %>% select(date, TDS)
# 284 obs


df_t<- left_join(df_c, tds, by="date")

# Fit linear model
lm_model <- lm(TDS ~ Conductivity, data = df_t)
summary(lm_model)


# Predict TDS where it's missing
df_t <- df_t %>%
  mutate(Predicted_TDS = ifelse(is.na(TDS), predict(lm_model, newdata = df_t), TDS))


ann <-aggregate(df_t$Conductivity,  
                by=list(df_t$year),  
                FUN=mean,
                na.rm=TRUE) 
ann<- ann %>% dplyr::rename(mean = x)


ann1 <-aggregate(df_t$Conductivity,  
                 by=list(df_t$year),  
                 FUN = max,
                 na.rm=TRUE)
ann1 <- ann1 %>% dplyr::rename(max = x)


ann2 <-aggregate(df_t$Conductivity,  
                 by=list(df_t$year),  
                 FUN = min,
                 na.rm=TRUE)
ann2 <- ann2 %>% dplyr::rename(min = x)

ann_tds<- left_join(ann, ann1,  by= "Group.1")
ann_tds<- left_join(ann_tds, ann2, by= "Group.1")


p_anTD<- ann_tds %>%
  ggplot(aes(x = Group.1, y = mean)) +
  geom_point()+
  geom_ribbon(aes(ymin = min, ymax = max), fill = '#165459B2', alpha = 1/2) +
  geom_line()+
  labs(fill = "SOI_PDO Phase")+
  ggtitle(expression(paste("Total dissolved solids")))+
  xlab("Year") + ylab(expression(paste("TDS" ," (mg L"^-1*")")))+
  theme_bw(base_size = 14)+
  theme(axis.text.y = element_text(face = "bold"))+
  theme(axis.title.y = element_text(face = "bold")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank()) #+
#scale_y_log10()
p_anTD

p_anTD <- p_anTD+ 
  geom_rect(data = in_phase, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = phase),
            alpha = 0.2,  inherit.aes = FALSE) +                 # Transparency for shading
  scale_fill_manual(values = c("Warm" = "red", "Cold" = "blue")) 

p_anTD

# simpler plot
pTD<- df_t %>%
  ggplot(aes(x = year, y = Conductivity)) +
  geom_line()+
  ggtitle(expression(paste("Conductivity")))+
  xlab("Year") + ylab(expression(paste("Conductivity" ," (", mu, "S cm"^-1*")")))+
  theme_bw(base_size = 14)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank()) 
pTD

##----------------------------------------------------------------------------##
## 1. turbidity
##----------------------------------------------------------------------------##


## Read in turb data 
df_tu <- read_csv("data/bpgamdataCLEAN_Turb.csv")

df_tu$Date<- as.Date(df_tu$Date, "%m/%d/%Y")

# add in Year and nMonth for numeric month and a proper Date class

df_tu <- mutate(df_tu,
             year = as.numeric(format(Date,'%Y')),
             DOY = as.numeric(format(Date,'%j')),
             nMonth = as.numeric(format(Date,'%m')),
             week = as.numeric(format(Date, '%W')))%>% 
  filter(year %in% c(1990:2022))
# 974 obs


## order the data
df_tu <- df_tu[with(df_tu, order(Date)),]


# rename for ease of use
df_tu <- df_tu%>%
  dplyr::rename(
    date = "Date",
    turb = "Turbidity")


#Remove rows without turb data
df_tu <- df_tu[complete.cases(df_tu$turb),]
# 974 obs


ann <-aggregate(df_tu$turb,  
                by=list(df_tu$year),  
                FUN=mean,
                na.rm=TRUE) 
ann<- ann %>% dplyr::rename(mean = x)


ann1 <-aggregate(df_tu$turb,  
                 by=list(df_tu$year),  
                 FUN = max,
                 na.rm=TRUE)
ann1 <- ann1 %>% dplyr::rename(max = x)


ann2 <-aggregate(df_tu$turb,  
                 by=list(df_tu$year),  
                 FUN = min,
                 na.rm=TRUE)
ann2 <- ann2 %>% dplyr::rename(min = x)

anntu<- left_join(ann, ann1,  by= "Group.1")
anntu<- left_join(anntu, ann2, by= "Group.1")


p_antu<- anntu %>%
  ggplot(aes(x = Group.1, y = mean)) +
  geom_point()+
  geom_ribbon(aes(ymin = min, ymax = max), fill = '#165459B2', alpha = 1/2) +
  geom_line()+
  ggtitle(expression(paste("Turbidity")))+
  xlab("Year") + ylab(expression(paste("Turbidity (NTU)")))+
  labs(fill = "SOI_PDO Phase")+
  theme_bw(base_size = 14)+
  theme(axis.text.y = element_text(face = "bold"))+
  theme(axis.title.y = element_text(face = "bold")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank()) #+
#scale_y_log10()
p_antu

p_antu <- p_antu + 
  geom_rect(data = in_phase, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = phase),
            alpha = 0.2,  inherit.aes = FALSE) +                 # Transparency for shading
  scale_fill_manual(values = c("Warm" = "red", "Cold" = "blue")) 

p_antu

ptu<- df_tu %>%
  ggplot(aes(x = year, y = turb)) +
  geom_line()+
  ggtitle(expression(paste("Turbidity")))+
  xlab("Year") + ylab(expression(paste("Turbidity (NTU)")))+
  theme_bw(base_size = 14)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank()) 
ptu
#---------------------------------------- plot annual average DOC


df_d <- read_csv("data/bpgamdataCLEAN_DOC.csv")


df_d$Date<- as.Date(df_d$Date, "%m/%d/%Y")

# add in Year and nMonth for numeric month and a proper Date class

df_d <- mutate(df_d,
             year = as.numeric(format(Date,'%Y')),
             DOY = as.numeric(format(Date,'%j')),
             nMonth = as.numeric(format(Date,'%m')),
             week = as.numeric(format(Date, '%W')))%>% 
  filter(year %in% c(1990:2022))

# 1670 obs


## order the data
df_d <- df_d[with(df_d, order(Date)),]


# rename for ease of use
df_d <- df_d%>%
  dplyr::rename(
    date = "Date"
  )

#Remove rows without DOC data
df_d1 <- df_d[complete.cases(df_d$DOC),]
# 1655 obs



# plot annual chla totals
ann <-aggregate(df_d1$DOC,  
                by=list(df_d1$year),  
                FUN=mean,
                na.rm=TRUE) 
ann<- ann %>% dplyr::rename(mean = x)

ann1 <-aggregate(df_d1$DOC,  
                 by=list(df_d1$year),  
                 FUN = max,
                 na.rm=TRUE)
ann1 <- ann1 %>% dplyr::rename(max = x)

ann2 <-aggregate(df_d1$DOC,  
                 by=list(df_d1$year),  
                 FUN = min,
                 na.rm=TRUE)
ann2 <- ann2 %>% dplyr::rename(min = x)

annDOC<- left_join(ann, ann1,  by= "Group.1")
annDOC<- left_join(annDOC, ann2, by= "Group.1")

p_anDOC<- annDOC %>%
  ggplot(aes(x = Group.1, y = mean)) +
  geom_point()+
  geom_ribbon(aes(ymin = min, ymax = max), fill = '#165459B2', alpha = 1/2) +
  geom_line()+
  ggtitle("DOC")+
  xlab("Year") + ylab(expression(paste("DOC  ", "(mg L"^-1*")")))+
  theme_bw(base_size = 14)+
  labs(fill = "SOI_PDO Phase")+
  theme(axis.text.y = element_text(face = "bold"))+
  theme(axis.title.y = element_text(face = "bold")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank()) #+
#scale_y_log10()
p_anDOC

p_anDOC <- p_anDOC + 
  geom_rect(data = in_phase, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = phase),
            alpha = 0.2,  inherit.aes = FALSE) +                 # Transparency for shading
  scale_fill_manual(values = c("Warm" = "red", "Cold" = "blue")) 
p_anDOC


pDOC<- df_d1 %>%
  ggplot(aes(x = year, y = DOC)) +
  geom_line()+
  ggtitle("DOC")+
  xlab("Year") + ylab(expression(paste("DOC  ", "(mg L"^-1*")")))+
  theme_bw(base_size = 14)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank()) 
pDOC


##---------- plot annual TON


df_o <- read_csv("data/bpgamdataCLEAN_TON.csv")

df_o$Date<- as.Date(df_o$Date, "%m/%d/%Y")

# add in Year and nMonth for numeric month and a proper Date class

df_o <- mutate(df_o,
             year = as.numeric(format(Date,'%Y')),
             DOY = as.numeric(format(Date,'%j')),
             nMonth = as.numeric(format(Date,'%m')),
             week = as.numeric(format(Date, '%W')))%>% 
  filter(year %in% c(1990:2022))
# 1435 obs

## order the data
df_o <- df_o[with(df_o, order(Date)),]


# rename for ease of use
df_o <- df_o%>%
  dplyr::rename(
    date = "Date",
    Odour = "Odour_TON")


#Remove rows without Odour data
df_o <- df_o[complete.cases(df_o$Odour),]

ann <-aggregate(df_o$Odour,  
                by=list(df_o$year),  
                FUN=mean,
                na.rm=TRUE) 
ann<- ann %>% dplyr::rename(mean = x)

ann1 <-aggregate(df_o$Odour,  
                 by=list(df_o$year),  
                 FUN = max,
                 na.rm=TRUE)
ann1 <- ann1 %>% dplyr::rename(max = x)

ann2 <-aggregate(df_o$Odour,  
                 by=list(df_o$year),  
                 FUN = min,
                 na.rm=TRUE)
ann2 <- ann2 %>% dplyr::rename(min = x)

annTON<- left_join(ann, ann1,  by= "Group.1")
annTON<- left_join(annTON, ann2, by= "Group.1")

p_anTON<- annTON %>%
  ggplot(aes(x = Group.1, y = mean)) +
  geom_point()+
  geom_ribbon(aes(ymin = min, ymax = max), fill = '#165459B2', alpha = 1/2) +
  geom_line()+
  ggtitle("Odour")+
  xlab("Year") + ylab(expression(paste("Odour (T.O.N.)")))+
  labs(fill = "SOI_PDO Phase")+
  theme_bw(base_size = 14)+
  theme(axis.text.y = element_text(face = "bold"))+
  theme(axis.title.y = element_text(face = "bold")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())
p_anTON

p_anTON <- p_anTON + 
  geom_rect(data = in_phase, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = phase),
            alpha = 0.2,  inherit.aes = FALSE) +                 # Transparency for shading
  scale_fill_manual(values = c("Warm" = "red", "Cold" = "blue")) 
p_anTON


pTON<- df_o %>%
  ggplot(aes(x = year, y = Odour)) +
  geom_line()+
  ggtitle("Odour")+
  xlab("Year") + ylab("Odour (T.O.N.)")+
  theme_bw(base_size = 14)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank()) 
pTON


#--------------------------------------------------- plot annual QLD totals


flow <- read_csv("Data/daily flow.csv")

flow<- flow %>%
  select(date, combined_05JG004.cms, SK05JG006.cms, RC_IC_cms)%>%
  dplyr::rename("Date" = date,
                QLD = "SK05JG006.cms",
                QWS = "RC_IC_cms",)
flow <- mutate(flow,
               year = as.numeric(format(Date,'%Y')),
               DOY = as.numeric(format(Date,'%j')),
               nMonth = as.numeric(format(Date,'%m')))%>% 
  filter(year %in% c(1990:2022))



ann <-aggregate(flow$QLD,  
                by=list(flow$year),  
                FUN=mean,
                na.rm=TRUE) 
ann<- ann %>% dplyr::rename(mean = x)

ann1 <-aggregate(flow$QLD,  
                 by=list(flow$year),  
                 FUN = max,
                 na.rm=TRUE)
ann1 <- ann1 %>% dplyr::rename(max = x)

ann2 <-aggregate(flow$QLD,  
                 by=list(flow$year),  
                 FUN = min,
                 na.rm=TRUE)
ann2 <- ann2 %>% dplyr::rename(min = x)

annLD<- left_join(ann, ann1,  by= "Group.1")
annLD<- left_join(annLD, ann2, by= "Group.1")

p_an<- annLD %>%
  ggplot(aes(x = Group.1, y = mean)) +
  geom_point()+
  geom_ribbon(aes(ymin = min, ymax = max), fill = '#165459B2', alpha = 1/2) +
  geom_line()+
  ggtitle("QLD")+
  xlab("Year") + ylab(expression(paste(
    "QLD (",
    m^3,  s^-1,
    ")")))+
  labs(fill = "SOI_PDO Phase")+
  theme_bw(base_size = 14)+
  theme(axis.text.y = element_text(face = "bold"))+
  theme(axis.title.y = element_text(face = "bold"))+
theme(axis.text.x = element_text(face = "bold", size=14))+
  theme(axis.title.x = element_text(face = "bold", size=16))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())
p_an

p_an <- p_an+ 
  geom_rect(data = in_phase, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = phase),
            alpha = 0.2,  inherit.aes = FALSE) +                 # Transparency for shading
  scale_fill_manual(values = c("Warm" = "red", "Cold" = "blue")) 

p_an

#------------------------------------------- plot annual total WS contributions


ann <-aggregate(flow$QWS,  
                by=list(flow$year),  
                FUN=mean,
                na.rm=TRUE) 
ann<- ann %>% dplyr::rename(mean = x)

ann1 <-aggregate(flow$QWS,  
                 by=list(flow$year),  
                 FUN = max,
                 na.rm=TRUE)
ann1 <- ann1 %>% dplyr::rename(max = x)

ann2 <-aggregate(flow$QWS,  
                 by=list(flow$year),  
                 FUN = min,
                 na.rm=TRUE)
ann2 <- ann2 %>% dplyr::rename(min = x)

annWS<- left_join(ann, ann1,  by= "Group.1")
annWS<- left_join(annWS, ann2, by= "Group.1")


p_anWS<- annWS %>%
  ggplot(aes(x = Group.1, y = mean)) +
  geom_point()+
  geom_ribbon(aes(ymin = min, ymax = max), fill = '#165459B2', alpha = 1/2) +
  geom_line()+
  ggtitle("QWS")+
  xlab("Year") + ylab(expression(paste(
    "QWS (",
    m^3, s^-1,
    ")")))+
  labs(fill = "SOI_PDO Phase")+
  theme_bw(base_size = 14)+
  theme(axis.text.y = element_text(face = "bold"))+
  theme(axis.title.y = element_text(face = "bold"))+
  theme(axis.text.x = element_text(face = "bold", size=14))+
  theme(axis.title.x = element_text(face = "bold", size=16)) #+
  #theme(axis.title.x = element_blank(),
        #axis.text.x = element_blank())
p_anWS

p_anWS <- p_anWS + 
  geom_rect(data = in_phase, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = phase),
            alpha = 0.2,  inherit.aes = FALSE) +                 # Transparency for shading
  scale_fill_manual(values = c("Warm" = "red", "Cold" = "blue")) 
p_anWS

## put it all together

p_all<- (p_anDOC / p_anTON/ p_anTD / p_antu)+ plot_layout(guides = "collect") 
#&  theme(legend.position = "top")
p_all

# p_all<- plot_grid(p_antu, p_anDOC, p_anTON, p_anCL, ncol = 1, align = "v")
# p_all

ggsave('output/annual totals DW indicators MAX and MIN SOI_PDO phase.png', p_all, height = 15, width  = 12, dpi = 320)

# Combine plots with one shared legend

p_all<- (p_anDOC / p_anTON/ p_anTD /p_antu/p_an/ p_anWS )+ plot_layout(guides = "collect") &  theme(legend.position = "bottom")
p_all

#  p_anSRP/ p_anDIN/ 
# p_all<- plot_grid(p_antu, p_anDOC, p_anTON, p_anCL, p_an, p_anWS, ncol = 1, align = "v")
# p_all

ggsave('output/annual totals nutreitns added MAX and MIN SOI_PDO phase .png', p_all, height = 15, width  = 12, dpi = 320)



## adding predictors.....

#Remove rows without complete SRP data
df_tus <- df_tu[complete.cases(df_tu$SRP_ug.L),]

ann <-aggregate(df_tus$SRP_ug.L,  
                by=list(df_tus$year),  
                FUN=mean,
                na.rm=TRUE) 
ann<- ann %>% dplyr::rename(mean = x)

ann1 <-aggregate(df_tus$SRP_ug.L,  
                 by=list(df_tus$year),  
                 FUN = max,
                 na.rm=TRUE)
ann1 <- ann1 %>% dplyr::rename(max = x)

ann2 <-aggregate(df_tus$SRP_ug.L,  
                 by=list(df_tus$year),  
                 FUN = min,
                 na.rm=TRUE)
ann2 <- ann2 %>% dplyr::rename(min = x)

annSRP<- left_join(ann, ann1,  by= "Group.1")
annSRP<- left_join(annSRP, ann2, by= "Group.1")

p_anSRP<- annSRP %>%
  ggplot(aes(x = Group.1, y = mean)) +
  geom_point()+
  geom_ribbon(aes(ymin = min, ymax = max), fill = '#165459B2', alpha = 1/2) +
  geom_line()+
  ggtitle("SRP")+
  xlab("Year") + ylab(expression(paste("SRP ", "(", mu, "g L"^-1*")")))+
  labs(fill = "SOI_PDO Phase")+
  theme_bw(base_size = 14)+
  theme(axis.text.y = element_text(face = "bold"))+
  theme(axis.title.y = element_text(face = "bold"))+
  theme(axis.text.x = element_text(face = "bold", size=14))+
  theme(axis.title.x = element_text(face = "bold", size=16)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())
p_anSRP

p_anSRP <- p_anSRP + 
  geom_rect(data = in_phase, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = phase),
            alpha = 0.2,  inherit.aes = FALSE) +                 # Transparency for shading
  scale_fill_manual(values = c("Warm" = "red", "Cold" = "blue")) 
p_anSRP

#------------------------------------------- plot annual TP loads

ann <-aggregate(df_ds$TP,  
                by=list(df_d$year),  
                FUN=mean,
                na.rm=TRUE) 
ann<- ann %>% dplyr::rename(mean = x)

ann1 <-aggregate(df_ds$TP,  
                 by=list(df_d$year),  
                 FUN = max,
                 na.rm=TRUE)
ann1 <- ann1 %>% dplyr::rename(max = x)

ann2 <-aggregate(df_ds$TP,  
                 by=list(df_d$year),  
                 FUN = min,
                 na.rm=TRUE)
ann2 <- ann2 %>% dplyr::rename(min = x)

annTP<- left_join(ann, ann1,  by= "Group.1")
annTP<- left_join(annTP, ann2, by= "Group.1")

p_anTP<- annTP %>%
  ggplot(aes(x = Group.1, y = mean)) +
  geom_point()+
  geom_ribbon(aes(ymin = min, ymax = max), fill = '#165459B2', alpha = 1/2) +
  geom_line()+
  theme_bw(base_size = 14)+
  theme(axis.text.y = element_text(face = "bold"))+
  theme(axis.title.y = element_text(face = "bold")) +
  ggtitle("TP")+
  xlab("Year") + ylab(expression(paste("TP  ", "(", mu, "g L"^-1*")")))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())
p_anTP

#ggsave('output/TP Annual total.png', p_anTP, height = 8, width  = 10)



#---------------------------------------------- plot annual DIN totals



df_tu$DIN <- rowSums(df_tu[,c("NO3_mg.L", "NH3_mg.L")], na.rm=TRUE)

df_tuN <- df_tu[complete.cases(df_tu$DIN),]

ann <-aggregate(df_tu$DIN,  
                by=list(df_tuN$year),  
                FUN=mean,
                na.rm=TRUE) 
ann<- ann %>% dplyr::rename(mean = x)

ann1 <-aggregate(df_tu$DIN,  
                 by=list(df_tuN$year),  
                 FUN = max,
                 na.rm=TRUE)
ann1 <- ann1 %>% dplyr::rename(max = x)

ann2 <-aggregate(df_tu$DIN,  
                 by=list(df_tuN$year),  
                 FUN = min,
                 na.rm=TRUE)
ann2 <- ann2 %>% dplyr::rename(min = x)

annDIN<- left_join(ann, ann1,  by= "Group.1")
annDIN<- left_join(annDIN, ann2, by= "Group.1")


p_anDIN<- annDIN %>%
  ggplot(aes(x = Group.1, y = mean)) +
  geom_point()+
  geom_ribbon(aes(ymin = min, ymax = max), fill = '#165459B2', alpha = 1/2) +
  geom_line()+
  ggtitle("DIN")+
  xlab("Year") + ylab(expression(paste("DIN ", "(mg L"^-1*")")))+
  labs(fill = "SOI_PDO Phase")+
  theme_bw(base_size = 14)+
  theme(axis.text.y = element_text(face = "bold"))+
  theme(axis.title.y = element_text(face = "bold"))+
  theme(axis.text.x = element_text(face = "bold", size=14))+
  theme(axis.title.x = element_text(face = "bold", size=16))  +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())
p_anDIN

p_anDIN <- p_anDIN + 
  geom_rect(data = in_phase, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = phase),
            alpha = 0.2,  inherit.aes = FALSE) +                 # Transparency for shading
  scale_fill_manual(values = c("Warm" = "red", "Cold" = "blue")) 
p_anDIN

#--------------------------------plot annual TN totals

df1$TN <- rowSums(df1[,c("NO3", "NH3", "orgN")], na.rm=TRUE)

df1 <- df1[complete.cases(df1$TN),]

ann <-aggregate(df1$TN,  
                by=list(df1$year),  
                FUN=mean,
                na.rm=TRUE) 
ann<- ann %>% dplyr::rename(mean = x)

ann1 <-aggregate(df1$TN,  
                 by=list(df1$year),  
                 FUN = max,
                 na.rm=TRUE)
ann1 <- ann1 %>%dplyr::rename(max = x)

ann2 <-aggregate(df1$TN,  
                 by=list(df1$year),  
                 FUN = min,
                 na.rm=TRUE)
ann2 <- ann2 %>% dplyr::rename(min = x)

annTN<- left_join(ann, ann1,  by= "Group.1")
annTN<- left_join(annTN, ann2, by= "Group.1")

p_anTN<- annTN %>%
  ggplot(aes(x = Group.1, y = mean)) +
  geom_point()+
  geom_ribbon(aes(ymin = min, ymax = max), fill = '#165459B2', alpha = 1/2) +
  geom_line()+
  ggtitle("TN")+
  xlab("Year") + ylab(expression(paste("TN ", "(mg L"^-1*")")))+
  theme_bw(base_size = 14)+
  theme(axis.text.y = element_text(face = "bold"))+
  theme(axis.title.y = element_text(face = "bold")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())
p_anTN


p_all<- plot_grid(p_antu, p_anDOC, p_anTON, p_anCL,  p_anSRP, p_anTP, p_anDIN, p_anTN, p_an,p_anWS, ncol = 1, align = "v")
p_all

ggsave('output/annual totals FLOW nutrients added MAX and MIN.png', p_all, height = 18, width  = 16, dpi = 320)



# Combine plots with one shared legend

p_all<- (p_anSRP / p_anDIN/ p_an/ p_anWS )+ plot_layout(guides = "collect") 
#&  theme(legend.position = "top")
p_all


ggsave('output/annual totals predictors SOI_PDO phase.png', p_all, height = 15, width  = 12, dpi = 320)


###----------------------------------------------------------------------------
### Decadal averages-water temperature

df$date <- as.POSIXct(df$date)

p1<- df %>%
  group_by(year) %>%
  summarize(max = max(W_temp, na.rm = TRUE), Date = median(date)) %>% {
    ggplot(., aes(Date, max)) +
      geom_line(color = 'gray') +
      geom_point(color = 'gray75') +
      geom_textsegment(aes(x = as.POSIXct('1984-01-01'), 
                           xend = as.POSIXct('1989-12-31'),
                           y = mean(max), yend = mean(max), color = '1980s',
                           label = '1980s'), vjust = -0.2, size = 6,
                       data = .[.$Date < as.POSIXct('1990-01-01'),], linetype = 2) +
      geom_textsegment(aes(x = as.POSIXct('1990-01-01'), 
                           xend = as.POSIXct('1999-12-31'),
                           y = mean(max), yend = mean(max), color = '1990s',
                           label = '1990s'), vjust = -0.2, size = 6,
                       data = .[.$Date < as.POSIXct('2000-01-01'),], linetype = 2) +
      geom_textsegment(aes(x = as.POSIXct('2000-01-01'), 
                           xend = as.POSIXct('2009-12-31'),
                           y = mean(max), yend = mean(max), color = '2000s',
                           label = '2000s'), vjust = -0.2, size = 6,
                       data = .[.$Date < as.POSIXct('2010-01-01') &
                                  .$Date > as.POSIXct('1999-12-31'),], linetype = 2) +
      geom_textsegment(aes(x = as.POSIXct('2010-01-01'), 
                           xend = as.POSIXct('2019-12-31'),
                           y = mean(max), yend = mean(max), color = '2010s',
                           label = '2010s'), vjust = -0.2, size = 6,
                       data = .[.$Date < as.POSIXct('2020-01-01') &
                                  .$Date > as.POSIXct('2009-12-31'),], linetype = 2) +
      geom_textsegment(aes(x = as.POSIXct('2020-01-01'), 
                           xend = as.POSIXct('2022-12-31'),
                           y = mean(max), yend = mean(max), color = '2020s',
                           label = '2020s'), vjust = -0.2, size = 6,
                       data = .[.$Date < as.POSIXct('2022-01-01') &
                                  .$Date > as.POSIXct('2019-12-31'),], linetype = 2) +
      theme_bw(base_size = 14) +
      theme(legend.position = 'none') +
      labs(title = 'Max water temperature', x = 'Year', y = expression(paste("Annual max water temperature (  " , degree*C,")")))
  }

p1<- p1 + scale_color_cosmic("hallmarks_dark")
p1

ggsave('output/max annual water temp.png', p1, height = 8, width  = 10)


p2<- df %>%
  group_by(year) %>%
  summarize(max = mean(W_temp, na.rm = TRUE), Date = median(date)) %>% {
    ggplot(., aes(Date, max)) +
      geom_line(color = 'gray') +
      geom_point(color = 'gray75') +
      geom_textsegment(aes(x = as.POSIXct('1984-01-01'), 
                           xend = as.POSIXct('1989-12-31'),
                           y = mean(max), yend = mean(max), color = '1980s',
                           label = '1980s'), vjust = -0.2, size = 6,
                       data = .[.$Date < as.POSIXct('1990-01-01'),], linetype = 2) +
      geom_textsegment(aes(x = as.POSIXct('1990-01-01'), 
                           xend = as.POSIXct('1999-12-31'),
                           y = mean(max), yend = mean(max), color = '1990s',
                           label = '1990s'), vjust = -0.2, size = 6,
                       data = .[.$Date < as.POSIXct('2000-01-01'),], linetype = 2) +
      geom_textsegment(aes(x = as.POSIXct('2000-01-01'), 
                           xend = as.POSIXct('2009-12-31'),
                           y = mean(max), yend = mean(max), color = '2000s',
                           label = '2000s'), vjust = -0.2, size = 6,
                       data = .[.$Date < as.POSIXct('2010-01-01') &
                                  .$Date > as.POSIXct('1999-12-31'),], linetype = 2) +
      geom_textsegment(aes(x = as.POSIXct('2010-01-01'), 
                           xend = as.POSIXct('2019-12-31'),
                           y = mean(max), yend = mean(max), color = '2010s',
                           label = '2010s'), vjust = -0.2, size = 6,
                       data = .[.$Date < as.POSIXct('2020-01-01') &
                                  .$Date > as.POSIXct('2009-12-31'),], linetype = 2) +
      geom_textsegment(aes(x = as.POSIXct('2020-01-01'), 
                           xend = as.POSIXct('2022-12-31'),
                           y = mean(max), yend = mean(max), color = '2020s',
                           label = '2020s'), vjust = -0.2, size = 6,
                       data = .[.$Date < as.POSIXct('2022-01-01') &
                                  .$Date > as.POSIXct('2019-12-31'),], linetype = 2) +
      theme_bw(base_size = 14) +
      theme(legend.position = 'none') +
      labs(title = 'Mean water temperature', x = 'Year', y = expression(paste("Annual mean water temperature (  " , degree*C,")")))
  }

p2<- p2 + scale_color_cosmic("hallmarks_dark")
p2

ggsave('output/mean annual water temp.png', p2, height = 8, width  = 10)


mean_annual <-aggregate(df$W_temp,  
                        by=list(df$year),  
                        FUN=mean,   
                        na.rm=TRUE) 

max <-aggregate(df$W_temp,  
                        by=list(df$year),  
                        FUN=max,   
                        na.rm=TRUE) 
# #------------------------------------------- plot annual TDS
# 
# # Read in TDS data 
# 
# tds <- read_csv("data/TDS.csv")
# 
# 
# tds$Date<- as.Date(tds$Date, "%m/%d/%Y")
# 
# tds <- mutate(tds,
#               year = as.numeric(format(Date,'%Y')),
#               DOY = as.numeric(format(Date,'%j')),
#               nMonth = as.numeric(format(Date,'%m')))%>% 
#   filter(year %in% c(1990:2022))
# 
# 
# ## order the data
# tds <- tds[with(tds, order(Date)),]
# 
# 
# # rename for ease of use
# tds <- tds%>%
#   dplyr::rename(
#     date = "Date")
# 
# 
# tds<-filter(tds, TDS != 0)
# 
# ann <-aggregate(tds$TDS,  
#                 by=list(tds$year),  
#                 FUN=mean,
#                 na.rm=TRUE) 
# ann<- ann %>% dplyr::rename(mean = x)
# 
# ann1 <-aggregate(tds$TDS,  
#                  by=list(tds$year),  
#                  FUN = max,
#                  na.rm=TRUE)
# ann1 <- ann1 %>% dplyr::rename(max = x)
# 
# ann2 <-aggregate(tds$TDS,  
#                  by=list(tds$year),  
#                  FUN = min,
#                  na.rm=TRUE)
# ann2 <- ann2 %>% dplyr::rename(min = x)
# 
# annTDS<- left_join(ann, ann1,  by= "Group.1")
# annTDS<- left_join(annTDS, ann2, by= "Group.1")
# 
# p_anTDS<- annTDS %>%
#   ggplot(aes(x = Group.1, y = mean)) +
#   geom_point()+
#   geom_ribbon(aes(ymin = min, ymax = max), fill = '#165459B2', alpha = 1/2) +
#   geom_line()+
#   theme_bw(base_size = 14)+
#   theme(axis.text.y = element_text(face = "bold"))+
#   theme(axis.title.y = element_text(face = "bold")) +
#   ggtitle("TDS")+
#   xlab("Year") + ylab(expression(paste("TDS  ", "(mg L"^-1*")")))+
#   xlim(1990,2022) +
#   theme_bw(base_size = 14)+
#   theme(axis.text.y = element_text(face = "bold"))+
#   theme(axis.title.y = element_text(face = "bold"))+
#   theme(axis.text.x = element_text(face = "bold", size=14))+
#   theme(axis.title.x = element_text(face = "bold", size=16)) +
#  theme(axis.title.x = element_blank(),
#        axis.text.x = element_blank())
# p_anTDS
# 
# #ggsave('output/TDS Annual total.png', p_anTDS, height = 8, width  = 10)
# 
# pTDS<- tds %>%
#   ggplot(aes(x = year, y = TDS)) +
#   geom_line()+
#   xlim(1990,2022) +
#   ggtitle("TDS")+
#   xlab("Year") + ylab(expression(paste("TDS  ", "(mg L"^-1*")")))+
#   theme_bw(base_size = 14)+
#   theme(axis.text.y = element_text(face = "bold"))+
#   theme(axis.title.y = element_text(face = "bold"))+
#   theme(axis.text.x = element_text(face = "bold", size=14))+
#   theme(axis.title.x = element_text(face = "bold", size=16))
# pTDS
