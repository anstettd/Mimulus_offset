##################################################################################
## Daniel Anstett
## CALCULATE DROUGHT ANOMALIES
## Use modified weather and cliamte data 
## 
## Last Modified January 22, 2020
###################################################################################

###################################################################################
#Load Libraries
library(tidyverse)

###################################################################################
#Import climate data
wna <- read_csv("Data/climate.csv") %>%  #Oct to Setp data
  #Rename variables and take log of MAP
  select(Site=ID, MAT.clim=MAT,MAP.clim=MAP,CMD.clim=CMD) %>% mutate(log.MAP.clim = log10(MAP.clim))
wna$Site <- as.factor(wna$Site) #Make site a factor

#Import weather data
wna2 <- read_csv("Data/weather.csv") #Import

#Rename variables, take log of MAP
wna2 <- wna2 %>% #Selects MAT, MAP, CMD,
  select(ID_Year1,Latitude,Longitude,Elevation,MAT.weath=MAT,MAP.weath=MAP,CMD.weath=CMD) %>% 
  mutate(log.MAP.weath = log10(MAP.weath)) %>% #Take log of MAP
  separate(ID_Year1, into = c("Site", "Year"), sep = "_") #makes site/year variable
wna2$Site <- as.factor(wna2$Site) ; wna2$Year <- as.numeric(wna2$Year) #define variables
###################################################################################
# join climate and weather, calculate anomaly
wna_all <- left_join(wna2, wna, by="Site") %>% 
  mutate(CMD.anom = CMD.weath - CMD.clim, 
         MAT.anom = MAT.weath - MAT.clim,
         MAP.anom =  log.MAP.clim - log.MAP.weath) #reverted to log scale

write.csv(wna_all,'Data/wna_all.csv') #Export file
###################################################################################
