##################################################################################
## Daniel Anstett
## CALCULATE CLIMATE FOR 2010-2016 
## Convert monthly temp & precip data October (of the previous year) to September Climate Data
## Average across calcualted years 1980 to 2009 
## 
## Last Modified January 22, 2020
###################################################################################

# Monthly climate data obtained from Climate NA 
# http://climatena.ca/

###################################################################################
#Load Libraries
library(tidyverse)

###################################################################################

#Import datasets and add year_actual variable
weather_1979 <- read.csv("Climate/timeseries_monthly_1979.csv", header=T)
weather_1980 <- read.csv("Climate/timeseries_monthly_1980.csv", header=T)
weather_1981 <- read.csv("Climate/timeseries_monthly_1981.csv", header=T)
weather_1982 <- read.csv("Climate/timeseries_monthly_1982.csv", header=T)
weather_1983 <- read.csv("Climate/timeseries_monthly_1983.csv", header=T)
weather_1984 <- read.csv("Climate/timeseries_monthly_1984.csv", header=T)
weather_1985 <- read.csv("Climate/timeseries_monthly_1985.csv", header=T)
weather_1986 <- read.csv("Climate/timeseries_monthly_1986.csv", header=T)
weather_1987 <- read.csv("Climate/timeseries_monthly_1987.csv", header=T)
weather_1988 <- read.csv("Climate/timeseries_monthly_1988.csv", header=T)
weather_1989 <- read.csv("Climate/timeseries_monthly_1989.csv", header=T)
weather_1990 <- read.csv("Climate/timeseries_monthly_1990.csv", header=T)
weather_1991 <- read.csv("Climate/timeseries_monthly_1991.csv", header=T)
weather_1992 <- read.csv("Climate/timeseries_monthly_1992.csv", header=T)
weather_1993 <- read.csv("Climate/timeseries_monthly_1993.csv", header=T)
weather_1994 <- read.csv("Climate/timeseries_monthly_1994.csv", header=T)
weather_1995 <- read.csv("Climate/timeseries_monthly_1995.csv", header=T)
weather_1996 <- read.csv("Climate/timeseries_monthly_1996.csv", header=T)
weather_1997 <- read.csv("Climate/timeseries_monthly_1997.csv", header=T)
weather_1998 <- read.csv("Climate/timeseries_monthly_1998.csv", header=T)
weather_1999 <- read.csv("Climate/timeseries_monthly_1999.csv", header=T)
weather_2000 <- read.csv("Climate/timeseries_monthly_2000.csv", header=T)
weather_2001 <- read.csv("Climate/timeseries_monthly_2001.csv", header=T)
weather_2002 <- read.csv("Climate/timeseries_monthly_2002.csv", header=T)
weather_2003 <- read.csv("Climate/timeseries_monthly_2003.csv", header=T)
weather_2004 <- read.csv("Climate/timeseries_monthly_2004.csv", header=T)
weather_2005 <- read.csv("Climate/timeseries_monthly_2005.csv", header=T)
weather_2006 <- read.csv("Climate/timeseries_monthly_2006.csv", header=T)
weather_2007 <- read.csv("Climate/timeseries_monthly_2007.csv", header=T)
weather_2008 <- read.csv("Climate/timeseries_monthly_2008.csv", header=T)
weather_2009 <- read.csv("Climate/timeseries_monthly_2009.csv", header=T)

###################################################################################
###################################################################################
# For loop to get Oct to Sept data from all needed years
impact_summary <- data.frame()
for(i in 1980:2009){
  impact<- eval(parse(text=(paste("weather",i-1,sep="_")))) %>% select(ID,ID2,Latitude,Longitude)
  impact<-cbind(impact,c(rep(i,12)))
  #MAT  
  impact_aT <- eval(parse(text=(paste("weather",i-1,sep="_")))) %>% select(Tave10,Tave11,Tave12)
  impact_bT <- eval(parse(text=(paste("weather",i,sep="_"))))%>% select(Tave01,Tave02,Tave03,Tave04,Tave05,
                                                                        Tave06,Tave07,Tave08,Tave09)
  impact_T <- cbind(impact_bT,impact_aT)
  MAT <- rowMeans(impact_T, na.rm = FALSE, dims = 1)
  impact <- cbind(impact,MAT)
  #MAP  
  impact_aT <- eval(parse(text=(paste("weather",i-1,sep="_")))) %>% select(PPT10,PPT11,PPT12)
  impact_bT <- eval(parse(text=(paste("weather",i,sep="_"))))%>% select(PPT01,PPT02,PPT03,PPT04,PPT05,PPT06,PPT07,PPT08,PPT09)
  impact_T <- cbind(impact_bT,impact_aT)
  MAP <- rowSums(impact_T, na.rm = FALSE, dims = 1)
  impact <- cbind(impact,MAP)
  #CMD  
  impact_aT <- eval(parse(text=(paste("weather",i-1,sep="_")))) %>% select(CMD10,CMD11,CMD12)
  impact_bT <- eval(parse(text=(paste("weather",i,sep="_"))))%>% select(CMD01,CMD02,CMD03,CMD04,CMD05,CMD06,CMD07,CMD08,CMD09)
  impact_T <- cbind(impact_bT,impact_aT)
  CMD <- rowSums(impact_T, na.rm = FALSE, dims = 1)
  impact <- cbind(impact,CMD)
  impact_summary <- rbind(impact_summary,impact)
}
colnames(impact_summary)[5]<-"hist_year"
###################################################################################

#Write modifed yearly data
write.csv(impact_summary,'Data/m_year.csv') #Export file

#Write 30-year average
climate_81_10 <- impact_summary %>% group_by(ID,ID2,Latitude,Longitude) %>% summarise_at(c("MAT", "MAP", "CMD"), mean, na.rm=TRUE)
write.csv(climate_81_10,'Data/climate.csv') #Export file
###################################################################################
