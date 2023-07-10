##################################################################################
## Extract offset values for each timeseires population
## Author Daniel Anstett
## 
## 
## Last Modified September 19, 2022
###################################################################################

#Library install and import
library(tidyverse) # for data manipulation
library(raster) # for handeling rasters
library(tidyverse)
library(cowplot)
library(sf)
library(rgdal)
library(rnaturalearth)
library(rnaturalearthdata)
#library(devtools)
#devtools::install_github("ropensci/rnaturalearthhires")

##############################################################################

## Import Demography Data
#lambda <- read_csv("data/demography data/siteYear.lambda_responses_2010-2015.csv")
#names(lambda)[names(lambda) == 'Region'] <- 'Region_demo'


#lambda <- read_csv("data/demography data/siteYear.lambda_slopes_2010-2015.csv")
#demo <- read_csv("data/demography data/siteYear.lambda_2010-2015.csv")
#Get Lat/Long for each site
#demo_unique <- demo %>% dplyr::select(Site,Latitude,Longitude) 
#demo_unique <- unique(demo_unique) %>% filter(Site!="Deer Creek") %>% filter(Site!="Mill Creek") 
#demo_unique <- demo_unique %>% dplyr::select(-Site)

genomic_pops <- read_csv("data/genomic_data/Baseline_Timeseries_pops_final2.csv")
timeseries <- genomic_pops %>% filter(Paper_ID<13)
timeseries <- timeseries[order(timeseries$Paper_ID),]


## Import offset rasters
mask_offset_1215 <- raster("data/genomic_data/offset_1215_beagle.tif") #2012-2015 offset
clim_diff_1215 <- raster("data/genomic_data/offset_climate.tif") #climate difference offset
mask_offset_45 <- raster("data/genomic_data/offset_4.5_beagle.tif") #future climte change SSE245
mask_offset_85 <- raster("data/genomic_data/offset_8.5_beagle.tif") #future climte change SSE585

#Stack rasters
rasStack_gcc <- stack(mask_offset_1215,clim_diff_1215,mask_offset_45,mask_offset_85)

#Define CRS
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS
#############################################################################################################
#Extract offset datapoints for demography populations

#Setup demography lat/long in correct format
timeseries_pop <- timeseries %>% dplyr::select(Long,Lat,Paper_ID)
coordinates(timeseries_pop)=cbind(timeseries_pop$Long,timeseries_pop$Lat)
proj4string(timeseries_pop) <- EPSG4326

#Extract offset data from rasters
rasStack <- raster::extract(rasStack_gcc,timeseries_pop)
colnames(rasStack) <- c("offset_1215","offset_climate","offset_SSP245","offset_SSP585") #label conlumns

#Save offset data in dataframe
#timeseries_offset <- baseline_pop %>% filter(Paper_ID<13)
offset_pop <- cbind(timeseries,rasStack)

#Add in site codes
#demo_pop <- read_csv("data/genomic_data/demo_site_meta.csv")

#Merge new sites codes with offset_pop data
#offset_pop_meta <- cbind(offset_pop,demo_pop)

#write_csv(offset_pop,"/Users/daniel_anstett/Dropbox/AM_Workshop/snp_change/data/offset_pop_timeseries_beagle.csv")


