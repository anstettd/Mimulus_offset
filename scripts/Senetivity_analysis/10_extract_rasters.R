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
lambda <- read_csv("data/demography data/siteYear.lambda_responses_2010-2015.csv")
names(lambda)[names(lambda) == 'Region'] <- 'Region_demo'

#Import sensitivity rasters all for 1215
bf30 <- raster("data/genomic_data/offset_bf30.tif") #BF>30 only
env9 <- raster("data/genomic_data/offset_env9.tif") #new WZA
oldWZA <- raster("data/genomic_data/offset_oldWZA_snp_set.tif") #WZA snps only before update
old_snp_set <- raster("data/genomic_data/offset_old_snp_set.tif") #Original WZA snps + BF30
new_snp_set <- raster("data/genomic_data/offset_1215.tif") #New WZA snps + BF30

#Stack rasters
rasStack_gcc <- stack(bf30,env9,oldWZA,old_snp_set,new_snp_set)

#Define CRS
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS
#############################################################################################################
#Extract offset datapoints for demography populations

#Setup demography lat/long in correct format
demography_pop <- lambda %>% dplyr::select(Longitude,Latitude,Site)
coordinates(demography_pop)=cbind(demography_pop$Longitude,demography_pop$Latitude)
proj4string(demography_pop) <- EPSG4326

#Extract offset data from rasters
rasStack <- raster::extract(rasStack_gcc,demography_pop)
colnames(rasStack) <- c("offset_bf30","offset_env9","offset_oldWZA",
                        "offset_old_snp_set","offset_new_snp_set") #label conlumns

#Save offset data in dataframe
#timeseries_offset <- baseline_pop %>% filter(Paper_ID<13)
offset_pop <- cbind(lambda,rasStack)

#Add in site codes
demo_pop <- read_csv("data/genomic_data/demo_site_meta.csv")

#Merge new sites codes with offset_pop data
offset_pop_meta <- cbind(offset_pop,demo_pop)


#write_csv(offset_pop_meta,"data/genomic_data/offset_pop_sense.csv")


