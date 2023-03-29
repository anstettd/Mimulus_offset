#############################################################################################################
## Convert asc to TiFF files 
## Author Daniel Anstett
## 
## 
## Last Modified September 16, 2022
#############################################################################################################
#Import libraries
library(raster)
library(tidyverse)
library(sf)
library(rgdal)

########################################################################################################################
#Import files & Reproject

#2012
pas_2012 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2012Yearly/pas.asc")
ext_2012 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2012Yearly/ext.asc")
Tave_wt_2012 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2012Seasonal/Tave_wt.asc")
Tave_sm_2012 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2012Seasonal/Tave_sm.asc")
PPT_wt_2012 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2012Seasonal/PPT_wt.asc")
PPT_sm_2012 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2012Seasonal/PPT_sm.asc")

EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS
crs(pas_2012) <- EPSG4326
crs(ext_2012) <- EPSG4326
crs(Tave_wt_2012) <- EPSG4326
crs(Tave_sm_2012) <- EPSG4326
crs(PPT_wt_2012) <- EPSG4326
crs(PPT_sm_2012) <- EPSG4326


#2013
pas_2013 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2013Yearly/pas.asc")
ext_2013 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2013Yearly/ext.asc")
Tave_wt_2013 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2013Seasonal/Tave_wt.asc")
Tave_sm_2013 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2013Seasonal/Tave_sm.asc")
PPT_wt_2013 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2013Seasonal/PPT_wt.asc")
PPT_sm_2013 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2013Seasonal/PPT_sm.asc")

EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS
crs(pas_2013) <- EPSG4326
crs(ext_2013) <- EPSG4326
crs(Tave_wt_2013) <- EPSG4326
crs(Tave_sm_2013) <- EPSG4326
crs(PPT_wt_2013) <- EPSG4326
crs(PPT_sm_2013) <- EPSG4326


#2014
pas_2014 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2014Yearly/pas.asc")
ext_2014 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2014Yearly/ext.asc")
Tave_wt_2014 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2014Seasonal/Tave_wt.asc")
Tave_sm_2014 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2014Seasonal/Tave_sm.asc")
PPT_wt_2014 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2014Seasonal/PPT_wt.asc")
PPT_sm_2014 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2014Seasonal/PPT_sm.asc")

EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS
crs(pas_2014) <- EPSG4326
crs(ext_2014) <- EPSG4326
crs(Tave_wt_2014) <- EPSG4326
crs(Tave_sm_2014) <- EPSG4326
crs(PPT_wt_2014) <- EPSG4326
crs(PPT_sm_2014) <- EPSG4326


#2015
pas_2015 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2015Yearly/pas.asc")
ext_2015 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2015Yearly/ext.asc")
Tave_wt_2015 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2015Seasonal/Tave_wt.asc")
Tave_sm_2015 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2015Seasonal/Tave_sm.asc")
PPT_wt_2015 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2015Seasonal/PPT_wt.asc")
PPT_sm_2015 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year__2015Seasonal/PPT_sm.asc")

EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS
crs(pas_2015) <- EPSG4326
crs(ext_2015) <- EPSG4326
crs(Tave_wt_2015) <- EPSG4326
crs(Tave_sm_2015) <- EPSG4326
crs(PPT_wt_2015) <- EPSG4326
crs(PPT_sm_2015) <- EPSG4326


########################################################################################################################

# Average 2012 to 2015 rasters
pas_1215 <- (pas_2012 + pas_2013 + pas_2014 + pas_2015)/4 
ext_1215 <- (ext_2012 + ext_2013 + ext_2014 + ext_2015)/4 
Tave_wt_1215 <- (Tave_wt_2012 + Tave_wt_2013 + Tave_wt_2014 + Tave_wt_2015)/4 
Tave_sm_1215 <- (Tave_sm_2012 + Tave_sm_2013 + Tave_sm_2014 + Tave_sm_2015)/4 
PPT_wt_1215 <- (PPT_wt_2012 + PPT_wt_2013 + PPT_wt_2014 + PPT_wt_2015)/4 
PPT_sm_1215 <- (PPT_sm_2012 + PPT_sm_2013 + PPT_sm_2014 + PPT_sm_2015)/4 



#Export tif file
writeRaster(pas_1215, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_1215/PAS.tif",format="GTiff")
writeRaster(ext_1215, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_1215/EXT.tif",format="GTiff")
writeRaster(Tave_wt_1215, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_1215/Tave_wt.tif",format="GTiff")
writeRaster(Tave_sm_1215, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_1215/Tave_sm.tif",format="GTiff")
writeRaster(PPT_wt_1215, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_1215/PPT_wt.tif",format="GTiff")
writeRaster(PPT_sm_1215, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_1215/PPT_sm.tif",format="GTiff")




