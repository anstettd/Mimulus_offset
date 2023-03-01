#############################################################################################################
## Convert asc to TiFF files
## Water Year 2012-2015
## For Mean Annual Temperature
## Author Daniel Anstett
## 
## 
## Last Modified Feb 27, 2023
#############################################################################################################
#Import libraries
library(raster)
library(tidyverse)
library(sf)
library(rgdal)

#############################################################################################################
#Define CRS
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS


#Import asc file

#2012 Water Year
map_2012_10 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2011/PPT10.asc")
map_2012_11 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2011/PPT11.asc")
map_2012_12 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2011/PPT12.asc")

map_2012_01 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/PPT01.asc")
map_2012_02 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/PPT02.asc")
map_2012_03 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/PPT03.asc")
map_2012_04 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/PPT04.asc")
map_2012_05 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/PPT05.asc")
map_2012_06 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/PPT06.asc")
map_2012_07 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/PPT07.asc")
map_2012_08 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/PPT08.asc")
map_2012_09 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/PPT09.asc")

crs(map_2012_10) <- EPSG4326
crs(map_2012_11) <- EPSG4326
crs(map_2012_12) <- EPSG4326

crs(map_2012_01) <- EPSG4326
crs(map_2012_02) <- EPSG4326
crs(map_2012_03) <- EPSG4326
crs(map_2012_04) <- EPSG4326
crs(map_2012_05) <- EPSG4326
crs(map_2012_06) <- EPSG4326
crs(map_2012_07) <- EPSG4326
crs(map_2012_08) <- EPSG4326
crs(map_2012_09) <- EPSG4326

map_2012 <- (map_2012_10 + 
               map_2012_11 + 
               map_2012_12 + 
               map_2012_01 +
               map_2012_02 +
               map_2012_03 +
               map_2012_04 +
               map_2012_05 +
               map_2012_06 +
               map_2012_07 +
               map_2012_08 +
               map_2012_09)/12





#2013 Water Year
map_2013_10 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/PPT10.asc")
map_2013_11 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/PPT11.asc")
map_2013_12 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/PPT12.asc")

map_2013_01 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/PPT01.asc")
map_2013_02 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/PPT02.asc")
map_2013_03 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/PPT03.asc")
map_2013_04 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/PPT04.asc")
map_2013_05 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/PPT05.asc")
map_2013_06 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/PPT06.asc")
map_2013_07 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/PPT07.asc")
map_2013_08 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/PPT08.asc")
map_2013_09 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/PPT09.asc")

crs(map_2013_10) <- EPSG4326
crs(map_2013_11) <- EPSG4326
crs(map_2013_12) <- EPSG4326

crs(map_2013_01) <- EPSG4326
crs(map_2013_02) <- EPSG4326
crs(map_2013_03) <- EPSG4326
crs(map_2013_04) <- EPSG4326
crs(map_2013_05) <- EPSG4326
crs(map_2013_06) <- EPSG4326
crs(map_2013_07) <- EPSG4326
crs(map_2013_08) <- EPSG4326
crs(map_2013_09) <- EPSG4326

map_2013 <- (map_2013_10 + 
               map_2013_11 + 
               map_2013_12 + 
               map_2013_01 +
               map_2013_02 +
               map_2013_03 +
               map_2013_04 +
               map_2013_05 +
               map_2013_06 +
               map_2013_07 +
               map_2013_08 +
               map_2013_09)/12





#2014 Water Year
map_2014_10 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/PPT10.asc")
map_2014_11 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/PPT11.asc")
map_2014_12 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/PPT12.asc")

map_2014_01 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/PPT01.asc")
map_2014_02 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/PPT02.asc")
map_2014_03 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/PPT03.asc")
map_2014_04 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/PPT04.asc")
map_2014_05 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/PPT05.asc")
map_2014_06 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/PPT06.asc")
map_2014_07 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/PPT07.asc")
map_2014_08 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/PPT08.asc")
map_2014_09 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/PPT09.asc")

crs(map_2014_10) <- EPSG4326
crs(map_2014_11) <- EPSG4326
crs(map_2014_12) <- EPSG4326

crs(map_2014_01) <- EPSG4326
crs(map_2014_02) <- EPSG4326
crs(map_2014_03) <- EPSG4326
crs(map_2014_04) <- EPSG4326
crs(map_2014_05) <- EPSG4326
crs(map_2014_06) <- EPSG4326
crs(map_2014_07) <- EPSG4326
crs(map_2014_08) <- EPSG4326
crs(map_2014_09) <- EPSG4326

map_2014 <- (map_2014_10 + 
               map_2014_11 + 
               map_2014_12 + 
               map_2014_01 +
               map_2014_02 +
               map_2014_03 +
               map_2014_04 +
               map_2014_05 +
               map_2014_06 +
               map_2014_07 +
               map_2014_08 +
               map_2014_09)/12





#2015 Water Year
map_2015_10 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/PPT10.asc")
map_2015_11 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/PPT11.asc")
map_2015_12 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/PPT12.asc")

map_2015_01 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/PPT01.asc")
map_2015_02 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/PPT02.asc")
map_2015_03 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/PPT03.asc")
map_2015_04 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/PPT04.asc")
map_2015_05 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/PPT05.asc")
map_2015_06 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/PPT06.asc")
map_2015_07 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/PPT07.asc")
map_2015_08 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/PPT08.asc")
map_2015_09 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/PPT09.asc")

crs(map_2015_10) <- EPSG4326
crs(map_2015_11) <- EPSG4326
crs(map_2015_12) <- EPSG4326

crs(map_2015_01) <- EPSG4326
crs(map_2015_02) <- EPSG4326
crs(map_2015_03) <- EPSG4326
crs(map_2015_04) <- EPSG4326
crs(map_2015_05) <- EPSG4326
crs(map_2015_06) <- EPSG4326
crs(map_2015_07) <- EPSG4326
crs(map_2015_08) <- EPSG4326
crs(map_2015_09) <- EPSG4326

map_2015 <- (map_2015_10 + 
               map_2015_11 + 
               map_2015_12 + 
               map_2015_01 +
               map_2015_02 +
               map_2015_03 +
               map_2015_04 +
               map_2015_05 +
               map_2015_06 +
               map_2015_07 +
               map_2015_08 +
               map_2015_09)/12

# Average 2012 to 2015 rasters
map_1215 <- (map_2012 + map_2013 + map_2014 + map_2015)/4 

writeRaster(map_1215, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_1215/MAP.tif",formap="GTiff")


