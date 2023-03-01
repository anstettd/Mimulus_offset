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
cmd_2012_10 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2011/CMD10.asc")
cmd_2012_11 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2011/CMD11.asc")
cmd_2012_12 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2011/CMD12.asc")

cmd_2012_01 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/CMD01.asc")
cmd_2012_02 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/CMD02.asc")
cmd_2012_03 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/CMD03.asc")
cmd_2012_04 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/CMD04.asc")
cmd_2012_05 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/CMD05.asc")
cmd_2012_06 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/CMD06.asc")
cmd_2012_07 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/CMD07.asc")
cmd_2012_08 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/CMD08.asc")
cmd_2012_09 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/CMD09.asc")

crs(cmd_2012_10) <- EPSG4326
crs(cmd_2012_11) <- EPSG4326
crs(cmd_2012_12) <- EPSG4326

crs(cmd_2012_01) <- EPSG4326
crs(cmd_2012_02) <- EPSG4326
crs(cmd_2012_03) <- EPSG4326
crs(cmd_2012_04) <- EPSG4326
crs(cmd_2012_05) <- EPSG4326
crs(cmd_2012_06) <- EPSG4326
crs(cmd_2012_07) <- EPSG4326
crs(cmd_2012_08) <- EPSG4326
crs(cmd_2012_09) <- EPSG4326

cmd_2012 <- (cmd_2012_10 + 
               cmd_2012_11 + 
               cmd_2012_12 + 
               cmd_2012_01 +
               cmd_2012_02 +
               cmd_2012_03 +
               cmd_2012_04 +
               cmd_2012_05 +
               cmd_2012_06 +
               cmd_2012_07 +
               cmd_2012_08 +
               cmd_2012_09)





#2013 Water Year
cmd_2013_10 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/CMD10.asc")
cmd_2013_11 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/CMD11.asc")
cmd_2013_12 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/CMD12.asc")

cmd_2013_01 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/CMD01.asc")
cmd_2013_02 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/CMD02.asc")
cmd_2013_03 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/CMD03.asc")
cmd_2013_04 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/CMD04.asc")
cmd_2013_05 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/CMD05.asc")
cmd_2013_06 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/CMD06.asc")
cmd_2013_07 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/CMD07.asc")
cmd_2013_08 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/CMD08.asc")
cmd_2013_09 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/CMD09.asc")

crs(cmd_2013_10) <- EPSG4326
crs(cmd_2013_11) <- EPSG4326
crs(cmd_2013_12) <- EPSG4326

crs(cmd_2013_01) <- EPSG4326
crs(cmd_2013_02) <- EPSG4326
crs(cmd_2013_03) <- EPSG4326
crs(cmd_2013_04) <- EPSG4326
crs(cmd_2013_05) <- EPSG4326
crs(cmd_2013_06) <- EPSG4326
crs(cmd_2013_07) <- EPSG4326
crs(cmd_2013_08) <- EPSG4326
crs(cmd_2013_09) <- EPSG4326

cmd_2013 <- (cmd_2013_10 + 
               cmd_2013_11 + 
               cmd_2013_12 + 
               cmd_2013_01 +
               cmd_2013_02 +
               cmd_2013_03 +
               cmd_2013_04 +
               cmd_2013_05 +
               cmd_2013_06 +
               cmd_2013_07 +
               cmd_2013_08 +
               cmd_2013_09)





#2014 Water Year
cmd_2014_10 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/CMD10.asc")
cmd_2014_11 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/CMD11.asc")
cmd_2014_12 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/CMD12.asc")

cmd_2014_01 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/CMD01.asc")
cmd_2014_02 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/CMD02.asc")
cmd_2014_03 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/CMD03.asc")
cmd_2014_04 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/CMD04.asc")
cmd_2014_05 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/CMD05.asc")
cmd_2014_06 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/CMD06.asc")
cmd_2014_07 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/CMD07.asc")
cmd_2014_08 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/CMD08.asc")
cmd_2014_09 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/CMD09.asc")

crs(cmd_2014_10) <- EPSG4326
crs(cmd_2014_11) <- EPSG4326
crs(cmd_2014_12) <- EPSG4326

crs(cmd_2014_01) <- EPSG4326
crs(cmd_2014_02) <- EPSG4326
crs(cmd_2014_03) <- EPSG4326
crs(cmd_2014_04) <- EPSG4326
crs(cmd_2014_05) <- EPSG4326
crs(cmd_2014_06) <- EPSG4326
crs(cmd_2014_07) <- EPSG4326
crs(cmd_2014_08) <- EPSG4326
crs(cmd_2014_09) <- EPSG4326

cmd_2014 <- (cmd_2014_10 + 
               cmd_2014_11 + 
               cmd_2014_12 + 
               cmd_2014_01 +
               cmd_2014_02 +
               cmd_2014_03 +
               cmd_2014_04 +
               cmd_2014_05 +
               cmd_2014_06 +
               cmd_2014_07 +
               cmd_2014_08 +
               cmd_2014_09)





#2015 Water Year
cmd_2015_10 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/CMD10.asc")
cmd_2015_11 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/CMD11.asc")
cmd_2015_12 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/CMD12.asc")

cmd_2015_01 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/CMD01.asc")
cmd_2015_02 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/CMD02.asc")
cmd_2015_03 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/CMD03.asc")
cmd_2015_04 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/CMD04.asc")
cmd_2015_05 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/CMD05.asc")
cmd_2015_06 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/CMD06.asc")
cmd_2015_07 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/CMD07.asc")
cmd_2015_08 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/CMD08.asc")
cmd_2015_09 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/CMD09.asc")

crs(cmd_2015_10) <- EPSG4326
crs(cmd_2015_11) <- EPSG4326
crs(cmd_2015_12) <- EPSG4326

crs(cmd_2015_01) <- EPSG4326
crs(cmd_2015_02) <- EPSG4326
crs(cmd_2015_03) <- EPSG4326
crs(cmd_2015_04) <- EPSG4326
crs(cmd_2015_05) <- EPSG4326
crs(cmd_2015_06) <- EPSG4326
crs(cmd_2015_07) <- EPSG4326
crs(cmd_2015_08) <- EPSG4326
crs(cmd_2015_09) <- EPSG4326

cmd_2015 <- (cmd_2015_10 + 
               cmd_2015_11 + 
               cmd_2015_12 + 
               cmd_2015_01 +
               cmd_2015_02 +
               cmd_2015_03 +
               cmd_2015_04 +
               cmd_2015_05 +
               cmd_2015_06 +
               cmd_2015_07 +
               cmd_2015_08 +
               cmd_2015_09)

# Average 2012 to 2015 rasters
cmd_1215 <- (cmd_2012 + cmd_2013 + cmd_2014 + cmd_2015)/4 

writeRaster(cmd_1215, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_1215/CMD.tif",forcmd="GTiff")


