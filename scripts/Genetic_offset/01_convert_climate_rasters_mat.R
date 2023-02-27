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
mat_2012_10 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2011/Tave10.asc")
mat_2012_11 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2011/Tave11.asc")
mat_2012_12 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2011/Tave12.asc")

mat_2012_01 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/Tave01.asc")
mat_2012_02 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/Tave02.asc")
mat_2012_03 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/Tave03.asc")
mat_2012_04 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/Tave04.asc")
mat_2012_05 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/Tave05.asc")
mat_2012_06 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/Tave06.asc")
mat_2012_07 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/Tave07.asc")
mat_2012_08 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/Tave08.asc")
mat_2012_09 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/Tave09.asc")

crs(mat_2012_10) <- EPSG4326
crs(mat_2012_11) <- EPSG4326
crs(mat_2012_12) <- EPSG4326

crs(mat_2012_01) <- EPSG4326
crs(mat_2012_02) <- EPSG4326
crs(mat_2012_03) <- EPSG4326
crs(mat_2012_04) <- EPSG4326
crs(mat_2012_05) <- EPSG4326
crs(mat_2012_06) <- EPSG4326
crs(mat_2012_07) <- EPSG4326
crs(mat_2012_08) <- EPSG4326
crs(mat_2012_09) <- EPSG4326

mat_2012 <- (mat_2012_10 + 
               mat_2012_11 + 
               mat_2012_12 + 
               mat_2012_01 +
               mat_2012_02 +
               mat_2012_03 +
               mat_2012_04 +
               mat_2012_05 +
               mat_2012_06 +
               mat_2012_07 +
               mat_2012_08 +
               mat_2012_09)/12





#2013 Water Year
mat_2013_10 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/Tave10.asc")
mat_2013_11 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/Tave11.asc")
mat_2013_12 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2012/Tave12.asc")

mat_2013_01 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/Tave01.asc")
mat_2013_02 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/Tave02.asc")
mat_2013_03 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/Tave03.asc")
mat_2013_04 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/Tave04.asc")
mat_2013_05 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/Tave05.asc")
mat_2013_06 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/Tave06.asc")
mat_2013_07 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/Tave07.asc")
mat_2013_08 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/Tave08.asc")
mat_2013_09 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/Tave09.asc")

crs(mat_2013_10) <- EPSG4326
crs(mat_2013_11) <- EPSG4326
crs(mat_2013_12) <- EPSG4326

crs(mat_2013_01) <- EPSG4326
crs(mat_2013_02) <- EPSG4326
crs(mat_2013_03) <- EPSG4326
crs(mat_2013_04) <- EPSG4326
crs(mat_2013_05) <- EPSG4326
crs(mat_2013_06) <- EPSG4326
crs(mat_2013_07) <- EPSG4326
crs(mat_2013_08) <- EPSG4326
crs(mat_2013_09) <- EPSG4326

mat_2013 <- (mat_2013_10 + 
               mat_2013_11 + 
               mat_2013_12 + 
               mat_2013_01 +
               mat_2013_02 +
               mat_2013_03 +
               mat_2013_04 +
               mat_2013_05 +
               mat_2013_06 +
               mat_2013_07 +
               mat_2013_08 +
               mat_2013_09)/12





#2014 Water Year
mat_2014_10 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/Tave10.asc")
mat_2014_11 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/Tave11.asc")
mat_2014_12 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2013/Tave12.asc")

mat_2014_01 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/Tave01.asc")
mat_2014_02 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/Tave02.asc")
mat_2014_03 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/Tave03.asc")
mat_2014_04 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/Tave04.asc")
mat_2014_05 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/Tave05.asc")
mat_2014_06 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/Tave06.asc")
mat_2014_07 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/Tave07.asc")
mat_2014_08 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/Tave08.asc")
mat_2014_09 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/Tave09.asc")

crs(mat_2014_10) <- EPSG4326
crs(mat_2014_11) <- EPSG4326
crs(mat_2014_12) <- EPSG4326

crs(mat_2014_01) <- EPSG4326
crs(mat_2014_02) <- EPSG4326
crs(mat_2014_03) <- EPSG4326
crs(mat_2014_04) <- EPSG4326
crs(mat_2014_05) <- EPSG4326
crs(mat_2014_06) <- EPSG4326
crs(mat_2014_07) <- EPSG4326
crs(mat_2014_08) <- EPSG4326
crs(mat_2014_09) <- EPSG4326

mat_2014 <- (mat_2014_10 + 
               mat_2014_11 + 
               mat_2014_12 + 
               mat_2014_01 +
               mat_2014_02 +
               mat_2014_03 +
               mat_2014_04 +
               mat_2014_05 +
               mat_2014_06 +
               mat_2014_07 +
               mat_2014_08 +
               mat_2014_09)/12





#2015 Water Year
mat_2015_10 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/Tave10.asc")
mat_2015_11 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/Tave11.asc")
mat_2015_12 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2014/Tave12.asc")

mat_2015_01 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/Tave01.asc")
mat_2015_02 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/Tave02.asc")
mat_2015_03 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/Tave03.asc")
mat_2015_04 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/Tave04.asc")
mat_2015_05 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/Tave05.asc")
mat_2015_06 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/Tave06.asc")
mat_2015_07 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/Tave07.asc")
mat_2015_08 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/Tave08.asc")
mat_2015_09 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Year_2015/Tave09.asc")

crs(mat_2015_10) <- EPSG4326
crs(mat_2015_11) <- EPSG4326
crs(mat_2015_12) <- EPSG4326

crs(mat_2015_01) <- EPSG4326
crs(mat_2015_02) <- EPSG4326
crs(mat_2015_03) <- EPSG4326
crs(mat_2015_04) <- EPSG4326
crs(mat_2015_05) <- EPSG4326
crs(mat_2015_06) <- EPSG4326
crs(mat_2015_07) <- EPSG4326
crs(mat_2015_08) <- EPSG4326
crs(mat_2015_09) <- EPSG4326

mat_2015 <- (mat_2015_10 + 
               mat_2015_11 + 
               mat_2015_12 + 
               mat_2015_01 +
               mat_2015_02 +
               mat_2015_03 +
               mat_2015_04 +
               mat_2015_05 +
               mat_2015_06 +
               mat_2015_07 +
               mat_2015_08 +
               mat_2015_09)/12

# Average 2012 to 2015 rasters
mat_1215 <- (mat_2012 + mat_2013 + mat_2014 + mat_2015)/4 

writeRaster(mat_1215, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_1215/MAT.tif",format="GTiff")


