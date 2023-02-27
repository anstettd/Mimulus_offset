#############################################################################################################
## Downscale asc raster from Climate NA to cover range extent
## Author Daniel Anstett
## 
## 
## Last Modified September 16, 2022
#############################################################################################################
#Import libraries
library(raster)
library(tidyverse)
library(sf)

#############################################################################################################
#Import asc file
#Range extent all of North America
na_asc <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/na800.asc")
#na_asc <- raster("C:/Users/anstett3/Downloads/Climatena_v730/InputFiles/na800.asc")

# Import M.cardinalis ensamble range extent as sf polygon
#c_range_old <- st_read("C:/Users/anstett3/Documents/Genomics/AM_Workshop/SDM/Output/c_range_2.shp")
#extent(c_range_old)

#Define range extent
ext <- as(extent(-124.6, -115.6, 31.4, 44.5), 'SpatialPolygons')
crs(ext) <- "+proj=longlat +datum=WGS84 +no_defs"



#Crop asc file by M. cardinalis range extent
#c_800 <- raster::crop(na_asc, extent(c_range_old))
c_800_all <- raster::crop(na_asc, ext)
#plot(c_800_all)


#Save c_800 raster as an asc file
writeRaster(c_800_all, "C:/Users/anstett3/Downloads/Climatena_v730/InputFiles/c_800_all.asc",format="ascii", overwrite = T)
writeRaster(c_800_all, "C:/Users/anstett3/Documents/Genomics/Large_files/c_800_all.asc",format="ascii", overwrite = T)

#Download all wanted climate layers from ClimateNA using ClimateNA_v7.30.exe, with c_800_all.asc as input file

