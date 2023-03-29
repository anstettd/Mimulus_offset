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

#############################################################################################################

#Import asc file

mat_8110 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Normal_1981_2010Y/mat.asc")
map_8110 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Normal_1981_2010Y/map.asc")
cmd_8110 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Normal_1981_2010Y/cmd.asc")
pas_8110 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Normal_1981_2010Y/pas.asc")
ext_8110 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Normal_1981_2010Y/ext.asc")
Tave_wt_8110 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Normal_1981_2010S/Tave_wt.asc")
Tave_sm_8110 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Normal_1981_2010S/Tave_sm.asc")
PPT_wt_8110 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Normal_1981_2010S/PPT_wt.asc")
PPT_sm_8110 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/Normal_1981_2010S/PPT_sm.asc")

mat_4170_45 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/13GCMs_ensemble_ssp245_2041-2070Y/mat.asc")
map_4170_45 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/13GCMs_ensemble_ssp245_2041-2070Y/map.asc")
cmd_4170_45 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/13GCMs_ensemble_ssp245_2041-2070Y/cmd.asc")
pas_4170_45 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/13GCMs_ensemble_ssp245_2041-2070Y/pas.asc")
ext_4170_45 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/13GCMs_ensemble_ssp245_2041-2070Y/ext.asc")
Tave_wt_4170_45 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/13GCMs_ensemble_ssp245_2041-2070S/Tave_wt.asc")
Tave_sm_4170_45 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/13GCMs_ensemble_ssp245_2041-2070S/Tave_sm.asc")
PPT_wt_4170_45 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/13GCMs_ensemble_ssp245_2041-2070S/PPT_wt.asc")
PPT_sm_4170_45 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/13GCMs_ensemble_ssp245_2041-2070S/PPT_sm.asc")

mat_4170_85 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/13GCMs_ensemble_ssp585_2041-2070Y/mat.asc")
map_4170_85 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/13GCMs_ensemble_ssp585_2041-2070Y/map.asc")
cmd_4170_85 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/13GCMs_ensemble_ssp585_2041-2070Y/cmd.asc")
pas_4170_85 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/13GCMs_ensemble_ssp585_2041-2070Y/pas.asc")
ext_4170_85 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/13GCMs_ensemble_ssp585_2041-2070Y/ext.asc")
Tave_wt_4170_85 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/13GCMs_ensemble_ssp585_2041-2070S/Tave_wt.asc")
Tave_sm_4170_85 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/13GCMs_ensemble_ssp585_2041-2070S/Tave_sm.asc")
PPT_wt_4170_85 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/13GCMs_ensemble_ssp585_2041-2070S/PPT_wt.asc")
PPT_sm_4170_85 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc_updated/13GCMs_ensemble_ssp585_2041-2070S/PPT_sm.asc")


#Reproject to WGS 1984 (EPSG4326)
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS

crs(mat_8110) <- EPSG4326
crs(map_8110) <- EPSG4326
crs(cmd_8110) <- EPSG4326
crs(pas_8110) <- EPSG4326
crs(ext_8110) <- EPSG4326
crs(Tave_wt_8110) <- EPSG4326
crs(Tave_sm_8110) <- EPSG4326
crs(PPT_wt_8110) <- EPSG4326
crs(PPT_sm_8110) <- EPSG4326

crs(mat_4170_45) <- EPSG4326
crs(map_4170_45) <- EPSG4326
crs(cmd_4170_45) <- EPSG4326
crs(pas_4170_45) <- EPSG4326
crs(ext_4170_45) <- EPSG4326
crs(Tave_wt_4170_45) <- EPSG4326
crs(Tave_sm_4170_45) <- EPSG4326
crs(PPT_wt_4170_45) <- EPSG4326
crs(PPT_sm_4170_45) <- EPSG4326

crs(mat_4170_85) <- EPSG4326
crs(map_4170_85) <- EPSG4326
crs(cmd_4170_85) <- EPSG4326
crs(pas_4170_85) <- EPSG4326
crs(ext_4170_85) <- EPSG4326
crs(Tave_wt_4170_85) <- EPSG4326
crs(Tave_sm_4170_85) <- EPSG4326
crs(PPT_wt_4170_85) <- EPSG4326
crs(PPT_sm_4170_85) <- EPSG4326



#Export tif file

#writeRaster(mat_8110, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_8110/MAT.tif",format="GTiff")
#writeRaster(map_8110, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_8110/MAP.tif",format="GTiff")
#writeRaster(cmd_8110, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_8110/CMD.tif",format="GTiff")
#writeRaster(pas_8110, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_8110/PAS.tif",format="GTiff")
#writeRaster(ext_8110, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_8110/EXT.tif",format="GTiff")
writeRaster(Tave_wt_8110, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_8110/Tave_wt.tif",format="GTiff")
writeRaster(Tave_sm_8110, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_8110/Tave_sm.tif",format="GTiff")
writeRaster(PPT_wt_8110, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_8110/PPT_wt.tif",format="GTiff")
writeRaster(PPT_sm_8110, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_8110/PPT_sm.tif",format="GTiff")


#writeRaster(mat_4170_45, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_4170_45/MAT.tif",format="GTiff")
#writeRaster(map_4170_45, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_4170_45/MAP.tif",format="GTiff")
#writeRaster(cmd_4170_45, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_4170_45/CMD.tif",format="GTiff")
#writeRaster(pas_4170_45, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_4170_45/PAS.tif",format="GTiff")
#writeRaster(ext_4170_45, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_4170_45/EXT.tif",format="GTiff")
writeRaster(Tave_wt_4170_45, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_4170_45/Tave_wt.tif",format="GTiff")
writeRaster(Tave_sm_4170_45, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_4170_45/Tave_sm.tif",format="GTiff")
writeRaster(PPT_wt_4170_45, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_4170_45/PPT_wt.tif",format="GTiff")
writeRaster(PPT_sm_4170_45, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_4170_45/PPT_sm.tif",format="GTiff")


#writeRaster(mat_4170_85, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_4170_85/MAT.tif",format="GTiff")
#writeRaster(map_4170_85, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_4170_85/MAP.tif",format="GTiff")
#writeRaster(cmd_4170_85, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_4170_85/CMD.tif",format="GTiff")
#writeRaster(pas_4170_85, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_4170_85/PAS.tif",format="GTiff")
#writeRaster(ext_4170_85, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_4170_85/EXT.tif",format="GTiff")
writeRaster(Tave_wt_4170_85, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_4170_85/Tave_wt.tif",format="GTiff")
writeRaster(Tave_sm_4170_85, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_4170_85/Tave_sm.tif",format="GTiff")
writeRaster(PPT_wt_4170_85, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_4170_85/PPT_wt.tif",format="GTiff")
writeRaster(PPT_sm_4170_85, "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_4170_85/PPT_sm.tif",format="GTiff")




