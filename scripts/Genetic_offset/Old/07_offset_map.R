##################################################################################
## Make gradient forest plots
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
library(tmap)
library(rnaturalearth)
#library(rnaturalearthdata)
#library(devtools)
#devtools::install_github("ropensci/rnaturalearthhires")
library(rnaturalearthhires)
library(rgeos)
library(RColorBrewer)
##############################################################################

## INPUTS

#Timeseries offset raster
mask_offset_1215 <- raster("data/genomic_data/offset_1215.tif") #pop data
clim_diff_1215 <- raster("data/genomic_data/offset_climate.tif") #pop data

#Future cliamte change offset raster
mask_offset_45_grain <- raster("data/genomic_data/offset_4.5_peakbf2_grain.tif") #pop data 
mask_offset_85_grain <- raster("data/genomic_data/offset_8.5_peakbf2_grain.tif") #pop data

#Labled as bf2, but is actually bf5

#Define CRS
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS

# California & Oregon Map Setup
states<-ne_states(country=c("canada","united states of america"),returnclass= "sf")
calo <- states %>%
  filter(name_en=="Oregon" | name_en=="California" | name_en=="Nevada")

#Baseline & Timeseries
#baseline_pop <- read_csv("data/genomic_data/paper_ID_site_select.csv")
#timeseries_pop <- baseline_pop %>% filter(Paper_ID<13) %>% dplyr::select(Long,Lat)
#timeseries_pop_sf <- st_as_sf(timeseries_pop,coords=c("Long","Lat"), crs=EPSG4326)
#baseline_pop <- baseline_pop  %>% dplyr::select(Long,Lat)
#baseline_pop_sf <- st_as_sf(baseline_pop,coords=c("Long","Lat"), crs=EPSG4326)

#Demography Pop
demography_pop <- read_csv("data/genomic_data/offset_pop_9var.csv")
demography_pop <- demography_pop  %>% dplyr::select(Longitude,Latitude)
demography_pop_sf <- st_as_sf(demography_pop,coords=c("Longitude","Latitude"), crs=EPSG4326)


##############################################################################

#Future climate change
#2041-2070 offset
#off_pallet <- c("#2166AC","#67A9CF","#D1E5F0","#f7c1c8","#f21836","#A50F15")
#off_palletas <- c("#004C99","#004C99","#67A9CF","#67A9CF","#D1E5F0","#D1E5F0","#f7c1c8")
#off_pallet45 <- c("#67A9CF","#67A9CF","#D1E5F0","#D1E5F0","#f7c1c8","#f7c1c8")
#off_pallet45 <- c("#D1E5F0","#f7c1c8","#f7c1c8","#f21836","#f21836","#A50F15")
#off_pallet45 <- c("#f7c1c8","#f21836","#f21836","#A50F15","#A50F15","#5c0915")
off_pallet45 <- c("#D1E5F0","#f7c1c8","#f7c1c8","#f21836","#f21836")
#Plot offset SSP245 (RCP 4.5)
tmap_mode("plot")
#tmap_mode("view")
offset45 <- tm_shape(mask_offset_45_grain, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = off_pallet45)+
  #  tm_raster(palette = rev(brewer.pal(6, "RdBu")))+
  #  tm_raster(palette = "Reds")+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(demography_pop_sf)+
  tm_dots(size=0.5,shape=1,col="black",border.lwd = 2.5)+
  #  tm_shape(baseline_pop_sf)+
  #  tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.62, 0.48),legend.title.size = 0.001)
offset45
#tmap_save(offset45, filename = "Graphs/offset45_9var.pdf",width=4, height=7)

#off_pallet2 <- c("#2166AC","#67A9CF","#D1E5F0","#f7c1c8","#f21836","#A50F15","#5c0915")
#off_pallet3 <- c("#2166AC","#67A9CF","#D1E5F0","#f7c1c8","#f21836","#A50F15")
#off_palleta <- c("#004C99","#67A9CF","#D1E5F0","#f7c1c8")
#off_pallet85 <- c("#67A9CF","#67A9CF","#D1E5F0","#D1E5F0","#f7c1c8","#f7c1c8","#f21836")
#off_pallet85 <- c("#f7c1c8","#f21836","#f21836","#A50F15","#A50F15","black","black")
#off_pallet85 <- c("#f21836","#A50F15","#A50F15","magenta2","magenta2","black","black")
#off_pallet85 <- c("#f21836","#A50F15","#A50F15","#5c0915","#5c0915","magenta2","magenta2")
off_pallet85 <- c("#f7c1c8","#f7c1c8","#f21836","#f21836","#A50F15","#A50F15")

#Plot offset SSP585 (RCP 8.5)
tmap_mode("plot")
#tmap_mode("view")
offset85 <- tm_shape(mask_offset_85_grain, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = off_pallet85)+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(demography_pop_sf)+
  tm_dots(size=0.5,shape=1,col="black",border.lwd = 2.5)+
  #  tm_shape(baseline_pop_sf)+
  #  tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.62, 0.48),legend.title.size = 0.001)
offset85
#tmap_save(offset85, filename = "Graphs/offset85_9var.pdf",width=4, height=7)


##############################################################################
#2012 to 2015 offset
#off_pallet <- c("#2166AC","#67A9CF","#D1E5F0","#f7c1c8","#f21836","#A50F15")
#off_pallet1215 <- c("#67A9CF","#D1E5F0","#f7c1c8","#f21836","#A50F15")
#off_pallet1215 <- c("#67A9CF","#67A9CF","#D1E5F0","#D1E5F0","#f7c1c8","#f7c1c8","#f21836")
off_pallet1215 <- c("#67A9CF","#D1E5F0","#f7c1c8","#f21836","#A50F15","#5c0915")
#Plot offset 
tmap_mode("plot")
#tmap_mode("view")
offset_1215 <- tm_shape(mask_offset_1215, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = off_pallet1215)+
  #  tm_raster(palette = rev(brewer.pal(6, "RdBu")))+
  #  tm_raster(palette = "Reds")+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(demography_pop_sf)+
  tm_dots(size=0.5,shape=1,col="black",border.lwd = 2.5)+
  #  tm_shape(baseline_pop_sf)+
  #  tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.62, 0.48),legend.title.size = 0.001)
offset_1215
#tmap_save(offset_1215, filename = "Graphs/offset1215_9var.pdf",width=4, height=7)


##############################################################################
#Clim diff
off_pallet <- c("#2166AC","#67A9CF","#D1E5F0","#f7c1c8","#f21836","#A50F15")
env_pallet <- c("#67A9CF","#D1E5F0","#f7c1c8","#f21836","#A50F15")
off_pallet1215 <- c("#D1E5F0","#D1E5F0","#f7c1c8","#f7c1c8","#f21836","#f21836","#A50F15")
#Plot offset 
tmap_mode("plot")
#tmap_mode("view")
env_1215 <- tm_shape(clim_diff_1215, bbox=st_bbox(calo)) + #legal boundires
  #tm_raster()+
  tm_raster(palette = env_pallet)+
  #  tm_raster(palette = rev(brewer.pal(6, "RdBu")))+
  #  tm_raster(palette = "Reds")+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(demography_pop_sf)+
  tm_dots(size=0.5,shape=1,col="black",border.lwd = 2.5)+
  #  tm_shape(baseline_pop_sf)+
  #  tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.62, 0.48),legend.title.size = 0.001)
env_1215
#tmap_save(env_1215, filename = "Graphs/clim_diff_1215_9var.pdf",width=4, height=7)







