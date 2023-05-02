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

#Define CRS
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS

# California & Oregon Map Setup
states<-ne_states(country=c("canada","united states of america"),returnclass= "sf")
calo <- states %>%
  filter(name_en=="Oregon" | name_en=="California" | name_en=="Nevada")

#Demography Pop
demography_pop <- read_csv("data/genomic_data/offset_pop_9var.csv")
demography_pop <- demography_pop  %>% dplyr::select(Longitude,Latitude)
demography_pop_sf <- st_as_sf(demography_pop,coords=c("Longitude","Latitude"), crs=EPSG4326)

#Import sensitivity rasters all for 1215
bf30 <- raster("data/genomic_data/offset_bf30.tif") #BF>30 only
env9 <- raster("data/genomic_data/offset_env9.tif") #new WZA
oldWZA <- raster("data/genomic_data/offset_oldWZA_snp_set.tif") #WZA snps only before update
old_snp_set <- raster("data/genomic_data/offset_old_snp_set.tif") #Original WZA snps + BF30
new_snp_set <- raster("data/genomic_data/offset_1215.tif") #New WZA snps + BF30


##############################################################################
#BF30
#off_pallet <- c("#2166AC","#67A9CF","#D1E5F0","#f7c1c8","#f21836","#A50F15")
#off_pallet1215 <- c("#67A9CF","#D1E5F0","#f7c1c8","#f21836","#A50F15")
#off_pallet1215 <- c("#67A9CF","#67A9CF","#D1E5F0","#D1E5F0","#f7c1c8","#f7c1c8","#f21836")
off_pallet1215 <- c("#67A9CF","#D1E5F0","#f7c1c8","#f21836","#A50F15","#5c0915")
#Plot offset 
tmap_mode("plot")
#tmap_mode("view")
offset_bf30 <- tm_shape(bf30, bbox=st_bbox(calo)) + #legal boundires
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
offset_bf30
tmap_save(offset_bf30, filename = "Graphs/Sensetivity/offset1215_bf30.pdf",width=4, height=7)


#env9
off_pallet1215 <- c("#67A9CF","#D1E5F0","#f7c1c8","#f21836","#A50F15","#5c0915")
#Plot offset 
tmap_mode("plot")
#tmap_mode("view")
offset_env9 <- tm_shape(env9, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = off_pallet1215)+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(demography_pop_sf)+
  tm_dots(size=0.5,shape=1,col="black",border.lwd = 2.5)+
  tm_layout(legend.position = c(0.62, 0.48),legend.title.size = 0.001)
offset_env9
tmap_save(offset_env9, filename = "Graphs/Sensetivity/offset1215_9var.pdf",width=4, height=7)


#oldWZA
off_pallet1215 <- c("#67A9CF","#D1E5F0","#f7c1c8","#f21836","#A50F15","#5c0915")
#Plot offset 
tmap_mode("plot")
#tmap_mode("view")
offset_oldWZA <- tm_shape(oldWZA, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = off_pallet1215)+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(demography_pop_sf)+
  tm_dots(size=0.5,shape=1,col="black",border.lwd = 2.5)+
  tm_layout(legend.position = c(0.62, 0.48),legend.title.size = 0.001)
offset_oldWZA
tmap_save(offset_oldWZA, filename = "Graphs/Sensetivity/offset1215_oldWZA.pdf",width=4, height=7)


#old_snp_set
off_pallet1215 <- c("#67A9CF","#D1E5F0","#f7c1c8","#f21836","#A50F15","#5c0915")
#Plot offset 
tmap_mode("plot")
#tmap_mode("view")
offset_old_snp_set <- tm_shape(old_snp_set, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = off_pallet1215)+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(demography_pop_sf)+
  tm_dots(size=0.5,shape=1,col="black",border.lwd = 2.5)+
  tm_layout(legend.position = c(0.62, 0.48),legend.title.size = 0.001)
offset_old_snp_set
tmap_save(offset_old_snp_set, filename = "Graphs/Sensetivity/offset1215_old_snp_set.pdf",width=4, height=7)


#new_snp_set
off_pallet1215 <- c("#67A9CF","#D1E5F0","#f7c1c8","#f21836","#A50F15","#5c0915")
#Plot offset 
tmap_mode("plot")
#tmap_mode("view")
offset_new_snp_set <- tm_shape(new_snp_set, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = off_pallet1215)+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(demography_pop_sf)+
  tm_dots(size=0.5,shape=1,col="black",border.lwd = 2.5)+
  tm_layout(legend.position = c(0.62, 0.48),legend.title.size = 0.001)
offset_new_snp_set
tmap_save(offset_new_snp_set, filename = "Graphs/Sensetivity/offset1215_new_snp_set.pdf",width=4, height=7)





