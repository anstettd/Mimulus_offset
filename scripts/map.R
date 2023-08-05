##################################################################################
## Make population mapes of baseline and timeseires
## Author Daniel Anstett
## 
## 
## Last Modified Jan 20, 2023
###################################################################################

#Library install and import
library(tidyverse) # for data manipulation
library(raster) # for handeling rasters
library(tidyverse)
library(cowplot)
library(sf)
library(tmap)
library(rnaturalearth)
library(rnaturalearthdata)
#library(devtools)
#devtools::install_github("ropensci/rnaturalearthhires")
library(rnaturalearthhires)
library(rgeos)
library(RColorBrewer)
##############################################################################
#Define CRS
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS

#Baseline & Timeseries
all_pop <- read_csv("data/genomic_data/Baseline_Timeseries_pops_final2.csv")%>% filter(Paper_ID<56)
all_pop_sf <- st_as_sf(all_pop,coords=c("Long","Lat"), crs=EPSG4326)
timeseries_pop <- all_pop %>% filter(Paper_ID<13) %>% dplyr::select(Long,Lat)
timeseries_pop_sf <- st_as_sf(timeseries_pop,coords=c("Long","Lat"), crs=EPSG4326)
baseline_pop <- all_pop %>% filter(Paper_ID>12) %>% dplyr::select(Long,Lat)
baseline_pop_sf <- st_as_sf(baseline_pop,coords=c("Long","Lat"), crs=EPSG4326)

#Demography
demo_pop <- read_csv("data/genomic_data/offset_pop_beagle.csv")
demo_pop_sf <- st_as_sf(demo_pop ,coords=c("Longitude","Latitude"), crs=EPSG4326)

#Admixture
admix <- read_csv("data/genomic_data/Pop_admixture.csv") %>% 
  filter(!is.na(Latitude),
         !is.na(Longitude))
admix_sf <- st_as_sf(admix ,coords=c("Longitude","Latitude"), crs=EPSG4326)

# California & Oregon Map Setup
states<-ne_states(country=c("canada","united states of america"),returnclass= "sf")
calo <- states %>%
  filter(name_en=="Oregon" | name_en=="California") #| name_en=="Nevada")

#Baseline + Timeseries
tmap_mode("plot")
#tmap_mode("view")
mim <-
  tm_shape(calo)+
  tm_borders()+
  tm_shape(baseline_pop_sf)+
  tm_dots(size=0.4,shape=21,col="blue",border.col ="black")+
  tm_shape(timeseries_pop_sf)+
  tm_dots(size=0.6,shape=21,col="red",border.col ="black")+
  tm_layout(legend.position = c(1.03, 0.73),legend.title.size = 0.001,frame.lwd = NA)
mim
tmap_save(mim, filename = "Graphs/Maps/base_time.png",width=4, height=7)


#Plot Baseline
tmap_mode("plot")
#tmap_mode("view")
mim_base <-
  tm_shape(calo)+
  tm_borders()+
  tm_shape(all_pop_sf)+
  tm_dots(size=0.6,shape=21,col="red2",border.col ="black")+
  tm_layout(legend.position = c(1.03, 0.73),legend.title.size = 0.001,frame.lwd = NA)
mim_base
tmap_save(mim_base, filename = "Graphs/Maps/baseline.png",width=4, height=7)



#Plot Timeseries
tmap_mode("plot")
#tmap_mode("view")
mim_time <-
  tm_shape(calo)+
  tm_borders()+
  tm_shape(timeseries_pop_sf)+
  tm_dots(size=0.9,shape=21,col="blue",border.col ="black")+
  tm_layout(legend.position = c(1.03, 0.73),legend.title.size = 0.001,frame.lwd = NA)
mim_time
tmap_save(mim_time, filename = "Graphs/Maps/timeseries.png",width=4, height=7)


#Plot Timeseries
tmap_mode("plot")
#tmap_mode("view")
mim_time <-
  tm_shape(calo)+
  tm_borders()+
  tm_shape(demo_pop_sf)+
  tm_dots(size=0.9,shape=21,col="magenta2",border.col ="black")+
  tm_layout(legend.position = c(1.03, 0.73),legend.title.size = 0.001,frame.lwd = NA)
mim_time
tmap_save(mim_time, filename = "Graphs/Maps/demography.png",width=4, height=7)


#Plot Admixture
west_coast <- map_data("state", c("oregon", "california","nevada")) 

admix_map <- ggplot(west_coast, aes(x=long, y=lat, group=group)) +
  xlab("Longitude") +
  ylab("Latitude") +
  geom_polygon(fill = "white", colour = "grey50") + 
  geom_scatterpie(aes(x = Longitude, y = Latitude, group=Site), data = admix, cols=colnames(admix[,c(3:6)]), size = 0.1)+
  scale_fill_manual(values=c("red2","green","yellow","deepskyblue")) +
  theme_classic()

ggsave(admix_map, filename = "Graphs/Maps/admixture.png",width=6, height=7)

