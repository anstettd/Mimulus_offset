##################################################################################
## Gradient forest for full SNP dataset
## Author Daniel Anstett
## 
## Modified from Keller & Fitzpatric 2015
## Last Modified September 19, 2022
###################################################################################

#Library install and import

# Clear environment
rm(list = ls())

# Install packages
#install.packages("extendedForest", repos="http://R-Forge.R-project.org")
#install.packages("gradientForest", repos="http://R-Forge.R-project.org")

library(randomForest)
library(gradientForest)
library(tidyverse)
library(dismo) # for biovars function
library(raster) # for raster grids
library(rgdal) # for transforming projections
library(sf)
library(vegan)
library(gdm)

############################################################################################################

#Function 1 : go from PCA to raster

# builds RGB raster from transformed environment
# snpPreds = dataframe of transformed variables from gf or gdm model
# rast = a raster mask to which RGB values are to be mapped
# cellNums = cell IDs to which RGB values should be assigned
pcaToRaster <- function(snpPreds, rast, mapCells){
  require(raster)
  
  pca <- prcomp(snpPreds, center=TRUE, scale.=FALSE)
  
  ##assigns to colors, edit as needed to maximize color contrast, etc.
  a1 <- pca$x[,1]
  a2 <- pca$x[,2]
  a3 <- pca$x[,3]
  
  r <- a1+a2 
  g <- -a2 
  b <- a3+a2-a1
  
  ##scales colors
  scalR <- (r-min(r))/(max(r)-min(r))*255
  scalG <- (g-min(g))/(max(g)-min(g))*255
  scalB <- (b-min(b))/(max(b)-min(b))*255
  
  ##assigns color to raster
  rast1 <- rast2 <- rast3 <- rast
  rast1[mapCells] <- scalR
  rast2[mapCells] <- scalG
  rast3[mapCells] <- scalB
  ##stacks color rasters
  outRast <- stack(rast1, rast2, rast3)
  return(outRast)
}


#Function 2: Map difference between spatial genetic predictions

# predMap1 = dataframe of transformed variables from gf or gdm model for first set of SNPs
# predMap2 = dataframe of transformed variables from gf or gdm model for second set of SNPs
# rast = a raster mask to which Procrustes residuals are to be mapped
# mapCells = cell IDs to which Procrustes residuals values should be assigned
RGBdiffMap <- function(predMap1, predMap2, rast, mapCells){
  require(vegan)
  PCA1 <- prcomp(predMap1, center=TRUE, scale.=FALSE)
  PCA2 <- prcomp(predMap2, center=TRUE, scale.=FALSE)
  diffProcrust <- procrustes(PCA1, PCA2, scale=TRUE, symmetrical=FALSE)
  residMap <- residuals(diffProcrust)
  rast[mapCells] <- residMap
  return(list(max(residMap), rast))
}

#Notes from Fitzpatrick & Keller
# OK, on to mapping. Script assumes:
# (1) a dataframe named env_trns containing extracted raster data (w/ cell IDs)
# and env. variables used in the models & with columns as follows: cell, bio1, bio2, etc.
#
# (2) a raster mask of the study region to which the RGB data will be written



#Function 3: Generating a raster stack

#From Tongli Wang

# x is the wrodkig directory
# varList: varList <- c('mat','map','td')
# rType='tif' or 'grid'
# vConvert=T: 1/10 for mat,mwmt, ...
rasterStack <- function(x,varList,rType='tif',vConvert=T){
  library(raster)
  for(var in varList){
    if(rType=='tif'){inF <- paste0(x,'/',var,'.tif')}
    if(rType=='grid'){inF <- paste0(x,'/',var)}
    r <- raster(inF)
    if(vConvert==T){
      if(grepl('mat|MAT|mwmt|MWMT|mcmt|MCMT|td|TD|emt|EMT|ext|EXT             
                |ahm|AHM|shm|SHM|mar|MAR|ta|Ta|tm|Tm',names(r))){r <- r/10}
    }
    if(var==varList[1]){stk=r} else{stk=stack(stk,r)}
  }
  return(stk)
}

#Examples
#varList <- c('mat','map','td','tmin_sp','ahm','tmax_wt')
#wd <- 'G:/ClimateXX_out/BC2/BC800v600/Normal_1961_1990MSY'
#stk <- rasterStack(wd,varList,rType='tif',vConvert=T);stk

############################################################################################################
############################################################################################################

## Import SNP data and arrange for gradient forest
#Import SNP Data & and reformat
snp_clim_bf20NA <- read_csv("data/genomic_data/snp_clim_peakbf10_noNA.csv") #pop data
test_snp <- snp_clim_bf20NA %>% dplyr::select(-Site_Name, -Paper_ID, -Latitude, -Longitude, -Elevation, -MAT, -MAP, -CMD,
                                       -PAS, -EXT, -Tave_wt, -Tave_sm, -PPT_wt, -PPT_sm)

## Generate specific dataframes for GF model
env_site <- snp_clim_bf20NA %>% dplyr::select(MAT,MAP,PAS,EXT,CMD,Tave_wt,Tave_sm,PPT_wt,PPT_sm)


#Convert to data frame
test_snp<-as.data.frame(test_snp)
env_site<-as.data.frame(env_site)


#Merge and make labels
df_in_1<-data.frame(env_site, test_snp)
resp<-colnames(test_snp)
pred<-colnames(env_site)

############################################################################################################

## Range wide polygon
# Import M.cardinalis ensamble range extent as sf polygon
#c_range <- st_read("SDM/Output/c_range_2.shp")
c_range <- st_read("data/genomic_data/Shape/c_range50.shp") 
c_range <- st_transform(c_range, crs = 4326) # reproject to WGS 1984 (EPSG 4326)


## Raster import and manipulation
#Import 1981-2010 raster data for West NA & and stack them
wd <- "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_8110"
vlist <- c("MAT","MAP","PAS","EXT","CMD","Tave_wt","Tave_sm","PPT_wt","PPT_sm")
stk <- rasterStack(wd,vlist,rType='tif',vConvert=F)

#Reproject to WGS 1984 (EPSG4326)
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS
stk <- projectRaster(stk, crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)


#Check CRS of Polygon and raster
crs(c_range) 
crs(stk)


#Clip raster using range-extent polygon
stk.clip <- raster::crop(stk, extent(c_range))
stk.mask <- mask(stk.clip, c_range)

#Extract point from raster stack
stk.df <- data.frame(rasterToPoints(stk.mask))
stk.df <- na.omit(stk.df)

#Convert xy coordinates into cell ID
stk.df.cell<-cellFromXY(stk.mask, cbind(stk.df$x, stk.df$y))

############################################################################################################
##Import 2012-2015 raster average

#######################
#Import 2012-2015 raster data for West NA & and stack them
wd <- "C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_1215"
vlist <- c("MAT","MAP","PAS","EXT","CMD","Tave_wt","Tave_sm","PPT_wt","PPT_sm")
stk_2012 <- rasterStack(wd,vlist,rType='tif',vConvert=F)

#Reproject to WGS 1984 (EPSG4326)
stk_2012 <- projectRaster(stk_2012, crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)
crs(stk_2012)

#Clip raster using range-extent polygon
stk_2012.mask <- raster::crop(stk_2012, extent(c_range))
stk_2012.mask <- mask(stk_2012.mask, c_range)

#Extract point from raster stack
stk_2012.df <- data.frame(rasterToPoints(stk_2012.mask))
stk_2012.df <- na.omit(stk_2012.df)
colnames(stk_2012.df)[3]<-"MAT"
colnames(stk_2012.df)[4]<-"MAP"
colnames(stk_2012.df)[5]<-"PAS"
colnames(stk_2012.df)[6]<-"EXT"
colnames(stk_2012.df)[7]<-"CMD"
colnames(stk_2012.df)[8]<-"Tave_wt"
colnames(stk_2012.df)[9]<-"Tave_sm"
colnames(stk_2012.df)[10]<-"PPT_wt"
colnames(stk_2012.df)[11]<-"PPT_sm"
#Convert xy coordinates into cell ID
stk_2012.df.cell<-cellFromXY(stk_2012.mask, cbind(stk_2012.df$x, stk_2012.df$y))





#Get mask for RGB
MAT.clip <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/Raster_updated/Year_8110/MAT.tif")
MAT.clip <- projectRaster(MAT.clip, crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)
rbg_mask <- raster::crop(MAT.clip, extent(c_range))
rbg_mask <- mask(rbg_mask, c_range)
rbg_mask <- rbg_mask * 0
mask_offset_2012 <- rbg_mask
mask_offset_2012_dist <- rbg_mask
rbg_2012 <- rbg_mask




############################################################################################################
############################################################################################################
#Train Gradient Forest Model

maxLevel <- log2(0.368*nrow(df_in_1)/2) #account for correlations, see ?gradientForest 

# Gradient Forest Model
gf <- gradientForest(df_in_1,predictor.vars = pred, 
                     response.vars = resp,ntree = 500, 
                     maxLevel = maxLevel, trace=T, corr.threshold = 0.5)
plot(gf) #Importance Plot


# Conpact version for larger datasets
#gf <- gradientForest(df_in_1,predictor.vars = pred, 
#                     response.vars = resp,ntree = 500, 
#                     transform = NULL, compact = T, nbin = 201 , corr.threshold = 0.5)


# transform env using gf models, see ?predict.gradientForest
predBF20 <- predict(gf, stk.df[,-1:-2]) # remove cell column before transforming

############################################################################################################
#Part 2, use gradient forest model on rasters

# Calculate and map "genetic offset" under climate change ----------------------
# Script assumes:
# (1) a dataframe of transformed env. variables for CURRENT climate 
# (e.g., predGI5 from above).
#
# (2) a dataframe named env_trns_future containing extracted raster data of 
# env. variables for FUTURE a climate scenario, same structure as env_trns

#2011 to 2016 climate variables

# first transform FUTURE env. variables
projBF20_2012 <- predict(gf, stk_2012.df[,-1:-2])


# calculate euclidean distance between current and future genetic spaces  
#offset_BF20_2012 <- sqrt((projBF20_2012[,1]-predBF20[,1])^2+(projBF20_2012[,2]-predBF20[,2])^2
#                         +(projBF20_2012[,3]-predBF20[,3])^2)

offset_BF20_2012 <- sqrt((projBF20_2012[,1]-predBF20[,1])^2
                        +(projBF20_2012[,2]-predBF20[,2])^2
                        +(projBF20_2012[,3]-predBF20[,3])^2
                        +(projBF20_2012[,4]-predBF20[,4])^2
                        +(projBF20_2012[,5]-predBF20[,5])^2
                        +(projBF20_2012[,6]-predBF20[,6])^2
                        +(projBF20_2012[,7]-predBF20[,7])^2
                        +(projBF20_2012[,8]-predBF20[,8])^2
                        +(projBF20_2012[,9]-predBF20[,9])^2)


# assign values to raster - can be tricky if current/future climate
# rasters are not identical in terms of # cells, extent, etc.

mask_offset_2012[stk_2012.df.cell] <- offset_BF20_2012
plot(mask_offset_2012)

#writeRaster(mask_offset_2012,"data/genomic_data/offset_1215.tif", format="GTiff", overwrite=TRUE)













#Get Climate Distance

#Make PCAs
past_pca <- prcomp(stk.df[3:5])
             
#Predict using past PCA
pred_past_env <- predict(past_pca, stk.df[,-1:-2]) # remove cell column before transforming
pred_1215_env <- predict(past_pca, stk_2012.df[,-1:-2])


clim_distance_1215 <- sqrt((pred_past_env[,1]-pred_1215_env[,1])^2+(pred_past_env[,2]-pred_1215_env[,2])^2
                           +(pred_past_env[,3]-pred_1215_env[,3])^2)


# calculate euclidean distance between current and future genetic spaces  
#clim_distance_1215 <- sqrt((stk.df[,3]-stk_2012.df[,3])^2+(stk.df[,4]-stk_2012.df[,4])^2
#                           +(stk.df[,5]-stk_2012.df[,5])^2)

mask_offset_2012_dist[stk_2012.df.cell] <- clim_distance_1215
plot(mask_offset_2012_dist)

#writeRaster(mask_offset_2012_dist,"data/genomic_data/clim_distance.tif", format="GTiff", overwrite=TRUE)






























