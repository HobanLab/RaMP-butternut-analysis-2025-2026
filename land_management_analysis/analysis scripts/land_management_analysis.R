#Big Thanks To Rebecca!!!!
#Credit her Github https://github.com/HobanLab/QUBR_GenGeoEcoDemoCorr/blob/main/analyses/Data_Processing_Script.R and https://github.com/HobanLab/QUBR_GenGeoEcoDemoCorr/blob/main/analyses/hypothesis_1_clumping_CLEANED.R

#### Loading libraries and relevant data ####

library(tidyverse) 
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc)
library(ggplot2)
library(PMCMRplus) # for Dunn test
library(geomtextpath) # for PCA graphing
library(spatstat) # to run the Ripley's K function: Kest
library(stars) # for sf_rasterize function
library(raster) #to use crop
library(starsExtra) #to use dist_to_nearest
library(geostatsp) # To successfully use as.im
library(tmaptools)
library(mgcv) #to use GAM function 
library(plotly) #to 3d plot variables
library(MuMIn) #to be able to use dredge
library(visreg) # to be able to plot Aspect/categorical variables with GAM
library(spatialEco) # for heat load index function
library(geosphere) # for finding the distance of populations to the trails

setwd("~/land_management_analysis_cfd")
data_ILM <- read.csv("data_ILM_clean_processed.csv")

data_ILM_sf <- st_as_sf(data_ILM, coords = c("gps_w", "gps_n"), crs = 4326)
plot(data_ILM_sf)

data_ILM_sf_coords <- st_coordinates(data_ILM_sf) #creates a dataframe with separate x and y columns
 data_ILM_sf_coords <- data_ILM_sf %>%
  cbind(data_ILM_sf_coords) #combines the x and y coordinate data frame with the transformed sf dataframe

#can transform health response variables if necessary---refer back to https://github.com/HobanLab/QUBR_GenGeoEcoDemoCorr/blob/main/analyses/Data_Processing_Script.R

#adding the ILM trails into our R space from the ArcGIS SHP file
setwd("~/land_management_analysis_cfd/trails_all_ILM")
trails_ILM <- st_read("BufferedFeatures.shp")
trails_ILM <- trails_ILM$geometry[1]
trails_ILM <- st_transform(trails_ILM, crs = 26918)
plot(trails_ILM)

#creating a buffer around the trails of 20 meters---these include the trees not ON the trail, but close to them
trail_buffer_ILM <- st_buffer(trails_ILM, 20)

ggplot() +
  geom_sf(data = trails_ILM) +
  geom_sf(data = trail_buffer_ILM) +
  geom_sf(data = data_ILM_sf)



#this comes from https://github.com/HobanLab/QUBR_GenGeoEcoDemoCorr/blob/main/analyses/hypothesis_1_clumping_CLEANED.R

#finding minimum and maximum lat and long values for visualizing the tree locations with *1.002 wiggle room
min_ILM_long <- min(data_ILM_sf_coords$X)*1.0002
max_ILM_long <- max(data_ILM_sf_coords$X) - (max(data_ILM_sf_coords$X) *.0002)
min_ILM_lat <- min(data_ILM_sf_coords$Y)*1.02
max_ILM_lat <- max(data_ILM_sf_coords$Y) - (max(data_ILM_sf_coords$Y)*.02)

#transforming to UTM 18N https://spatialreference.org/ref/epsg/26918/
data_ILM_sf_transformed <- st_transform(data_ILM_sf, crs = 26918)
#making the boundary box
data_ILM_sf_box <- data_ILM_sf_transformed %>%
  st_bbox %>%
  st_as_sfc()

win <- as.owin(data_ILM_sf_box) #turning the box into a window...I hope? this is where I think the code breaks down
ppp <- as.ppp(st_coordinates(data_ILM_sf_transformed), W = win)
#apparently every single point is outside the window.  If I use data_ILM_sf, the data contains duplicated points.  Not sure what either of those things mean lol

#some debugging
# get coordinates from sf object
coords_test <- matrix(unlist(data_ILM_sf_transformed$geometry),ncol=2,byrow=T)
# inspect
coords_test
# create ppp object
test_ppp <- ppp(x=coords_test[,1],y=coords_test[,2],window=win,check=T)
#duplicated points

########################################
#FIX THE ABOVE ERRORS BEFORE CONTINUING#
########################################

#finding duplicates
anyDuplicated(data_ILM$number)
anyDuplicated(data_ILM$gps_n & data_ILM$gps_w)
duplicates <- which(duplicated((data_ILM$number)))
#filtering out duplicates
data_ILM_no_duplicates <- data_ILM %>%
  filter(data_ILM$number == data_ILM$number[c(duplicates)])

plot(ppp, pch = 16, cex = 0.5) # plotting the randomized point pattern
K <- Kest(ppp, correction = "Ripley") # Ripley's K function, using the isotropic/Ripley correction
plot(K, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) # plotting Ripley's output

# Buffer River (100 m)
ILM_trail_buffer_wind <- as.owin(trail_buffer_ILM) #turning the buffer into a window
LM_ppp_buffer <- as.ppp(st_coordinates(data_ILM_sf), W = ILM_trail_buffer_wind) #creating the poisson point pattern for lm
plot(LM_ppp_buffer, pch = 16, cex = 0.5) # plotting the randomized point pattern
LM_k_buffer <- Kest(LM_ppp_buffer, correction = "Ripley") # Ripley's K function, using the isotropic/Ripley correction
plot(LM_k_buffer, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE))  # plotting Ripley's output

##########################
#AVERAGE NEAREST NEIGHBOR#
##########################

trails_ILM_trans <- st_transform(trails_ILM, crs = 26918)
trails_ILM_trans <- st_sf(geometry = trails_ILM_trans)

#turning the river polygon into a linestring object and then into a raster, to be able to later calculate the distances
trails_ILM_trans_points <- st_cast(trails_ILM_trans, "MULTILINESTRING") #turning the polyline of the river into a linestring object
trails_ILM_trans_point_raster <- st_rasterize(trails_ILM_trans_points) #creating a raster out of the river linestring object
plot(trails_ILM_trans_point_raster) #plotting the river linestring object

#turning the river buffer polygon into a linestring object and then into a raster to be able to later calculate the distances
trails_ILM_buffer_trans_outline <- st_cast(trail_buffer_ILM, "MULTILINESTRING") #turning the polygon of the river buffer into a linestring object
trails_buffer_point_raster <- st_rasterize(trails_ILM_buffer_trans_outline) #creating a raster of river buffer linestring object
plot(trails_buffer_point_raster) #plotting the river buffer linestring object

trails_buffer_point_raster[is.na(trails_buffer_point_raster[])] <- 0  #making sure the cell that are not the river buffer linestring raster have a 0 value
dist_near_trails_buffer <- dist_to_nearest(trails_buffer_point_raster, trails_ILM_trans_points, progress = T) #creating a raster of the distances of each cell in the buffer raster to the linestring object of the river polygon, this can take a while to run

trail_vect <- project(vect(trails_ILM), rast(dist_near_trails_buffer))
