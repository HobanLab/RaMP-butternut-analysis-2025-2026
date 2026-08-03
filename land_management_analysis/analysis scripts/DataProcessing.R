#Catherine dell'Olio
#Land management analysis
#Data processing script

#First off, a big thanks To Rebecca!!!!
#Credit to her Github https://github.com/HobanLab/QUBR_GenGeoEcoDemoCorr/blob/main/analyses/Data_Processing_Script.R, upon which this is based.  I copied this script over to adapt this for my butternut RaMP work.
#This script (and all scripts in this folder) is essentially Catherine's land management edits to Rebecca's work.  This code uses the framework Rebecca laid out to analyze Quercus brandegii trees along three riverbed sites.  I am adapting it for use of butternuts in two sites that have been cleared for trail building.
#This script processed our health assessment and trail polygon data into spatial objects that we can plot and perform analyses on.

#### Loading libraries and relevant data ####

library(tidyverse) 
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc)
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
library(geosphere) # for finding the distance of populations to the coast
library(tmap) # for making beautiful maps and graphics

#### Loading and processing relevant data ####

setwd("~/land_management_analysis/data")
ILM_fixed_field_data_processed <- read.csv("data_ILM.csv") #renames the csv created from the 2025 health data. This is the SAME data_ILM sheet (just saved as a csv) created from my script DataSplitter.R, which is found in a separate part of this same repository, so feel free to check that out, but I have the csv downloaded here for ease of use.
#I'm essentially just renaming data_ILM.csv here to match Rebecca's original naming schema.

#I was receiving the error that some trees in this dataset had the same gps points (such as 67 and 67a, two trunks growing out of the same base), and that will cause problems with our Poisson Point analysis later on, so I remove those here
anyDuplicated(ILM_fixed_field_data_processed[,c('gps_n','gps_w')])
duplicates <- which(duplicated(ILM_fixed_field_data_processed[,c('gps_n','gps_w')]) == FALSE)
ILM_fixed_field_data_processed <- ILM_fixed_field_data_processed[c(duplicates),]
anyDuplicated(ILM_fixed_field_data_processed[,c('gps_n','gps_w')]) #checking that all the duplicates are gone!

# creating the point shapefiles of the tree locations for each population in UTM 18 N

#creating a point shapefile of all points with lat lon coordinates and other attributes in WGS 1984
#sf objects are dataframes with rows representing simple features with attributes and a simple feature geometry list-column (sfc)
ILM_fixed_field_data_processed_sf <- st_as_sf(ILM_fixed_field_data_processed, 
                                          coords = c("gps_w", "gps_n"), crs = 4326)

#creating a transformed point shapefile with UTM 18 N an equal area projection
ILM_fixed_field_data_processed_sf_transformed <- st_transform(ILM_fixed_field_data_processed_sf, crs = 26918)

#storing point shapefiles for the trees by population
ILM_fixed_field_data_processed_sf <- ILM_fixed_field_data_processed_sf_transformed %>%
  st_as_sf()

#create dataframe with X and Y UTM coordinates
ILM_fixed_field_data_processed_sf_trans_coords <- st_coordinates(ILM_fixed_field_data_processed_sf_transformed) #creates a dataframe with separate x and y coordinate columns from the UTM 18N transformation
ILM_fixed_field_data_processed_sf_trans_coordinates <- ILM_fixed_field_data_processed_sf_transformed %>%
  cbind(ILM_fixed_field_data_processed_sf_trans_coords) #combines the x and y coordinate data frame with the transformed sf dataframe

# Checking that the site is correct!
ILM_fixed_field_data_processed <- ILM_fixed_field_data_processed_sf_trans_coordinates %>%
  filter(site == "ILM") 

#### Loading in the trail lines and ILM polygon, originally from ArcGIS ####
setwd("~/land_management_analysis/data/trails_all_ILM")
trails_ILM <- st_read("BufferedFeatures.shp")
trails_ILM <- trails_ILM$geometry[1]
plot(trails_ILM)

#changing the coordinate reference system of the river polygons to be the same equal area projection (UTM 18N), uses meters as distance measurement
trails_ILM_trans <- st_transform(trails_ILM, crs = 26918)

#ensuring the river outlines are shapefiles for the distance measurements
trails_ILM_trans <- st_sf(geometry = trails_ILM_trans)

#turning the ILM into a shapefile, to be able to visualize the point locations
setwd("~/land_management_analysis/data/ILM_boundaries")
ILM_polygon <- read_sf("ILM_boundaries.shp")
ILM_polygon <- st_as_sf(ILM_polygon)
plot(ILM_polygon)

#creating buffers around the trails, which I can use to confirm the trees closest to the trails
trails_buffer_ILM <- st_buffer(trails_ILM, 10) #10 m buffer

#creating bounding boxes for the space around the trails (for Poisson point modelling later)
ILM_box <- st_bbox(trails_ILM_trans)
#I also use the whole of ILM as a bounding box later on, when creating distance-to-trail rasters

#let's plot what we have so far!
ggplot()+  #plotting the trail shapefile, the buffer, and the tree points
  geom_sf(data = trails_buffer_ILM)+
  geom_sf(data = trails_ILM_trans)+ #is this working???????
  geom_sf(data = ILM_fixed_field_data_processed_sf)

#This will make a much prettier version of the map above, because it is specifically designed to work with maps rather than plots
tm_shape(ILM_polygon) + 
  tm_polygons(fill = "bisque", col = "bisque") +
  tm_shape(trails_ILM) +
  tm_polygons(fill = "burlywood4", col = "burlywood4") +
  tm_shape(ILM_fixed_field_data_processed_sf) +
  tm_dots(fill = "green4") +
  tm_title("ILM butternuts (green) and ILM trails (brown)") +
  tm_grid() +
  tm_graticules()

#creating a boundary box for LM with the UTM 18 N min and max lat lon values and then turning it into a simple feature geometry
ILM_fixed_field_data_processed_box <- ILM_fixed_field_data_processed_sf_transformed %>%
  filter(site == "ILM") %>%
  st_bbox %>%
  st_as_sfc()

#histograms
ggplot(ILM_fixed_field_data_processed) + #DBH
  geom_histogram(aes(x = dbh_cm))+
  xlab("Diameter at breast height")+
  ylab("Frequency")

#---------------------------------
#Distance-to-trail raster building
#THIS TAKES A LONG TIME TO RUN!!!!

#Creating distance to trail calculations!
#turning the trail polygons into multilinestring objects and then into a raster, to be able to later calculate the distances
trails_ILM_trans_points <- st_cast(trails_ILM_trans, "MULTILINESTRING") #turning the polyline of the trails into a multilinestring object
trails_ILM_trans_point_raster <- st_rasterize(trails_ILM_trans_points) #creating a raster out of the trail linestring object
plot(trails_ILM_trans_point_raster, breaks = "equal") #plotting the trail multilinestring object

#doing the same with the buffers around the trails, created just above
trails_buffer_ILM_trans <- st_transform(trails_buffer_ILM, crs = 26918)
trails_buffer_ILM_trans <- st_sf(geometry = trails_buffer_ILM_trans)
trails_ILM_buffer_trans_outline <- st_cast(trails_buffer_ILM_trans, "MULTILINESTRING")
trails_buffer_ILM_point_raster <- st_rasterize(trails_ILM_buffer_trans_outline)
plot(trails_buffer_ILM_point_raster)

#generating a distance to trail raster with the distances of each cell in the buffer raster from the trail edge points, and the trail raster cells are set to a distance of 0 m
trails_buffer_ILM_point_raster[is.na(trails_buffer_ILM_point_raster[])] <- 0  #making sure the cell that are not the trail buffer multilinestring raster have a 0 value
load("dist_near_trails_buffer_ILM.rData")
#if the line above isn't working / you want to make the raster over again
#dist_near_trails_buffer_ILM <- dist_to_nearest(trails_buffer_ILM_point_raster, trails_ILM_trans_points, progress = T) #creating a raster of the distances of each cell in the buffer raster to the linestring object of the trail polygons, this can take a while to run
#save(dist_near_trails_buffer_ILM, file="dist_near_trails_buffer_ILM.rData")
plot(dist_near_trails_buffer_ILM) #plotting the distance to river raster

#trying instead of using the buffer to define how big the raster can be, use the whole box ILM_fixed_field_data_processed_box
trails_fixed_box_ILM_trans <- st_transform(ILM_fixed_field_data_processed_box, crs = 26918)
trails_fixed_box_ILM_trans <- st_sf(geometry = trails_fixed_box_ILM_trans)
trails_ILM_fixed_box_trans_outline <- st_cast(trails_fixed_box_ILM_trans, "MULTILINESTRING")
trails_fixed_box_ILM_point_raster <- st_rasterize(trails_ILM_fixed_box_trans_outline)
plot(trails_fixed_box_ILM_point_raster)

trails_fixed_box_ILM_point_raster[is.na(trails_fixed_box_ILM_point_raster[])] <- 0
load("dist_near_trails_fixed_box_ILM.rData")
#if the line above isn't working / you want to make the raster over again
#dist_near_trails_fixed_box_ILM <- dist_to_nearest(trails_fixed_box_ILM_point_raster, trails_ILM_trans_points, progress = T)
#save(dist_near_trails_fixed_box_ILM, file="dist_near_trails_fixed_box_ILM.rData")
plot(dist_near_trails_fixed_box_ILM)

## Making it so the cells in the distance raster within or overlapping with the river raster are assigned 1 
#Assigning points within and overlapping with the river to be "true"
ILM_points_intersects_trails <- st_intersects(ILM_fixed_field_data_processed, trails_ILM_trans, sparse = F) #creating a list of true or falses for whether points intersect the river shapefiles
ILM_fixed_field_data_processed_intersects_trails <- cbind(ILM_fixed_field_data_processed, ILM_points_intersects_trails) #binding the list of true or falses with the point data
#printing the river polygon and the tree points, colored by whether are or aren't within or overlapping with the river
ggplot()+
  geom_sf(data=trails_ILM_trans)+
  geom_sf(data=ILM_fixed_field_data_processed)+
  geom_sf(data=ILM_fixed_field_data_processed_intersects_trails, aes(color = ILM_points_intersects_trails))


## Extracting distance to trail for each tree using the distance to trail raster

ILM_distance_data_pts <- st_extract(dist_near_trails_buffer_ILM, ILM_fixed_field_data_processed) #extracting distance to trail for each tree
ILM_fixed_field_data_processed_distance  <- cbind(ILM_fixed_field_data_processed, ILM_distance_data_pts) #binding the distance to trail data for each point to the ILM point dataframe

## Assigning all points within/overlapping trail to distances of 0

ILM_fixed_field_data_processed_distance <- ILM_fixed_field_data_processed_distance %>% 
  mutate(d = case_when((ILM_fixed_field_data_processed_intersects_trails$ILM_points_intersects_trails == T) ~ 0,  #assigns 0 to points within trail
                       (ILM_fixed_field_data_processed_intersects_trails$ILM_points_intersects_trails == F) ~ d)) #to points outside of trail, it leaves the original distance value

#----
#CPVT
#----

#### Loading and processing relevant data ####

setwd("~/land_management_analysis/data")
CPVT_fixed_field_data_processed <- read.csv("data_CPVT.csv") #renames the csv created from the 2025 health data. This is the SAME data_CPVT sheet (just saved as a csv) created from my script DataSplitter.R, which is found in a separate part of this same repository, so feel free to check that out, but I have the csv downloaded here for ease of use.
#I'm essentially just renaming data_CPVT.csv here to match Rebecca's original naming schema.

#I was receiving the error that some trees in this dataset had the same gps points (such as 67 and 67a, two trunks growing out of the same base), and that will cause problems with our Poisson Point analysis later on, so I remove those here
anyDuplicated(CPVT_fixed_field_data_processed[,c('gps_n','gps_w')])
duplicates <- which(duplicated(CPVT_fixed_field_data_processed[,c('gps_n','gps_w')]) == FALSE)
CPVT_fixed_field_data_processed <- CPVT_fixed_field_data_processed[c(duplicates),]
anyDuplicated(CPVT_fixed_field_data_processed[,c('gps_n','gps_w')]) #checking that all the duplicates are gone!

# creating the point shapefiles of the tree locations for each population in UTM 18 N

#creating a point shapefile of all points with lat lon coordinates and other attributes in WGS 1984
#sf objects are dataframes with rows representing simple features with attributes and a simple feature geometry list-column (sfc)
CPVT_fixed_field_data_processed_sf <- st_as_sf(CPVT_fixed_field_data_processed, 
                                              coords = c("gps_w", "gps_n"), crs = 4326)

#creating a transformed point shapefile with UTM 18 N an equal area projection
CPVT_fixed_field_data_processed_sf_transformed <- st_transform(CPVT_fixed_field_data_processed_sf, crs = 26918)

#storing point shapefiles for the trees by population
CPVT_fixed_field_data_processed_sf <- CPVT_fixed_field_data_processed_sf_transformed %>%
  st_as_sf()

#create dataframe with X and Y UTM coordinates
CPVT_fixed_field_data_processed_sf_trans_coords <- st_coordinates(CPVT_fixed_field_data_processed_sf_transformed) #creates a dataframe with separate x and y coordinate columns from the UTM 18N transformation
CPVT_fixed_field_data_processed_sf_trans_coordinates <- CPVT_fixed_field_data_processed_sf_transformed %>%
  cbind(CPVT_fixed_field_data_processed_sf_trans_coords) #combines the x and y coordinate data frame with the transformed sf dataframe

# Checking that the site is correct!
CPVT_fixed_field_data_processed <- CPVT_fixed_field_data_processed_sf_trans_coordinates %>%
  filter(site == "CPVT") 

#### Loading in the trail lines and CPVT polygon, originally from ArcGIS ####
setwd("~/land_management_analysis/data/trails_all_CPVT")
trails_CPVT <- st_read("BufferedFeatures.shp")
trails_CPVT <- trails_CPVT$geometry[1]
plot(trails_CPVT)

#changing the coordinate reference system of the river polygons to be the same equal area projection (UTM 18N), uses meters as distance measurement
trails_CPVT_trans <- st_transform(trails_CPVT, crs = 26918)

#ensuring the river outlines are shapefiles for the distance measurements
trails_CPVT_trans <- st_sf(geometry = trails_CPVT_trans)

#turning the CPVT into a shapefile, to be able to visualize the point locations
setwd("~/land_management_analysis/data/CPVT_boundary")
CPVT_polygon <- read_sf("CPVT_boundary.shp")
CPVT_polygon <- st_as_sf(CPVT_polygon)
plot(CPVT_polygon)

#creating buffers around the trails, which I can use to confirm the trees closest to the trails
trails_buffer_CPVT <- st_buffer(trails_CPVT, 10) #10 m buffer

#creating bounding boxes for the space around the trails (for Poisson point modelling later)
CPVT_box <- st_bbox(trails_CPVT_trans)
#I also use the whole of CPVT as a bounding box later on, when creating distance-to-trail rasters

#let's plot what we have so far!
ggplot()+  #plotting the trail shapefile, the buffer, and the tree points
  geom_sf(data = trails_buffer_CPVT)+
  geom_sf(data = trails_CPVT_trans)+ #is this working???????
  geom_sf(data = CPVT_fixed_field_data_processed_sf)

#This will make a much prettier version of the map above, because it is specifically designed to work with maps rather than plots
tm_shape(CPVT_polygon) + 
  tm_polygons(fill = "bisque", col = "bisque") +
  tm_shape(trails_CPVT) +
  tm_polygons(fill = "burlywood4", col = "burlywood4") +
  tm_shape(CPVT_fixed_field_data_processed_sf) +
  tm_dots(fill = "green4") +
  tm_title("CPVT butternuts (green) and CPVT trails (brown)") +
  tm_grid() +
  tm_graticules()

#creating a boundary box for LM with the UTM 18 N min and max lat lon values and then turning it into a simple feature geometry
CPVT_fixed_field_data_processed_box <- CPVT_fixed_field_data_processed_sf_transformed %>%
  filter(site == "CPVT") %>%
  st_bbox %>%
  st_as_sfc()

#histograms
ggplot(CPVT_fixed_field_data_processed) + #DBH
  geom_histogram(aes(x = dbh_cm))+
  xlab("Diameter at breast height")+
  ylab("Frequency")

#---------------------------------
#Distance-to-trail raster building
#THIS TAKES A LONG TIME TO RUN!!!!

#Creating distance to trail calculations!
#turning the trail polygons into multilinestring objects and then into a raster, to be able to later calculate the distances
trails_CPVT_trans_points <- st_cast(trails_CPVT_trans, "MULTILINESTRING") #turning the polyline of the trails into a multilinestring object
trails_CPVT_trans_point_raster <- st_rasterize(trails_CPVT_trans_points) #creating a raster out of the trail linestring object
plot(trails_CPVT_trans_point_raster, breaks = "equal") #plotting the trail multilinestring object

#doing the same with the buffers around the trails, created just above
trails_buffer_CPVT_trans <- st_transform(trails_buffer_CPVT, crs = 26918)
trails_buffer_CPVT_trans <- st_sf(geometry = trails_buffer_CPVT_trans)
trails_CPVT_buffer_trans_outline <- st_cast(trails_buffer_CPVT_trans, "MULTILINESTRING")
trails_buffer_CPVT_point_raster <- st_rasterize(trails_CPVT_buffer_trans_outline)
plot(trails_buffer_CPVT_point_raster)

#generating a distance to trail raster with the distances of each cell in the buffer raster from the trail edge points, and the trail raster cells are set to a distance of 0 m
trails_buffer_CPVT_point_raster[is.na(trails_buffer_CPVT_point_raster[])] <- 0  #making sure the cell that are not the trail buffer multilinestring raster have a 0 value
load("dist_near_trails_buffer_CPVT.rData")
#dist_near_trails_buffer_CPVT <- dist_to_nearest(trails_buffer_CPVT_point_raster, trails_CPVT_trans_points, progress = T) #creating a raster of the distances of each cell in the buffer raster to the linestring object of the trail polygons, this can take a while to run
#save(dist_near_trails_buffer_CPVT, file="dist_near_trails_buffer_CPVT.rData")
plot(dist_near_trails_buffer_CPVT) #plotting the distance to river raster

#trying instead of using the buffer to define how big the raster can be, use the whole box CPVT_fixed_field_data_processed_box
trails_fixed_box_CPVT_trans <- st_transform(CPVT_fixed_field_data_processed_box, crs = 26918)
trails_fixed_box_CPVT_trans <- st_sf(geometry = trails_fixed_box_CPVT_trans)
trails_CPVT_fixed_box_trans_outline <- st_cast(trails_fixed_box_CPVT_trans, "MULTILINESTRING")
trails_fixed_box_CPVT_point_raster <- st_rasterize(trails_CPVT_fixed_box_trans_outline)
plot(trails_fixed_box_CPVT_point_raster)

trails_fixed_box_CPVT_point_raster[is.na(trails_fixed_box_CPVT_point_raster[])] <- 0
load("dist_near_trails_fixed_box_CPVT.rData")
#if the line above isn't working / you want to make the raster over again
#dist_near_trails_fixed_box_CPVT <- dist_to_nearest(trails_fixed_box_CPVT_point_raster, trails_CPVT_trans_points, progress = T)
#save(dist_near_trails_fixed_box_CPVT, file="dist_near_trails_fixed_box_CPVT.rData")
plot(dist_near_trails_fixed_box_CPVT)

## Making it so the cells in the distance raster within or overlapping with the river raster are assigned 1 
#Assigning points within and overlapping with the river to be "true"
CPVT_points_intersects_trails <- st_intersects(CPVT_fixed_field_data_processed, trails_CPVT_trans, sparse = F) #creating a list of true or falses for whether points intersect the river shapefiles
CPVT_fixed_field_data_processed_intersects_trails <- cbind(CPVT_fixed_field_data_processed, CPVT_points_intersects_trails) #binding the list of true or falses with the point data
#printing the river polygon and the tree points, colored by whether are or aren't within or overlapping with the river
ggplot()+
  geom_sf(data=trails_CPVT_trans)+
  geom_sf(data=CPVT_fixed_field_data_processed)+
  geom_sf(data=CPVT_fixed_field_data_processed_intersects_trails, aes(color = CPVT_points_intersects_trails))


## Extracting distance to trail for each tree using the distance to trail raster

CPVT_distance_data_pts <- st_extract(dist_near_trails_buffer_CPVT, CPVT_fixed_field_data_processed) #extracting distance to trail for each tree
CPVT_fixed_field_data_processed_distance  <- cbind(CPVT_fixed_field_data_processed, CPVT_distance_data_pts) #binding the distance to trail data for each point to the CPVT point dataframe

## Assigning all points within/overlapping trail to distances of 0

CPVT_fixed_field_data_processed_distance <- CPVT_fixed_field_data_processed_distance %>% 
  mutate(d = case_when((CPVT_fixed_field_data_processed_intersects_trails$CPVT_points_intersects_trails == T) ~ 0,  #assigns 0 to points within trail
                       (CPVT_fixed_field_data_processed_intersects_trails$CPVT_points_intersects_trails == F) ~ d)) #to points outside of trail, it leaves the original distance value