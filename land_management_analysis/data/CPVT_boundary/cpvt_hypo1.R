# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%Looking to see if butternut is more clumped or dispersed than by random%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#This script is for investigating whether Quercus brandegeei seeds are predominantly dispersed either by heavy rainfall or gravity, 
#which may impact both the spatial distribution and the genetic structure of a given population. 
#If they were dispersed by heavy rainfall, we would expect more dispersal than at random. 
#If they were dispersed by gravity, we would expect more clumping of trees than at random. 
#To test this, we used a Ripley's K and compared the known tree locations to randomly generated locations
#produced in either convex hulls, the river shapefiles, and buffers around the river 
#(all three were attempted separately and compared) to determine whether the 
#tree points were more clustered or dispersed than we would expect at random. 
#We then used an Average Nearest Neighbor Analysis (ANN) to support with p-values whether the points 
#seem more clustered or dispersed than at random. 
#Finally, we used Poisson Point Models, to see if models that take into account the effect that the inverse 
#distance to the river of the points has on the placement of the trees better explains the distribution of 
#the points than if they were distributed at random.

# The Ripley's K and ANN analyses both test whether the trees are more clumped or dispersed than at random, whereby 
#the ANN is the only test of the two that provides a p-value of the two. 
#The PPM test is the only test of the three that sees whether the river itself has an influence on the 
#points distribution.

# The script is broken into sections of 
#1) loading and processing the packages and spatial data for the trees in the Las Matancitas,
#San Dionisio, and La Cobriza populations, 
#2) running the Ripley's K analysis, 
#3) running the Average Nearest Neighber (ANN) analysis,
#4) running the Poisson Point Model Analysis. 

# NOTE: Uncomment and run line 51, sourcing Data_Processing_Script.R, if the line has not yet to be run across any of the scripts/the environment has been cleared 

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

# loading in the processed tree data 
# NOTE: uncomment and run line 51, sourcing Data_Processing_Script.R, if the line has not yet to be run across any of the scripts/the environment has been cleared 
# source("./analyses/Data_Processing_Script.R")

#ensuring there is a column from latitude and longitude in the populations transformed dataframe because those columns are needed in "hypothesis_1_clumping_CLEANED.R" 

#CPVT_fixed_field_data_processed <- fixed_field_data_processed %>%
#  filter(site == "CPVT")

#### Ripley's K Analysis (version with box, convex hull, and 20 m buffer around river) ####

## Plotting the Baja/tree Polygons and Creating Bounding Boxes using tree points of each population ##

# all trees

# #finding minimum and maximum lat and long values for visualizing the tree locations with *1.002 wiggle room
# min_all_locality_long <- min(fixed_field_data_processed$gps_w)*1.0002
# max_all_locality_long <- max(fixed_field_data_processed$gps_w) - (max(fixed_field_data_processed$gps_w) *.0002)
# min_all_locality_lat <- min(fixed_field_data_processed$gps_n)*1.02
# max_all_locality_lat <- max(fixed_field_data_processed$gps_n) - (max(fixed_field_data_processed$gps_n)*.02)
# 
# #plotting the BCS polygon with all of the tree points, colored by locality
# ggplot(data = CPVT_polygon) +
#   geom_sf() +
#   geom_sf(data = fixed_field_data_processed_sf_transformed, aes(color = site)) + 
#   coord_sf(xlim = c(min_all_locality_long, max_all_locality_long), 
#            ylim = c(min_all_locality_lat, max_all_locality_lat))+
#   theme_classic()
# 
# #creating BCS boundary shapefile, turning sf of all points into sfc
# fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
#   st_as_sfc()
# 
# #creating a boundary box of LM with the UTM 12 N min and max lat lon values and then 
# #turning it into a simple feature geometry
# LM_fixed_field_data_processed_box <- fixed_field_data_processed_sf_transformed %>%
#   filter(Locality == "LM") %>%
#   st_bbox %>%
#   st_as_sfc()

# plotting the points, river, and bounding box
ggplot(CPVT_fixed_field_data_processed_box)+
  geom_sf() +
  geom_sf(data = trails_CPVT_trans) +
  geom_sf(data = CPVT_fixed_field_data_processed_sf)

## Creating Convex Hulls using tree points of each population ##

# creating one geometry with the tree points using st_union and then creating a convex hull with st_convex_hull

# LM
trails_CPVT_convex_hull <- st_convex_hull(st_union(CPVT_fixed_field_data_processed_sf))  
ggplot(trails_CPVT_convex_hull)+
  geom_sf() +
  geom_sf(data = trails_CPVT_trans) +
  geom_sf(data = CPVT_fixed_field_data_processed_sf)
#wow. so here's the thing, is this supposed to be around the trails or around the trees? RN it is around the trees

## Calculating the Ripley's K ##

# For all Ripley's K:
# 1) Turn the bounding box into a window object
# 2) Create the Poisson Point Pattern within the window (a randomized version of the point locations)
# 3) Calculate the Ripley's K function (Kest), using the Ripley's correction (edge correction)
# 4) Plot the Ripley's K function output. Dotted red line is the randomized point locations. 
# r is the distance from a center point (radius)
# K(r) is the expected number of random points within the radius of a point
# If the K(r) of our known points (black line) is lower than the randomized points (red line) than the points are more 
#dispersed than expected at random
# If the K(r) of our known points (black line) is higher than the randomized points (red line) than the points are more 
#clustered than expected at random

# Ripley's K for LM 

# Bounding Box
win <- as.owin(CPVT_fixed_field_data_processed_box) #turning the box into a window
CPVT_ppp <- as.ppp(st_coordinates(CPVT_fixed_field_data_processed_sf), W = win) #creating the poisson point pattern for LM
plot(CPVT_ppp, pch = 16, cex = 0.5) # plotting the randomized point pattern
K <- Kest(CPVT_ppp, correction = "Ripley") # Ripley's K function, using the isotropic/Ripley correction
plot(K, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) # plotting Ripley's output

# Convex Hull
CPVT_win_convex <- as.owin(trails_CPVT_convex_hull)  #turning the convex hull into a window
CPVT_ppp_convex <- as.ppp(st_coordinates(CPVT_fixed_field_data_processed_sf), W = CPVT_win_convex) #creating the poisson point pattern for lm
plot(CPVT_ppp_convex, pch = 16, cex = 0.5) # plotting the randomized point pattern
CPVT_k_convex <- Kest(CPVT_ppp_convex, correction = "Ripley") # Ripley's K function, using the isotropic/Ripley correction
plot(CPVT_k_convex, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE))  # plotting Ripley's output

############################none of this works for me. All points are out of the window.#################
# # Buffer River (100 m)
# CPVT_win_buffer <- as.owin(trails_buffer_CPVT) #turning the buffer into a window
# CPVT_ppp_buffer <- as.ppp(st_coordinates(CPVT_fixed_field_data_processed_sf), W = CPVT_win_buffer) #creating the poisson point pattern for lm
# plot(CPVT_ppp_buffer, pch = 16, cex = 0.5) # plotting the randomized point pattern
# CPVT_k_buffer <- Kest(CPVT_ppp_buffer, correction = "Ripley") # Ripley's K function, using the isotropic/Ripley correction
# plot(CPVT_k_buffer, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE))  # plotting Ripley's output
# 
# # Making a nicer version of the plot for papers/presentations
# plot(CPVT_k_buffer, main=NULL, las=1, ylab = "", xlab = "", legendargs=list(cex=0.8, xpd=TRUE), 
#      yaxp = c(0, 220000, 10))#legend inside of the plot
# par(mar = c(6,6,6,6))
# title(ylab = bquote(italic("K(r), Las Matancitas trees")), cex = 1.2, line = 4)
# title(xlab = bquote(italic("r (m)")),)
###########################################################################################################



#### ANN Analysis (test for clustering/dispersion) ####

# For all ANN Analyses (for each population and for convex hulls/buffers/and river shapefiles):
# 1) Find Average Nearest Neighbor Value for the trees of a population
# 2) Simulate a randomized distribution of points, calculate and store the average nearest neighbor, 599 times
# 3) Make sure the randomized points have the correct coordinate reference system for mapping
# 4) Plot the Randomized Points, convex hull, and river shapefile
# 5) Plot the histogram of simulated randomized average nearest neighbor values with a line for the actuall ANN value
# 6) Calculate and store the p-value 

# The convex hull windows were used because they produced the most close to the original point distances as possible.

# In this analysis, we do not control for the presence of the river as something that would influence the distribution of points

# When the average nearest neighbor is to the left of the histogram the points are more clustered than expected at random
# When the average nearest neighbor is to the right of the histogram the points are more dispersed than expected at random

# creating the river rasters for the ANN analysis 


# LM

# creating the rasters that will be used for point generation later

#making sure the trails polygon and distance to trails raster have the same projection
trails_vect_CPVT <- project(vect(trails_CPVT), rast(dist_near_trails_buffer_CPVT))

#creating the corrected distance to trails column where the values inside the cells touching the polygon equal 0
dist_near_trails_buffer_CPVT_corrected <- rasterize(trails_vect_CPVT, rast(dist_near_trails_buffer_CPVT), field=0, update=TRUE, touches=TRUE)

#turning the distance to trails correct raster into a stars object
dist_near_trails_buffer_CPVT_corrected_stars <- st_as_stars(dist_near_trails_buffer_CPVT_corrected)

#creating the inverse of the distance raster so that the higher values are closer to the trails and the values are between 0-1
dist_near_trails_buffer_CPVT_inverse <- dist_near_trails_buffer_CPVT_corrected_stars %>%
  st_as_sf() %>%
  mutate(layer = case_when(layer >= 1 ~ 1/layer,
                           layer > 0 & layer < 1 ~ 1,
                           layer == 0 ~ 1)) %>%
  st_rasterize()

# creating the cropped versions ## 

#creating a raster out of the inverse distance stars object
dist_near_trails_buffer_CPVT_inverse_im_raster <- rast(dist_near_trails_buffer_CPVT_inverse)

#projecting the inverse distance raster to match the other crs
crs(dist_near_trails_buffer_CPVT_inverse_im_raster) <- crs(rast(dist_near_trails_buffer_CPVT_inverse))

#cropping the distance to trails raster
#dist_near_trails_buffer_CPVT_inverse_cropped <- raster::crop(dist_near_trails_buffer_CPVT_inverse_im_raster, trails_buffer_CPVT, mask = T)

#trimming off the NAs
#dist_near_trails_buffer_CPVT_inverse_cropped <- trim(dist_near_trails_buffer_CPVT_inverse_cropped)



#creating ANN Analysis function

ANN_analysis <- function(population, window) {
  if (population == "CPVT") {
    ppp <- CPVT_ppp #assigning poisson point pattern 
    dataframe <- CPVT_fixed_field_data_processed_sf #assigning dataframe
    
    #window selection
    if (window == "Convex Hull"){ #ANN without controlling for river
      selected_window <- trails_CPVT_convex_hull
    } else if (window == "Just River"){ #ANN with controlling for river
      selected_window <- trails_CPVT_trans_point_raster
    } else if (window == "Inside, On, and Outside River"){
      selected_window <- st_as_stars(dist_near_trails_buffer_CPVT_inverse_cropped)
    } else if (window == "On and Inside River"){
      selected_window <- st_rasterize(trails_CPVT_trans)
    }
  }
  
  #calculating the average nearest neighbor value for the entire population of trees
  ann.p <- mean(nndist(ppp, k=1))
  ann.p
  #12.54147
  
  #simulating the random points and calculating the average nearest neighbor for each 566 permutations
  if (window == "Convex Hull"){ 
    #simulation to create a list of ANN from randomly placed points
    n <- 566L #defines the number of simulations
    ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
    for (i in 1:n){
      rand.p <- rpoint(n=length(dataframe), win = as.owin(selected_window)) # generating the random points within the convex hull window
      ann.r[i] <- mean(nndist(rand.p, k=1)) #for each simulated random distribution of points it calculates the mean ANN across all of the trees
    } #for the number of points at LM, it assigns a random point within the convex hull window
  } else { 
    #ANN analysis controlling for river
    n <- 599L #defines the number of simulations
    ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
    for (i in 1:n){ 
      rand.p <- rpoint(n=length(dataframe), f = as.im(selected_window)) # generating the random points within the window
      ann.r[i] <- mean(nndist(rand.p, k=1)) #for each simulated random distribution of points it calculates the mean ANN across all of the trees
    } #for the length of the number of points at LM, it assigns a random point on top of the river's edge while controlling for the river's edge
  }
  
  #adding the UTM 18 crs to rand.p
  rand.p.crs <- rand.p %>% 
    st_as_sf()%>%
    st_set_crs(26918)
  
  #calculating pseudo p-value for 
  total = 0  #set empty vaue
  for (i in 1:length(ann.r)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
    if (ann.r[i] < ann.p){
      total = total + 1
    }
  } #add number of values of in the random set of ANN values that are less than our mean ANN
  p_value <- total / length(ann.r)
  
  plot(rand.p)
  print(paste0("Average Nearest Neighbor for Original Trees: ", ann.p))
  print(paste0("P-Value: ", p_value))
  return(list(random_points = rand.p.crs, observed_ANN = ann.p, ann.r = ann.r, p.value = p_value)) #the proportion of random ANNs that are less than our ANN (p-value)
  
}

## Convex Hull

# LM

CPVT_ANN_Anlysis <- ANN_analysis("CPVT", "Convex Hull")
CPVT_ANN_Anlysis #first index is the ANN value, the second is the left-tailed p-value

#plotting the randomly generated points, tree points, and the river
ggplot()+ 
  geom_sf(data=trails_CPVT_trans)+ #plotting the river
  geom_sf(data=CPVT_fixed_field_data_processed_sf, aes(col = "red")) + #plotting the tree points
  geom_sf(data=CPVT_ANN_Anlysis$random_points, fill = NA) #plotting the random points

#creating a histogram of the ANN Simulation Results
as_tibble(CPVT_ANN_Anlysis$ann.r) %>%  #turns the list of ann values from the simulations of random points and turns it into a tibble/dataframe
  ggplot()+
  geom_histogram(aes(x = value), fill = "dodgerblue1", color = "black", bins = 50) +
  xlim(range(CPVT_ANN_Anlysis$observed_ANN, CPVT_ANN_Anlysis$ann.r)) + #sets the limit of the xaxis to encompass the ANN for our trees and histogram of ANNs from the simulation
  geom_vline(xintercept=CPVT_ANN_Anlysis$observed_ANN, col = "red") + #adds a verticle line of our tree'\s' ANN
  xlab("ANN")+
  theme_classic()

#### ANN Analysis (test for clustering/dispersion) while controlling for the river ####

# The steps for this ANN are the same as previously, except we use three different versions of the windows in which we generate random points
# with varying levels of control for the river to see if the points still seem significantly clustered despite the 
# presence of the rivers (similar to the PPM analysis later).

# The three ways of controlling for the river include 
# a) controlling for the river border (using a river multipoint raster window), 
# b) controlling for on, inside, and around the river (using an inverse distance raster window), and 
# c) controlling for on and inside the river (using a raster of the river window)


# To do this, we add new steps in the beginning 
# 1) Create rasters of the river shapefile, river buffer, and create a raster with the inverse distances of 
#each cell to the river shapefile (closer cells are weighted higher)
#cells within a certain distance of the river equal 1 and the other points equals 1/distance
# 2) Run the simulations whereby the windows either use the river border raster, the inverse distance raster where 
#the randomized points are placed more likely based on the raster or the higher cell weights, and the river polygon raster

###  LM

## Version of ANN analysis controlling for the river with just the river multipoint 
CPVT_ANN_Anlysis_river <- ANN_analysis("CPVT", "Just River")
CPVT_ANN_Anlysis_river #first index is the ANN value, the second is the left-tailed p-value

#plotting the randomly generated points, tree points, and probability/distance raster
ggplot()+ 
  geom_stars(data=trails_CPVT_trans_point_raster)+ #plotting the river edge raster
  geom_sf(data=CPVT_fixed_field_data_processed_sf, aes(col = "red"))+ #plotting the tree points
  geom_sf(data=CPVT_ANN_Anlysis_river$random_points$geom, fill = NA) #plotting the random points

#graphing the histogram of simulated ANN values and the mean ANN from our trees
as_tibble(CPVT_ANN_Anlysis_river$ann.r) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "dodgerblue1", color = "black", bins = 50) + 
  xlim(range(CPVT_ANN_Anlysis_river$observed_ANN, CPVT_ANN_Anlysis_river$ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=CPVT_ANN_Anlysis_river$observed_ANN, col = "red") + #plotting our tree's mean ANN
  xlab("ANN") +
  theme_classic()

## Version of ANN analysis controlling for the river with inside, on, and outside the river
#This didn't work as I wasn't able to make the cropping work
CPVT_ANN_Anlysis_inside_on_outside_river <- ANN_analysis("CPVT", "Inside, On, and Outside River")
CPVT_ANN_Anlysis_inside_on_outside_river #first index is the ANN value, the second is the left-tailed p-value

#plotting the randomly generated points, tree points, and probability/distance raster
ggplot()+ 
  geom_stars(data=na.omit(st_as_stars(dist_near_river_buffer_LM_inverse_cropped), aes(fill = layer)))+ #plotting the distance inverse raster 
  scale_fill_distiller(palette = "Blues", na.value = "transparent", trans = "reverse")+
  geom_sf(data=st_cast(LM_ANN_Anlysis_inside_on_outside_river$random_points$geom, "POINT"), alpha = 0.5, aes(color = "Randomly Generated"), fill = NA, shape = 16) + #plotting the random points
  geom_sf(data=LM_fixed_field_data_processed_sf, aes(color = "Actual Trees"), shape = 16, alpha = 0.5)+ #plotting the tree points
  labs(color = "Actual Trees", fill = "Inverse Distance (1/m)", 
       x = "Longitude", 
       y = "Latitude")+
  scale_color_manual(
    name = "Trees",
    values = c("Actual Trees" = "red", 
               "Randomly Generated" = "black"))+
  theme_minimal()+
  # guides(color = guide_legend(override.aes = list(shape = c(16,16), linetype = 0)))+
  labs(title = "Las Matancitas")+
  theme_classic() +
  theme(title=element_text(size=15), 
        axis.text=element_text(size=15),  axis.title.x =element_text(size= 15),
        axis.title.y =element_text(size= 15),
        text = element_text(family = "serif"))

#graphing the histogram of simulated ANN values and the mean ANN from our trees
as_tibble(LM_ANN_Anlysis_inside_on_outside_river$ann.r) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "skyblue", color = "black", bins = 50) + 
  xlim(range(LM_ANN_Anlysis_inside_on_outside_river$observed_ANN, LM_ANN_Anlysis_inside_on_outside_river$ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=LM_ANN_Anlysis_inside_on_outside_river$observed_ANN, col = "red") + #plotting our tree's mean ANN
  xlab("Average Nearest Neighbor") +
  ylab("Frequency")+
  labs(title = "Las Matancitas")+
  geom_text(aes(label = round(LM_ANN_Anlysis_inside_on_outside_river$observed_ANN, 2)), x = 6.7, y = 40, color = "red") +
  theme_classic() +
  theme(title=element_text(size=15), 
        axis.text=element_text(size=15),  axis.title.x =element_text(size= 15),
        axis.title.y =element_text(size= 15),
        text = element_text(family = "serif"))

## Version of ANN analysis controlling for the river with on and inside the river 

LM_ANN_Anlysis_on_inside_river <- ANN_analysis("LM", "On and Inside River")
LM_ANN_Anlysis_on_inside_river #first index is the ANN value, the second is the left-tailed p-value

#plotting the randomly generated points, tree points, and river raster
ggplot()+ 
  geom_stars(data=st_rasterize(river_LM_trans))+ #plotting the river raster 
  geom_sf(data=LM_fixed_field_data_processed_sf, aes(col = "red"))+ #plotting the tree points
  geom_sf(data=LM_ANN_Anlysis_on_inside_river$random_points, fill = NA) #plotting the random points

#graphing the histogram of simulated ANN values and the mean ANN from our trees
as_tibble(LM_ANN_Anlysis_on_inside_river$ann.r) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "skyblue", color = "black", bins = 50) + 
  xlim(range(LM_ANN_Anlysis_on_inside_river$observed_ANN, LM_ANN_Anlysis_on_inside_river$ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=LM_ANN_Anlysis_on_inside_river$observed_ANN, col = "red") + #plotting our tree's mean ANN
  xlab("Average Nearest Neighbor") +
  ylab("Frequency")+
  theme_classic()


#### PPM analysis ####

# for every PPM analysis 

# 1) generate an image object of distance raster 
# 2) create the Poisson Point Model using ppm() function with the river influencing the location of the points (Alternative Hypothesis Model)
# 3) create the Poisson Point Model using ppm() function with the river not influencing the location of the points (Null Hypothesis Model)
# 4) Use an ANOVA likelihood Ratio Test to compare the Alternate and Null hypotheses
# 5) Plot the influence of the river as the distance to the river decreases (inverse distance)

#Test for CPVT

#creating the image of the distance to trails stars
dist_near_trails_buffer_CPVT_inverse_im <- as.im(st_as_stars(dist_near_trails_buffer_CPVT_inverse))

#Alternative hypothesis, seeing if the distance to the trails's edge influences the tree point placement
PPM1 <- ppm(Q = as.ppp(CPVT_fixed_field_data_processed_sf) ~ dist_near_trails_buffer_CPVT_inverse_im) 
PPM1

#null hypothesis, no change in the trend of the points
PPM0 <- ppm(as.ppp(CPVT_fixed_field_data_processed_sf) ~ 1)
PPM0
PPM0 <- ppm(CPVT_ppp ~ 1)
PPM0

#using a likelihood ratio test to compare the alternative and null models
anova(PPM0, PPM1, test="LRT")

#plotting the alternative model
plot(effectfun(PPM1, "dist_near_trails_buffer_CPVT_inverse_im", se.fit = TRUE), main = "Distance to trails of CPVT",
     ylab = "Butternut Trees", xlab = "Inverse Distance to trails", legend = FALSE)

#Test for LC

#creating the image of the distance to river stars
dist_near_river_buffer_LC_inverse_im <- as.im(dist_near_river_buffer_LC_inverse)

#Alternative hypothesis, seeing if the distance to the river's edge influences the tree point placement
PPM1 <- ppm(Q = as.ppp(LC_fixed_field_data_processed_sf) ~ dist_near_river_buffer_LC_inverse_im) 
PPM1

#null hypothesis, no change in the trend of the points
PPM0 <- ppm(as.ppp(LC_fixed_field_data_processed_sf) ~ 1)
PPM0

#using a likelihood ratio test to compare the alternative and null models
anova(PPM0, PPM1, test="LRT")

#plotting the alternative model
plot(effectfun(PPM1, "dist_near_river_buffer_LC_inverse_im", se.fit = TRUE), main = "Distance to River of La Cobriza",
     ylab = "Quercus brandegeei Trees", xlab = "Inverse Distance to River", legend = FALSE)

#Test for SD

#creating the image of the distance to river stars
dist_near_river_buffer_SD_inverse_im <- as.im(dist_near_river_buffer_SD_inverse)

#Alternative hypothesis, seeing if the distance to the river's edge influences the tree point placement
PPM1 <- ppm(Q = as.ppp(SD_fixed_field_data_processed_sf) ~ dist_near_river_buffer_SD_inverse_im) 
PPM1

#null hypothesis, no change in the trend of the points
PPM0 <- ppm(as.ppp(SD_fixed_field_data_processed_sf) ~ 1)
PPM0

#using a likelihood ratio test to compare the alternative and null models
anova(PPM0, PPM1, test="LRT")

#plotting the alternative model
plot(effectfun(PPM1, "dist_near_river_buffer_SD_inverse_im", se.fit = TRUE), main = "Distance to River of San Dionisio",
     ylab = "Quercus brandegeei Trees", xlab = "Inverse Distance to River", legend = FALSE)


# making examples of random point distributions vs. points only along the river's edge for presentation 

points_box = sf::st_sample(SD_box, size=50) #randomizing points onlu in population bbox
points_river = sf::st_sample(river_SD_trans_points, size=50) #randomizing points along river's edge

#plotting the randomized box points
ggplot()+
  geom_sf(data = river_SD_trans)+
  geom_sf(data = points_box, size = 2)+
  theme_classic()

#plotting the randomized river's edge points
ggplot()+
  geom_sf(data = river_SD_trans)+
  geom_sf(data = points_river, size = 2)+
  theme_classic()

## Using more cropped versions

# LM 

#creating the image of the distance to river stars
dist_near_river_buffer_LM_inverse_im_cropped <- as.im(st_as_stars(dist_near_river_buffer_LM_inverse_cropped))

#creating a window of the river buffer for cropping purposes
river_buffer_LM_W <- spatstat.geom::as.owin(st_as_sf(river_buffer_LM))

#creating a poison point model object of our known trees 
LM_fixed_field_data_processed_ppp <- as.ppp(LM_fixed_field_data_processed_sf)

#cropping the river with the river buffer window
LM_fixed_field_data_processed_ppp_crop <- LM_fixed_field_data_processed_ppp[river_buffer_LM_W]

#Alternative hypothesis, seeing if the distance to the river's edge influences the tree point placement
PPM1 <- ppm(Q = LM_fixed_field_data_processed_ppp_crop ~ dist_near_river_buffer_LM_inverse_im_cropped)
PPM1

#null hypothesis, no change in the trend of the points
PPM0 <- ppm(LM_fixed_field_data_processed_ppp_crop ~ 1)
PPM0

#using a likelihood ratio test to compare the alternative and null models
anova(PPM0, PPM1, test="LRT")

#plotting the alternative model
plot(effectfun(PPM1, "dist_near_river_buffer_LM_inverse_im_cropped", se.fit = TRUE), main = "Distance to River of Las Matancitas",
     ylab = "Quercus brandegeei Trees", xlab = "Inverse Distance to River", legend = FALSE)

#Test for LC

#creating the image of the distance to river stars
dist_near_river_buffer_LC_inverse_im_cropped <- as.im(st_as_stars(dist_near_river_buffer_LC_inverse_cropped))

#creating a window of the river buffer for cropping purposes
river_buffer_LC_W <- spatstat.geom::as.owin(st_as_sf(river_buffer_LC))

#creating a poison point model object of our known trees 
LC_fixed_field_data_processed_ppp <- as.ppp(LC_fixed_field_data_processed_sf)

#cropping the river with the river buffer window
LC_fixed_field_data_processed_ppp_crop <- LC_fixed_field_data_processed_ppp[river_buffer_LC_W]

#Alternative hypothesis, seeing if the distance to the river's edge influences the tree point placement
PPM1 <- ppm(Q = LC_fixed_field_data_processed_ppp_crop ~ dist_near_river_buffer_LC_inverse_im_cropped) 
PPM1

#null hypothesis, no change in the trend of the points
PPM0 <- ppm(LC_fixed_field_data_processed_ppp_crop ~ 1)
PPM0

#using a likelihood ratio test to compare the alternative and null models
anova(PPM0, PPM1, test="LRT")

#plotting the alternative model
plot(effectfun(PPM1, "dist_near_river_buffer_LC_inverse_im_cropped", se.fit = TRUE), main = "Distance to River of La Cobriza",
     ylab = "Quercus brandegeei Trees", xlab = "Inverse Distance to River", legend = FALSE)

#Test for SD

#creating the image of the distance to river stars
dist_near_river_buffer_SD_inverse_im_cropped <- as.im(st_as_stars(dist_near_river_buffer_SD_inverse_cropped))

#creating a window of the river buffer for cropping purposes
river_buffer_SD_W <- spatstat.geom::as.owin(st_as_sf(river_buffer_SD))

#creating a poison point model object of our known trees 
SD_fixed_field_data_processed_ppp <- as.ppp(SD_fixed_field_data_processed_sf)

#cropping the river with the river buffer window
SD_fixed_field_data_processed_ppp_crop <- SD_fixed_field_data_processed_ppp[river_buffer_SD_W]

#Alternative hypothesis, seeing if the distance to the river's edge influences the tree point placement
PPM1 <- ppm(Q = SD_fixed_field_data_processed_ppp_crop ~ dist_near_river_buffer_SD_inverse_im_cropped) 
PPM1

#null hypothesis, no change in the trend of the points
PPM0 <- ppm(SD_fixed_field_data_processed_ppp_crop ~ 1)
PPM0

#using a likelihood ratio test to compare the alternative and null models
anova(PPM0, PPM1, test="LRT")

#plotting the alternative model
plot(effectfun(PPM1, "dist_near_river_buffer_SD_inverse_im_cropped", se.fit = TRUE), main = "Distance to River of San Dionisio",
     ylab = "Quercus brandegeei Trees", xlab = "Inverse Distance to River", legend = FALSE)



## Using the non-inverse distance to river raster ##


#### Creating the distance to river rasters where everything inside the river equals 1 ####

#LM

#making sure the river polygon and distance to river raster have the same projection
river_vect_LM <- project(vect(st_as_sf(river_LM)), rast(dist_near_river_buffer_LM))

#creating the corrected distance to river column where the values inside the cells touching the polygon equal 0
dist_near_river_buffer_LM_corrected <- rasterize(river_vect_LM, rast(dist_near_river_buffer_LM), field=0, update=TRUE, touches=TRUE)

#making sure the projections are identifical for the buffer and the raster
river_buffer_LM <- project(river_buffer_LM, rast(dist_near_river_buffer_LM_corrected))

#cropping the distance to river raster
dist_near_river_buffer_LM_corrected_cropped <- raster::crop(dist_near_river_buffer_LM_corrected, river_buffer_LM, mask = T)

#trimming off the NAs
dist_near_river_buffer_LM_corrected_cropped <- trim(dist_near_river_buffer_LM_corrected_cropped)

#LC

#making sure the river polygon and distance to river raster have the same projection
river_vect_LC <- project(vect(st_as_sf(river_LC)), rast(dist_near_river_buffer_LC))

#creating the corrected distance to river column where the values inside the cells touching the polygon equal 0
dist_near_river_buffer_LC_corrected <- rasterize(river_vect_LC, rast(dist_near_river_buffer_LC), field=0, update=TRUE, touches=TRUE)

#making sure the projections are identifical for the buffer and the raster
river_buffer_LC <- project(vect(river_buffer_LC), rast(dist_near_river_buffer_LC_corrected))

#cropping the distance to river raster
dist_near_river_buffer_LC_corrected_cropped <- raster::crop(dist_near_river_buffer_LC_corrected, river_buffer_LC, mask = T)

#trimming off the NAs
dist_near_river_buffer_LC_corrected_cropped <- trim(dist_near_river_buffer_LC_corrected_cropped)

#SD

#making sure the river polygon and distance to river raster have the same projection
river_vect_SD <- project(vect(st_as_sf(river_SD)), rast(dist_near_river_buffer_SD))

#creating the corrected distance to river column where the values inside the cells touching the polygon equal 0
dist_near_river_buffer_SD_corrected <- rasterize(river_vect_SD, rast(dist_near_river_buffer_SD), field=0, update=TRUE, touches=TRUE)

#making sure the projections are identifical for the buffer and the raster
river_buffer_SD <- project(vect(river_buffer_SD), rast(dist_near_river_buffer_SD_corrected))

#cropping the distance to river raster
dist_near_river_buffer_SD_corrected_cropped <- raster::crop(dist_near_river_buffer_SD_corrected, river_buffer_SD, mask = T)

#trimming off the NAs
dist_near_river_buffer_SD_corrected_cropped <- trim(dist_near_river_buffer_SD_corrected_cropped)

ggplot()+
  geom_raster(data = as.data.frame(dist_near_river_buffer_LM_corrected_cropped, xy=T), aes(x=x, y=y, fill = layer))+
  #geom_sf(data = river_LM_trans)+
  geom_sf(data = LM_fixed_field_data_processed_sf)

#Test for LM

#creating the image of the distance to river stars
dist_near_river_buffer_LM_corrected_im <- as.im(st_as_stars(dist_near_river_buffer_LM_corrected_cropped))

#creating a window of the river buffer for cropping purposes
river_buffer_LM_W <- spatstat.geom::as.owin(st_as_sf(river_buffer_LM))

#creating a poison point model object of our known trees 
LM_fixed_field_data_processed_ppp <- as.ppp(LM_fixed_field_data_processed_sf)

#cropping the river with the river buffer window
LM_fixed_field_data_processed_ppp_crop <- LM_fixed_field_data_processed_ppp[river_buffer_LM_W]

#Alternative hypothesis, seeing if the distance to the river's edge influences the tree point placement
PPM1 <- ppm(Q = LM_fixed_field_data_processed_ppp_crop ~ dist_near_river_buffer_LM_corrected_im, na.rm = TRUE) 
PPM1

#null hypothesis, no change in the trend of the points
PPM0 <- ppm(LM_fixed_field_data_processed_ppp_crop ~ 1)
PPM0

#using a likelihood ratio test to compare the alternative and null models
anova(PPM0, PPM1, test="LRT")

#plotting the alternative model
plot(effectfun(PPM1, "dist_near_river_buffer_LM_corrected_im", se.fit = TRUE), main = "Distance to River of Las Matancitas",
     ylab = "Quercus brandegeei Trees", xlab = "Distance to River", legend = FALSE)

#Test for LC

#creating the image of the distance to river stars
dist_near_river_buffer_LC_corrected_im <- as.im(st_as_stars(dist_near_river_buffer_LC_corrected_cropped))

#creating a window of the river buffer for cropping purposes
river_buffer_LC_W <- spatstat.geom::as.owin(st_as_sf(river_buffer_LC))

#creating a poison point model object of our known trees 
LC_fixed_field_data_processed_ppp <- as.ppp(LC_fixed_field_data_processed_sf)

#cropping the river with the river buffer window
LC_fixed_field_data_processed_ppp_crop <- LC_fixed_field_data_processed_ppp[river_buffer_LC_W]

#Alternative hypothesis, seeing if the distance to the river's edge influences the tree point placement
PPM1 <- ppm(Q = LC_fixed_field_data_processed_ppp_crop ~ dist_near_river_buffer_LC_corrected_im) 
PPM1

#null hypothesis, no change in the trend of the points
PPM0 <- ppm(LC_fixed_field_data_processed_ppp_crop ~ 1)
PPM0

#using a likelihood ratio test to compare the alternative and null models
anova(PPM0, PPM1, test="LRT")

#plotting the alternative model
plot(effectfun(PPM1, "dist_near_river_buffer_LC_corrected_im", se.fit = TRUE), main = "Distance to River of La Cobriza",
     ylab = "Quercus brandegeei Trees", xlab = "Distance to River", legend = FALSE)

#Test for SD

#creating the image of the distance to river stars
dist_near_river_buffer_SD_corrected_im <- as.im(st_as_stars(dist_near_river_buffer_SD_corrected_cropped))

#creating a window of the river buffer for cropping purposes
river_buffer_SD_W <- spatstat.geom::as.owin(st_as_sf(river_buffer_SD))

#creating a poison point model object of our known trees 
SD_fixed_field_data_processed_ppp <- as.ppp(SD_fixed_field_data_processed_sf)

#cropping the river with the river buffer window
SD_fixed_field_data_processed_ppp_crop <- SD_fixed_field_data_processed_ppp[river_buffer_SD_W]

#Alternative hypothesis, seeing if the distance to the river's edge influences the tree point placement
PPM1 <- ppm(Q = SD_fixed_field_data_processed_ppp_crop ~ dist_near_river_buffer_SD_corrected_im) 
PPM1

#null hypothesis, no change in the trend of the points
PPM0 <- ppm(SD_fixed_field_data_processed_ppp_crop ~ 1)
PPM0

#using a likelihood ratio test to compare the alternative and null models
anova(PPM0, PPM1, test="LRT")

#plotting the alternative model
plot(effectfun(PPM1, "dist_near_river_buffer_SD_corrected_im", se.fit = TRUE), main = "Distance to River of San Dionisio",
     ylab = "Quercus brandegeei Trees", xlab = "Distance to River", legend = FALSE)


# making examples of random point distributions vs. points only along the river's edge for presentation 

points_box = sf::st_sample(SD_box, size=50) #randomizing points onlu in population bbox
points_river = sf::st_sample(river_SD_trans_points, size=50) #randomizing points along river's edge

#plotting the randomized box points
ggplot()+
  geom_sf(data = river_SD_trans)+
  geom_sf(data = points_box, size = 2)+
  theme_classic()

#plotting the randomized river's edge points
ggplot()+
  geom_sf(data = river_SD_trans)+
  geom_sf(data = points_river, size = 2)+
  theme_classic()



#### Session Info ####
# 
# R version 4.4.3 (2025-02-28)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.2
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] lmtest_0.9-40          zoo_1.8-15             gdalUtilities_1.2.5   
# [4] car_3.1-5              carData_3.0-6          googledrive_2.1.1     
# [7] spdep_1.3-10           spData_2.3.4           performance_0.16.0    
# [10] ggsoiltexture_1.1.1    whitebox_2.4.3         spatialEco_2.0-3      
# [13] visreg_2.7.0           MuMIn_1.48.11          plotly_4.10.4         
# [16] mgcv_1.9-1             tmaptools_3.2          geostatsp_2.0.8       
# [19] terra_1.8-29           Matrix_1.7-2           starsExtra_0.2.8      
# [22] stars_0.6-8            abind_1.4-8            ggnewscale_0.5.1      
# [25] rstatix_0.7.2          raster_3.6-32          sp_2.2-1              
# [28] spatstat_3.3-1         spatstat.linnet_3.2-5  spatstat.model_3.3-4  
# [31] rpart_4.1.24           spatstat.explore_3.5-2 nlme_3.1-167          
# [34] spatstat.random_3.4-1  spatstat.geom_3.5-0    spatstat.univar_3.1-4 
# [37] spatstat.data_3.1-6    geomtextpath_0.1.5     PMCMRplus_1.9.12      
# [40] ggpmisc_0.6.1          ggpp_0.5.8-1           smatr_3.4-8           
# [43] sf_1.0-24              moments_0.14.1         lubridate_1.9.4       
# [46] forcats_1.0.0          stringr_1.6.0          dplyr_1.2.0           
# [49] purrr_1.2.1            readr_2.2.0            tidyr_1.3.2           
# [52] tibble_3.3.1           ggplot2_4.0.2          tidyverse_2.0.0       
# [55] emmeans_2.0.2         
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3    wk_0.9.4              rstudioapi_0.18.0    
# [4] jsonlite_2.0.0        magrittr_2.0.4        TH.data_1.1-5        
# [7] estimability_1.5.1    spatstat.utils_3.1-5  SuppDists_1.1-9.8    
# [10] farver_2.1.2          fs_1.6.7              vctrs_0.7.2          
# [13] memoise_2.0.1         usethis_3.2.1         htmltools_0.5.9      
# [16] polynom_1.4-1         curl_6.2.1            broom_1.0.12         
# [19] s2_1.1.7              BWStest_0.2.3         Formula_1.2-5        
# [22] KernSmooth_2.23-26    htmlwidgets_1.6.4     sandwich_3.1-1       
# [25] cachem_1.1.0          igraph_2.2.2          lifecycle_1.0.5      
# [28] pkgconfig_2.0.3       R6_2.6.1              fastmap_1.2.0        
# [31] digest_0.6.39         numDeriv_2016.8-1.1   colorspace_2.1-2     
# [34] tensor_1.5            pkgload_1.4.1         nngeo_0.4.8          
# [37] textshaping_1.0.0     labeling_0.4.3        lwgeom_0.2-14        
# [40] spatstat.sparse_3.1-0 timechange_0.3.0      httr_1.4.7           
# [43] polyclip_1.10-7       compiler_4.4.3        gargle_1.5.2         
# [46] remotes_2.5.0         proxy_0.4-27          withr_3.0.2          
# [49] backports_1.5.0       DBI_1.2.3             pkgbuild_1.4.8       
# [52] MASS_7.3-64           quantreg_6.1          sessioninfo_1.2.3    
# [55] classInt_0.4-11       tools_4.4.3           units_0.8-6          
# [58] goftest_1.2-3         glue_1.8.0            grid_4.4.3           
# [61] generics_0.1.4        gtable_0.3.6          tzdb_0.5.0           
# [64] class_7.3-23          data.table_1.17.0     hms_1.1.4            
# [67] pillar_1.11.1         splines_4.4.3         lattice_0.22-9       
# [70] survival_3.8-3        gmp_0.7-5             deldir_2.0-4         
# [73] SparseM_1.84-2        tidyselect_1.2.1      stats4_4.4.3         
# [76] devtools_2.4.6        stringi_1.8.7         boot_1.3-31          
# [79] lazyeval_0.2.2        codetools_0.2-20      kSamples_1.2-10      
# [82] multcompView_0.1-10   cli_3.6.5             xtable_1.8-4         
# [85] systemfonts_1.2.3     munsell_0.5.1         dichromat_2.0-0.1    
# [88] Rcpp_1.1.1            coda_0.19-4.1         XML_3.99-0.20        
# [91] parallel_4.4.3        ellipsis_0.3.2        MatrixModels_0.5-4   
# [94] Rmpfr_1.0-0           viridisLite_0.4.3     mvtnorm_1.3-6        
# [97] scales_1.4.0          e1071_1.7-16          insight_1.4.6        
# [100] crayon_1.5.3          rlang_1.1.7           multcomp_1.4-30  
# 