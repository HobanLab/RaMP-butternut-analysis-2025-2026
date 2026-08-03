# Catherine dell'Olio
# Testing distribution of butternut

#Again, I owe everything and anything to Rebecca Wanger's work on Butternut.  The specific script I draw inspiration from here i https://github.com/HobanLab/QUBR_GenGeoEcoDemoCorr/blob/b85cbd70de1591ef81ef08558da2d392a1040a42/analyses/hypothesis_1_clumping_CLEANED.R

#This script will investigate whether butternut trees follow a clumped or dispersed distribution, or if there distribution could be explained by randomness.  Clumping of butternut populations might indicate clumping of seeds when they are cached by wildlife, or indicate some survival bias in the areas where butternut seedlings are able to germinate and butternut trees able to grow (due to light or nutrient availability, for example).

#To test this, we used a Ripley's K and compared the known tree locations to randomly generated locations produced in either convex hulls, the river shapefiles, and buffers around the river (all three were attempted separately and compared) to determine whether the tree points were more clustered or dispersed than we would expect at random. 
#We then used an Average Nearest Neighbor Analysis (ANN) to support with p-values whether the points seem more clustered or dispersed than at random. 

#Finally, we used Poisson Point Models, to see if models that take into account the effect that the inverse distance to the river of the points has on the placement of the trees better explains the distribution of the points than if they were distributed at random.

# The Ripley's K and ANN analyses both test whether the trees are more clumped or dispersed than at random, whereby the ANN is the only test of the two that provides a p-value of the two. 
#The PPM test is the only test of the three that sees whether the river itself has an influence on the points distribution.

# The script is broken into sections of 
#1) loading and processing the packages and spatial data for the butternuts
#2) running the Ripley's K analysis, 
#3) running the Average Nearest Neighber (ANN) analysis,
#4) running the Poisson Point Model Analysis. 

# NOTE: Uncomment and run the line sourcing Data_Processing_Script.R (line 42 here), if the line has not yet to be run across any of the scripts/the environment has been cleared 

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
# NOTE: uncomment and run line 42, sourcing Data_Processing_Script.R, if the line has not yet to be run across any of the scripts/the environment has been cleared 
# source("./analyses/Data_Processing_Script.R")

#ensuring there is a column from latitude and longitude in the populations transformed dataframe because those columns are needed in "hypothesis_1_clumping_CLEANED.R" 

#### Ripley's K Analysis (version with box and convex hull) ####
### ILM ###

# plotting the points, river, and bounding box (encompassing all of the outlying trees)
ggplot(ILM_fixed_field_data_processed_box)+
  geom_sf() +
  geom_sf(data = trails_ILM_trans) +
  geom_sf(data = ILM_fixed_field_data_processed_sf)

# creating one geometry with the tree points using st_union and then creating a convex hull with st_convex_hull
trails_ILM_convex_hull <- st_convex_hull(st_union(ILM_fixed_field_data_processed_sf))
#let's check on that
ggplot(trails_ILM_convex_hull)+
  geom_sf() +
  geom_sf(data = trails_ILM_trans) +
  geom_sf(data = ILM_fixed_field_data_processed_sf)

## Calculating the Ripley's K: ##
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

# Ripley's K for ILM 

# Bounding Box
win <- as.owin(ILM_fixed_field_data_processed_box) #turning the box into a window
ILM_ppp <- as.ppp(st_coordinates(ILM_fixed_field_data_processed_sf), W = win) #creating the poisson point pattern for LM
plot(ILM_ppp, pch = 16, cex = 0.5) # plotting the randomized point pattern
K <- Kest(ILM_ppp, correction = "Ripley") # Ripley's K function, using the isotropic/Ripley correction
plot(K, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) # plotting Ripley's output

# Convex Hull
ILM_win_convex <- as.owin(trails_ILM_convex_hull)  #turning the convex hull into a window
ILM_ppp_convex <- as.ppp(st_coordinates(ILM_fixed_field_data_processed_sf), W = ILM_win_convex) #creating the poisson point pattern for lm
plot(ILM_ppp_convex, pch = 16, cex = 0.5) # plotting the randomized point pattern
ILM_k_convex <- Kest(ILM_ppp_convex, correction = "Ripley") # Ripley's K function, using the isotropic/Ripley correction
plot(ILM_k_convex, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE))  # plotting Ripley's output

#unlike Rebecca, I can't have the window set as a buffer around the trail network because those buffers overlap weirdly, so I skip that section

#### ANN Analysis (test for clustering/dispersion) ####

# For all ANN Analyses (for each convex hulls/buffers):
# 1) Find Average Nearest Neighbor Value for the trees of a population
# 2) Simulate a randomized distribution of points, calculate and store the average nearest neighbor, 599 times
# 3) Make sure the randomized points have the correct coordinate reference system for mapping
# 4) Plot the Randomized Points, convex hull, and trail shapefiles
# 5) Plot the histogram of simulated randomized average nearest neighbor values with a line for the actual ANN value
# 6) Calculate and store the p-value 

# The convex hull windows were used because they produced the most close-to-the-original point distances as possible.

# In this analysis, we do not control for the presence of the trails as something that would influence the distribution of points

# When the average nearest neighbor is to the left of the histogram the points are more clustered than expected at random
# When the average nearest neighbor is to the right of the histogram the points are more dispersed than expected at random

# creating the trail rasters for the ANN analysis that will be used for point generation later
#making sure the trails polygon and distance to trails raster have the same projection
trails_vect_ILM <- project(vect(trails_ILM), rast(dist_near_trails_fixed_box_ILM))

#creating the corrected distance to trails column where the values inside the cells touching the polygon equal 0
dist_near_trails_fixed_box_ILM_corrected <- rasterize(trails_vect_ILM, rast(dist_near_trails_fixed_box_ILM), field=0, update=TRUE, touches=TRUE)

#turning the distance to trails correct raster into a stars object
dist_near_trails_fixed_box_ILM_corrected_stars <- st_as_stars(dist_near_trails_fixed_box_ILM_corrected)

#creating the inverse of the distance raster so that the higher values are closer to the trails and the values are between 0-1
#this is the problem layer I believe for the ppm, because it excludes most of the values
dist_near_trails_fixed_box_ILM_inverse <- dist_near_trails_fixed_box_ILM_corrected_stars %>%
  st_as_sf() %>%
  mutate(layer = case_when(layer >= 1 ~ 1/layer,
                           layer > 0 & layer < 1 ~ 1,
                           layer == 0 ~ 1)) %>%
  st_rasterize()

#creating a raster out of the inverse distance stars object
dist_near_trails_fixed_box_ILM_inverse_im_raster <- rast(dist_near_trails_fixed_box_ILM_inverse)

#projecting the inverse distance raster to match the other crs
crs(dist_near_trails_fixed_box_ILM_inverse_im_raster) <- crs(rast(dist_near_trails_fixed_box_ILM_inverse))

#creating ANN Analysis function

ANN_analysis <- function(population, window) {
  if (population == "ILM") {
    ppp <- ILM_ppp #assigning poisson point pattern 
    dataframe <- ILM_fixed_field_data_processed_sf #assigning dataframe
    
    #window selection
    if (window == "Convex Hull"){ #ANN without controlling for trails
      selected_window <- trails_ILM_convex_hull
    } else if (window == "Just trails"){ #ANN with controlling for trails
      selected_window <- trails_ILM_trans_point_raster
    } else if (window == "On and Inside trails"){
      selected_window <- st_rasterize(trails_ILM_trans)
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
    #ANN analysis controlling for trails
    n <- 599L #defines the number of simulations
    ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
    for (i in 1:n){ 
      rand.p <- rpoint(n=length(dataframe), f = as.im(selected_window)) # generating the random points within the window
      ann.r[i] <- mean(nndist(rand.p, k=1)) #for each simulated random distribution of points it calculates the mean ANN across all of the trees
    } #for the length of the number of points at LM, it assigns a random point on top of the trails's edge while controlling for the trails's edge
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

ILM_ANN_Anlysis <- ANN_analysis("ILM", "Convex Hull")
ILM_ANN_Anlysis #first index is the ANN value, the second is the left-tailed p-value

#plotting the randomly generated points, tree points, and the trails
ggplot()+ 
  geom_sf(data=trails_ILM_trans)+ #plotting the trails
  geom_sf(data=ILM_fixed_field_data_processed_sf, aes(col = "red")) + #plotting the tree points
  geom_sf(data=ILM_ANN_Anlysis$random_points, fill = NA) #plotting the random points

#creating a histogram of the ANN Simulation Results
as_tibble(ILM_ANN_Anlysis$ann.r) %>%  #turns the list of ann values from the simulations of random points and turns it into a tibble/dataframe
  ggplot()+
  geom_histogram(aes(x = value), fill = "dodgerblue1", color = "black", bins = 50) +
  xlim(range(ILM_ANN_Anlysis$observed_ANN, ILM_ANN_Anlysis$ann.r)) + #sets the limit of the xaxis to encompass the ANN for our trees and histogram of ANNs from the simulation
  geom_vline(xintercept=ILM_ANN_Anlysis$observed_ANN, col = "red") + #adds a verticle line of our tree'\s' ANN
  xlab("ANN")+
  theme_classic()

#### ANN Analysis (test for clustering/dispersion) while controlling for the trails ####

# The steps for this ANN are the same as previously, except we use three different versions of the windows in which we generate random points
# with varying levels of control for the trails to see if the points still seem significantly clustered despite the 
# presence of the trailss (similar to the PPM analysis later).

# The three ways of controlling for the trails include 
# a) controlling for the trails border (using a trails multipoint raster window), 
# b) controlling for on, inside, and around the trails (using an inverse distance raster window), and 
# c) controlling for on and inside the trails (using a raster of the trails window)


# To do this, we add new steps in the beginning 
# 1) Create rasters of the trails shapefile, trails fixed_box, and create a raster with the inverse distances of 
#each cell to the trails shapefile (closer cells are weighted higher)
#cells within a certain distance of the trails equal 1 and the other points equals 1/distance
# 2) Run the simulations whereby the windows either use the trails border raster, the inverse distance raster where 
#the randomized points are placed more likely based on the raster or the higher cell weights, and the trails polygon raster

## Version of ANN analysis controlling for the trails with just the trails multipoint 
ILM_ANN_Anlysis_trails <- ANN_analysis("ILM", "Just trails")
ILM_ANN_Anlysis_trails #first index is the ANN value, the second is the left-tailed p-value

#plotting the randomly generated points, tree points, and probability/distance raster
ggplot()+ 
  geom_stars(data=trails_ILM_trans_point_raster)+ #plotting the trails edge raster
  geom_sf(data=ILM_fixed_field_data_processed_sf, aes(col = "red"))+ #plotting the tree points
  geom_sf(data=ILM_ANN_Anlysis_trails$random_points$geom, fill = NA) #plotting the random points

#graphing the histogram of simulated ANN values and the mean ANN from our trees
as_tibble(ILM_ANN_Anlysis_trails$ann.r) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "dodgerblue1", color = "black", bins = 50) + 
  xlim(range(ILM_ANN_Anlysis_trails$observed_ANN, ILM_ANN_Anlysis_trails$ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=ILM_ANN_Anlysis_trails$observed_ANN, col = "red") + #plotting our tree's mean ANN
  xlab("ANN") +
  theme_classic()

#plotting the randomly generated points, tree points, and probability/distance raster
ggplot()+ 
  geom_stars(data=na.omit(st_as_stars(dist_near_trails_fixed_box_ILM_inverse), aes(fill = layer)))+ #plotting the distance inverse raster 
  scale_fill_distiller(palette = "Blues", na.value = "transparent", trans = "reverse")+
  geom_sf(data=st_cast(ILM_ANN_Anlysis_trails$random_points$geom, "POINT"), alpha = 0.5, aes(color = "Randomly Generated"), fill = NA, shape = 16) + #plotting the random points
  geom_sf(data=ILM_fixed_field_data_processed_sf, aes(color = "Actual Trees"), shape = 16, alpha = 0.5)+ #plotting the tree points
  labs(color = "Actual Trees", fill = "Inverse Distance (1/m)", 
       x = "Longitude", 
       y = "Latitude")+
  scale_color_manual(
    name = "Trees",
    values = c("Actual Trees" = "red", 
               "Randomly Generated" = "black"))+
  theme_minimal()+
  # guides(color = guide_legend(override.aes = list(shape = c(16,16), linetype = 0)))+
  labs(title = "Isle la Motte")+
  theme_classic() +
  theme(title=element_text(size=15), 
        axis.text=element_text(size=15),  axis.title.x =element_text(size= 15),
        axis.title.y =element_text(size= 15),
        text = element_text(family = "serif"))

#graphing the histogram of simulated ANN values and the mean ANN from our trees
as_tibble(ILM_ANN_Anlysis_trails$ann.r) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "skyblue", color = "black", bins = 50) + 
  xlim(range(ILM_ANN_Anlysis_trails$observed_ANN, ILM_ANN_Anlysis_trails$ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=ILM_ANN_Anlysis$observed_ANN, col = "red") + #plotting our tree's mean ANN
  xlab("Average Nearest Neighbor") +
  ylab("Frequency")+
  labs(title = "Isle la Motte")+
  geom_text(aes(label = round(ILM_ANN_Anlysis_trails$observed_ANN, 2)), x = 6.7, y = 40, color = "red") +
  theme_classic() +
  theme(title=element_text(size=15), 
        axis.text=element_text(size=15),  axis.title.x =element_text(size= 15),
        axis.title.y =element_text(size= 15),
        text = element_text(family = "serif"))

## Version of ANN analysis controlling for the trails with on and inside the trails 

ILM_ANN_Anlysis_on_inside_trails <- ANN_analysis("ILM", "On and Inside trails")
ILM_ANN_Anlysis_on_inside_trails #first index is the ANN value, the second is the left-tailed p-value

#plotting the randomly generated points, tree points, and trails raster
ggplot()+ 
  geom_stars(data=st_rasterize(trails_ILM_trans))+ #plotting the trails raster 
  geom_sf(data=ILM_fixed_field_data_processed_sf, aes(col = "red"))+ #plotting the tree points
  geom_sf(data=ILM_ANN_Anlysis_on_inside_trails$random_points, fill = NA) #plotting the random points

#graphing the histogram of simulated ANN values and the mean ANN from our trees
as_tibble(ILM_ANN_Anlysis_on_inside_trails$ann.r) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "skyblue", color = "black", bins = 50) + 
  xlim(range(ILM_ANN_Anlysis_on_inside_trails$observed_ANN, ILM_ANN_Anlysis_on_inside_trails$ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=ILM_ANN_Anlysis_on_inside_trails$observed_ANN, col = "red") + #plotting our tree's mean ANN
  xlab("Average Nearest Neighbor") +
  ylab("Frequency")+
  theme_classic()


#### PPM analysis ####

# for every PPM analysis 

# 1) generate an image object of distance raster 
# 2) create the Poisson Point Model using ppm() function with the trails influencing the location of the points (Alternative Hypothesis Model)
# 3) create the Poisson Point Model using ppm() function with the trails not influencing the location of the points (Null Hypothesis Model)
# 4) Use an ANOVA likelihood Ratio Test to compare the Alternate and Null hypotheses
# 5) Plot the influence of the trails as the distance to the trails decreases (inverse distance)

#Test for ILM

#creating the image of the distance to trails stars
dist_near_trails_fixed_box_ILM_inverse_im <- as.im(st_as_stars(dist_near_trails_fixed_box_ILM_inverse))

ILM_ppp <- unmark(ILM_ppp)

#Alternative hypothesis, seeing if the distance to the trails's edge influences the tree point placement
PPM1 <- ppm(Q = ILM_ppp ~ dist_near_trails_fixed_box_ILM_inverse_im) 
PPM1

#null hypothesis, no change in the trend of the points
PPM0 <- ppm(ILM_ppp ~ 1)
PPM0

#using a likelihood ratio test to compare the alternative and null models
anova(PPM0, PPM1, test="LRT")

#plotting the alternative model
plot(effectfun(PPM1, "dist_near_trails_fixed_box_ILM_inverse_im", se.fit = TRUE), main = "Distance to trails of ILM",
     ylab = "Butternut Trees", xlab = "Inverse Distance to trails", legend = FALSE)

## Using the non-inverse distance to trails raster ##


#### Creating the distance to trails rasters where everything inside the trails equals 1 ####

#making sure the trails polygon and distance to trails raster have the same projection
trails_vect_ILM <- project(vect(st_as_sf(trails_ILM)), rast(dist_near_trails_fixed_box_ILM))

#creating the corrected distance to trails column where the values inside the cells touching the polygon equal 0
dist_near_trails_fixed_box_ILM_corrected <- rasterize(trails_vect_ILM, rast(dist_near_trails_fixed_box_ILM), field=0, update=TRUE, touches=TRUE)

#making sure the projections are identical for the buffer and the raster
#this for some reason doesn't work
#trails_fixed_box_ILM <- project(trails_buffer_ILM, rast(dist_near_trails_fixed_box_ILM_corrected))

#trimming off the NAs
dist_near_trails_fixed_box_ILM_corrected <- trim(dist_near_trails_fixed_box_ILM_corrected)

ggplot()+
  geom_raster(data = as.data.frame(dist_near_trails_fixed_box_ILM_corrected, xy=T), aes(x=x, y=y, fill = layer))+
  geom_sf(data = trails_ILM_trans)+
  geom_sf(data = ILM_fixed_field_data_processed_sf)

#Test for ILM

#creating the image of the distance to trails stars
dist_near_trails_fixed_box_ILM_corrected_im <- as.im(st_as_stars(dist_near_trails_fixed_box_ILM_corrected))

#creating a poison point model object of our known trees 
ILM_fixed_field_data_processed_ppp <- as.ppp(ILM_fixed_field_data_processed_sf)
ILM_fixed_field_data_processed_ppp <- unmark(ILM_fixed_field_data_processed_ppp)

#Alternative hypothesis, seeing if the distance to the trails's edge influences the tree point placement
PPM1 <- ppm(Q = ILM_fixed_field_data_processed_ppp ~ dist_near_trails_fixed_box_ILM_corrected_im, na.rm = TRUE) 
PPM1

#null hypothesis, no change in the trend of the points
PPM0 <- ppm(ILM_fixed_field_data_processed_ppp ~ 1)
PPM0

#using a likelihood ratio test to compare the alternative and null models
anova(PPM0, PPM1, test="LRT")

#plotting the alternative model
plot(effectfun(PPM1, "dist_near_trails_fixed_box_ILM_corrected_im", se.fit = TRUE), main = "Distance to trails of Isle la Motte",
     ylab = "Butternut Trees", xlab = "Distance to trails", legend = FALSE)

### CPVT ###

# plotting the points, river, and bounding box (encompassing all of the outlying trees)
ggplot(CPVT_fixed_field_data_processed_box)+
  geom_sf() +
  geom_sf(data = trails_CPVT_trans) +
  geom_sf(data = CPVT_fixed_field_data_processed_sf)

# creating one geometry with the tree points using st_union and then creating a convex hull with st_convex_hull
trails_CPVT_convex_hull <- st_convex_hull(st_union(CPVT_fixed_field_data_processed_sf))
#let's check on that
ggplot(trails_CPVT_convex_hull)+
  geom_sf() +
  geom_sf(data = trails_CPVT_trans) +
  geom_sf(data = CPVT_fixed_field_data_processed_sf)

## Calculating the Ripley's K: ##
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

# Ripley's K for CPVT 

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

#unlike Rebecca, I can't have the window set as a buffer around the trail network because those buffers overlap weirdly, so I skip that section

#### ANN Analysis (test for clustering/dispersion) ####

# For all ANN Analyses (for each convex hulls/buffers):
# 1) Find Average Nearest Neighbor Value for the trees of a population
# 2) Simulate a randomized distribution of points, calculate and store the average nearest neighbor, 599 times
# 3) Make sure the randomized points have the correct coordinate reference system for mapping
# 4) Plot the Randomized Points, convex hull, and trail shapefiles
# 5) Plot the histogram of simulated randomized average nearest neighbor values with a line for the actual ANN value
# 6) Calculate and store the p-value 

# The convex hull windows were used because they produced the most close-to-the-original point distances as possible.

# In this analysis, we do not control for the presence of the trails as something that would influence the distribution of points

# When the average nearest neighbor is to the left of the histogram the points are more clustered than expected at random
# When the average nearest neighbor is to the right of the histogram the points are more dispersed than expected at random

# creating the trail rasters for the ANN analysis that will be used for point generation later
#making sure the trails polygon and distance to trails raster have the same projection
trails_vect_CPVT <- project(vect(trails_CPVT), rast(dist_near_trails_fixed_box_CPVT))

#creating the corrected distance to trails column where the values inside the cells touching the polygon equal 0
dist_near_trails_fixed_box_CPVT_corrected <- rasterize(trails_vect_CPVT, rast(dist_near_trails_fixed_box_CPVT), field=0, update=TRUE, touches=TRUE)

#turning the distance to trails correct raster into a stars object
dist_near_trails_fixed_box_CPVT_corrected_stars <- st_as_stars(dist_near_trails_fixed_box_CPVT_corrected)

#creating the inverse of the distance raster so that the higher values are closer to the trails and the values are between 0-1
#this is the problem layer I believe for the ppm, because it excludes most of the values
dist_near_trails_fixed_box_CPVT_inverse <- dist_near_trails_fixed_box_CPVT_corrected_stars %>%
  st_as_sf() %>%
  mutate(layer = case_when(layer >= 1 ~ 1/layer,
                           layer > 0 & layer < 1 ~ 1,
                           layer == 0 ~ 1)) %>%
  st_rasterize()

#creating a raster out of the inverse distance stars object
dist_near_trails_fixed_box_CPVT_inverse_im_raster <- rast(dist_near_trails_fixed_box_CPVT_inverse)

#projecting the inverse distance raster to match the other crs
crs(dist_near_trails_fixed_box_CPVT_inverse_im_raster) <- crs(rast(dist_near_trails_fixed_box_CPVT_inverse))

#creating ANN Analysis function

ANN_analysis <- function(population, window) {
  if (population == "CPVT") {
    ppp <- CPVT_ppp #assigning poisson point pattern 
    dataframe <- CPVT_fixed_field_data_processed_sf #assigning dataframe
    
    #window selection
    if (window == "Convex Hull"){ #ANN without controlling for trails
      selected_window <- trails_CPVT_convex_hull
    } else if (window == "Just trails"){ #ANN with controlling for trails
      selected_window <- trails_CPVT_trans_point_raster
    } else if (window == "On and Inside trails"){
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
    #ANN analysis controlling for trails
    n <- 599L #defines the number of simulations
    ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
    for (i in 1:n){ 
      rand.p <- rpoint(n=length(dataframe), f = as.im(selected_window)) # generating the random points within the window
      ann.r[i] <- mean(nndist(rand.p, k=1)) #for each simulated random distribution of points it calculates the mean ANN across all of the trees
    } #for the length of the number of points at LM, it assigns a random point on top of the trails's edge while controlling for the trails's edge
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

CPVT_ANN_Anlysis <- ANN_analysis("CPVT", "Convex Hull")
CPVT_ANN_Anlysis #first index is the ANN value, the second is the left-tailed p-value

#plotting the randomly generated points, tree points, and the trails
ggplot()+ 
  geom_sf(data=trails_CPVT_trans)+ #plotting the trails
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

#### ANN Analysis (test for clustering/dispersion) while controlling for the trails ####

# The steps for this ANN are the same as previously, except we use three different versions of the windows in which we generate random points
# with varying levels of control for the trails to see if the points still seem significantly clustered despite the 
# presence of the trailss (similar to the PPM analysis later).

# The three ways of controlling for the trails include 
# a) controlling for the trails border (using a trails multipoint raster window), 
# b) controlling for on, inside, and around the trails (using an inverse distance raster window), and 
# c) controlling for on and inside the trails (using a raster of the trails window)


# To do this, we add new steps in the beginning 
# 1) Create rasters of the trails shapefile, trails fixed_box, and create a raster with the inverse distances of 
#each cell to the trails shapefile (closer cells are weighted higher)
#cells within a certain distance of the trails equal 1 and the other points equals 1/distance
# 2) Run the simulations whereby the windows either use the trails border raster, the inverse distance raster where 
#the randomized points are placed more likely based on the raster or the higher cell weights, and the trails polygon raster

## Version of ANN analysis controlling for the trails with just the trails multipoint 
CPVT_ANN_Anlysis_trails <- ANN_analysis("CPVT", "Just trails")
CPVT_ANN_Anlysis_trails #first index is the ANN value, the second is the left-tailed p-value

#plotting the randomly generated points, tree points, and probability/distance raster
ggplot()+ 
  geom_stars(data=trails_CPVT_trans_point_raster)+ #plotting the trails edge raster
  geom_sf(data=CPVT_fixed_field_data_processed_sf, aes(col = "red"))+ #plotting the tree points
  geom_sf(data=CPVT_ANN_Anlysis_trails$random_points$geom, fill = NA) #plotting the random points

#graphing the histogram of simulated ANN values and the mean ANN from our trees
as_tibble(CPVT_ANN_Anlysis_trails$ann.r) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "dodgerblue1", color = "black", bins = 50) + 
  xlim(range(CPVT_ANN_Anlysis_trails$observed_ANN, CPVT_ANN_Anlysis_trails$ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=CPVT_ANN_Anlysis_trails$observed_ANN, col = "red") + #plotting our tree's mean ANN
  xlab("ANN") +
  theme_classic()

#plotting the randomly generated points, tree points, and probability/distance raster
ggplot()+ 
  geom_stars(data=na.omit(st_as_stars(dist_near_trails_fixed_box_CPVT_inverse), aes(fill = layer)))+ #plotting the distance inverse raster 
  scale_fill_distiller(palette = "Blues", na.value = "transparent", trans = "reverse")+
  geom_sf(data=st_cast(CPVT_ANN_Anlysis_trails$random_points$geom, "POINT"), alpha = 0.5, aes(color = "Randomly Generated"), fill = NA, shape = 16) + #plotting the random points
  geom_sf(data=CPVT_fixed_field_data_processed_sf, aes(color = "Actual Trees"), shape = 16, alpha = 0.5)+ #plotting the tree points
  labs(color = "Actual Trees", fill = "Inverse Distance (1/m)", 
       x = "Longitude", 
       y = "Latitude")+
  scale_color_manual(
    name = "Trees",
    values = c("Actual Trees" = "red", 
               "Randomly Generated" = "black"))+
  theme_minimal()+
  # guides(color = guide_legend(override.aes = list(shape = c(16,16), linetype = 0)))+
  labs(title = "Charlotte Park")+
  theme_classic() +
  theme(title=element_text(size=15), 
        axis.text=element_text(size=15),  axis.title.x =element_text(size= 15),
        axis.title.y =element_text(size= 15),
        text = element_text(family = "serif"))

#graphing the histogram of simulated ANN values and the mean ANN from our trees
as_tibble(CPVT_ANN_Anlysis_trails$ann.r) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "skyblue", color = "black", bins = 50) + 
  xlim(range(CPVT_ANN_Anlysis_trails$observed_ANN, CPVT_ANN_Anlysis_trails$ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=CPVT_ANN_Anlysis$observed_ANN, col = "red") + #plotting our tree's mean ANN
  xlab("Average Nearest Neighbor") +
  ylab("Frequency")+
  labs(title = "Charlotte Park")+
  geom_text(aes(label = round(CPVT_ANN_Anlysis_trails$observed_ANN, 2)), x = 6.7, y = 40, color = "red") +
  theme_classic() +
  theme(title=element_text(size=15), 
        axis.text=element_text(size=15),  axis.title.x =element_text(size= 15),
        axis.title.y =element_text(size= 15),
        text = element_text(family = "serif"))

## Version of ANN analysis controlling for the trails with on and inside the trails 

CPVT_ANN_Anlysis_on_inside_trails <- ANN_analysis("CPVT", "On and Inside trails")
CPVT_ANN_Anlysis_on_inside_trails #first index is the ANN value, the second is the left-tailed p-value

#plotting the randomly generated points, tree points, and trails raster
ggplot()+ 
  geom_stars(data=st_rasterize(trails_CPVT_trans))+ #plotting the trails raster 
  geom_sf(data=CPVT_fixed_field_data_processed_sf, aes(col = "red"))+ #plotting the tree points
  geom_sf(data=CPVT_ANN_Anlysis_on_inside_trails$random_points, fill = NA) #plotting the random points

#graphing the histogram of simulated ANN values and the mean ANN from our trees
as_tibble(CPVT_ANN_Anlysis_on_inside_trails$ann.r) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "skyblue", color = "black", bins = 50) + 
  xlim(range(CPVT_ANN_Anlysis_on_inside_trails$observed_ANN, CPVT_ANN_Anlysis_on_inside_trails$ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=CPVT_ANN_Anlysis_on_inside_trails$observed_ANN, col = "red") + #plotting our tree's mean ANN
  xlab("Average Nearest Neighbor") +
  ylab("Frequency")+
  theme_classic()


#### PPM analysis ####

# for every PPM analysis 

# 1) generate an image object of distance raster 
# 2) create the Poisson Point Model using ppm() function with the trails influencing the location of the points (Alternative Hypothesis Model)
# 3) create the Poisson Point Model using ppm() function with the trails not influencing the location of the points (Null Hypothesis Model)
# 4) Use an ANOVA likelihood Ratio Test to compare the Alternate and Null hypotheses
# 5) Plot the influence of the trails as the distance to the trails decreases (inverse distance)

#Test for CPVT

#creating the image of the distance to trails stars
dist_near_trails_fixed_box_CPVT_inverse_im <- as.im(st_as_stars(dist_near_trails_fixed_box_CPVT_inverse))

CPVT_ppp <- unmark(CPVT_ppp)

#Alternative hypothesis, seeing if the distance to the trails's edge influences the tree point placement
PPM1 <- ppm(Q = CPVT_ppp ~ dist_near_trails_fixed_box_CPVT_inverse_im) 
PPM1

#null hypothesis, no change in the trend of the points
PPM0 <- ppm(CPVT_ppp ~ 1)
PPM0

#using a likelihood ratio test to compare the alternative and null models
anova(PPM0, PPM1, test="LRT")

#plotting the alternative model
plot(effectfun(PPM1, "dist_near_trails_fixed_box_CPVT_inverse_im", se.fit = TRUE), main = "Distance to trails of CPVT",
     ylab = "Butternut Trees", xlab = "Inverse Distance to trails", legend = FALSE)

## Using the non-inverse distance to trails raster ##


#### Creating the distance to trails rasters where everything inside the trails equals 1 ####

#making sure the trails polygon and distance to trails raster have the same projection
trails_vect_CPVT <- project(vect(st_as_sf(trails_CPVT)), rast(dist_near_trails_fixed_box_CPVT))

#creating the corrected distance to trails column where the values inside the cells touching the polygon equal 0
dist_near_trails_fixed_box_CPVT_corrected <- rasterize(trails_vect_CPVT, rast(dist_near_trails_fixed_box_CPVT), field=0, update=TRUE, touches=TRUE)

#making sure the projections are identical for the buffer and the raster
#this for some reason doesn't work
#trails_fixed_box_CPVT <- project(trails_buffer_CPVT, rast(dist_near_trails_fixed_box_CPVT_corrected))

#trimming off the NAs
dist_near_trails_fixed_box_CPVT_corrected <- trim(dist_near_trails_fixed_box_CPVT_corrected)

ggplot()+
  geom_raster(data = as.data.frame(dist_near_trails_fixed_box_CPVT_corrected, xy=T), aes(x=x, y=y, fill = layer))+
  geom_sf(data = trails_CPVT_trans)+
  geom_sf(data = CPVT_fixed_field_data_processed_sf)

#Test for CPVT

#creating the image of the distance to trails stars
dist_near_trails_fixed_box_CPVT_corrected_im <- as.im(st_as_stars(dist_near_trails_fixed_box_CPVT_corrected))

#creating a poison point model object of our known trees 
CPVT_fixed_field_data_processed_ppp <- as.ppp(CPVT_fixed_field_data_processed_sf)
CPVT_fixed_field_data_processed_ppp <- unmark(CPVT_fixed_field_data_processed_ppp)

#Alternative hypothesis, seeing if the distance to the trails's edge influences the tree point placement
PPM1 <- ppm(Q = CPVT_fixed_field_data_processed_ppp ~ dist_near_trails_fixed_box_CPVT_corrected_im, na.rm = TRUE) 
PPM1

#null hypothesis, no change in the trend of the points
PPM0 <- ppm(CPVT_fixed_field_data_processed_ppp ~ 1)
PPM0

#using a likelihood ratio test to compare the alternative and null models
anova(PPM0, PPM1, test="LRT")

#plotting the alternative model
plot(effectfun(PPM1, "dist_near_trails_fixed_box_CPVT_corrected_im", se.fit = TRUE), main = "Distance to trails of Charlotte Park",
     ylab = "Butternut Trees", xlab = "Distance to trails", legend = FALSE)