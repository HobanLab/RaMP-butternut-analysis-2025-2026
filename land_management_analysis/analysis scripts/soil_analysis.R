#Catherine dell'Olio
#Soil analysis

#### Loading libraries and relevant data ####

library(tidyverse)
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc) #ggplot extension
library(PMCMRplus) # for Dunn test
library(geomtextpath) # for PCA graphing
library(spatstat) # to run the nndist function
library(raster) #to plot rasters
library(rstatix) #to run the Games-Howell Test
library(ggnewscale) #to be able to assign different colors to different layered rasters

#soil data brought to you by https://websoilsurvey.nrcs.usda.gov/app/WebSoilSurvey.aspx
setwd("~/GitHub/butternut-health-assessment-2025/land_management_analysis/data/wss_aoi_2026_ILM/spatial")
soil_map_ILM <- st_read("soilmu_a_aoi.shp", crs = 4326)
soil_map_ILM_transformed <- st_crop(st_transform(soil_map_ILM, crs = 26918), ILM_fixed_field_data_processed_box)
plot(soil_map_ILM_transformed)
#note that the WSS system labels the soil categories as "MUSYM" and so that column is what denotes soil category for each butternut, determined by where that butternut's point geometry intersects with the soil polygons
ILM_data_w_soils <- st_join(ILM_fixed_field_data_processed_sf_transformed, soil_map_ILM_transformed, join = st_intersects)


#RANDOM POINT ANALYSIS

set.seed(25)

#creating the dataframe that stores how many points are in each soil type for each random distribution
stored_random_distributions <- data.frame(MUSYM = unique(soil_map_ILM$MUSYM))
#this loop will run 1000 times (so it takes a while to run!!) and do a few things:
#1) create a random 20 point distribution within the bounding box of ILM points
#2) figure out what soil type each random point is on
#3) count up the hits to each soil type
#4) store those hits on the above dataframe
for (y in 1:1000){ 
  #generating and extracting the randomly distributed population soil categories
  random_249 <- st_sample(ILM_fixed_field_data_processed_box, 249) #select random 20 points within the cropped ILM polygon
  random_249 <- random_249 %>%
    st_as_sf() %>% #making sure the random points are stored as simple features
    st_transform(random_249, crs = 26918) #making sure they are in the right CRS
  
  random_249_soil_categories <- st_join(random_249, soil_map_ILM_transformed, join = st_intersects) #linking each random point to its soil category
  
  random_249_soil_counts <- as.data.frame(random_249_soil_categories %>% count(MUSYM)) #counting up how many random hits we got to each soil category
  random_249_soil_counts <- subset(random_249_soil_counts, select = -x) #removing the geometry column

  stored_random_distributions <- stored_random_distributions %>% left_join(random_249_soil_counts, by = join_by(MUSYM)) #adding the number of points in each soil category to a storing dataframe (for each iteration of random points)
  stored_random_distributions[is.na(stored_random_distributions)] <- 0 #if no points were found in that soil type, it's listed as 0 points, not NA
}

stored_random_distributions_cleaned <- data.frame(t(stored_random_distributions))
colnames(stored_random_distributions_cleaned) <- stored_random_distributions_cleaned[1,]
stored_random_distributions_cleaned <- stored_random_distributions_cleaned[-1,]
rownames(stored_random_distributions_cleaned) <- NULL
stored_random_distributions_cleaned <- mutate_all(stored_random_distributions_cleaned, function(x) as.numeric((x)))

#now let's get the data on the real stuff
#creating the dataframe that stores how many points are in each soil type for each random distribution
stored_real_distribution <- data.frame(MUSYM = unique(soil_map_ILM$MUSYM))

real_249_soil_counts <- as.data.frame(ILM_data_w_soils %>% count(MUSYM))

real_249_soil_counts <- subset(real_249_soil_counts, select = -geometry) #removing the geometry column
#why the flippity flop is it called geometry here and x above????
  
stored_real_distribution <- stored_real_distribution %>% left_join(real_249_soil_counts, by = join_by(MUSYM)) #adding the number of points in each soil category to a storing dataframe (for each iteration of random points)

stored_real_distribution[is.na(stored_real_distribution)] <- 0 #if no points were found in that soil type, it's listed as 0 points, not NA

stored_real_distribution_cleaned <- data.frame(t(stored_real_distribution))
colnames(stored_real_distribution_cleaned) <- stored_real_distribution_cleaned[1,]
stored_real_distribution_cleaned <- stored_real_distribution_cleaned[-1,]
rownames(stored_real_distribution_cleaned) <- NULL
# Source - https://stackoverflow.com/a/48419207
# Posted by R. Prost, modified by community. See post 'Timeline' for change history
# Retrieved 2026-08-05, License - CC BY-SA 4.0
stored_real_distribution_cleaned <- mutate_all(stored_real_distribution_cleaned, function(x) as.numeric((x)))


#YOU HAVE TO RUN THIS AS ONE BIG BLOCK
real_vs_random_quantiles <- data.frame(MUSYM = unique(soil_map_ILM$MUSYM))
real_vs_random_quantiles <- data.frame(t(real_vs_random_quantiles))
colnames(real_vs_random_quantiles) <- real_vs_random_quantiles[1,]
for (i in 1:ncol(stored_random_distributions_cleaned)){
  cat <- real_vs_random_quantiles[[i]]
  real_vs_random_quantiles[1,i] <- ecdf(stored_random_distributions_cleaned[[cat]])(as.numeric(stored_real_distribution_cleaned[[cat]]))
}
#END OF BLOCK

#plotting the histogram of the randomly distributed means and our real mean
plot_out_histogram <- function(cat){
  # Source - https://stackoverflow.com/a/36015931
  ggplot(stored_random_distributions_cleaned, aes(x = {{cat}}))+
    geom_histogram(fill = "dodgerblue1", color = "black")+
    geom_vline(xintercept= as.numeric(stored_real_distribution_cleaned[[deparse(substitute(cat))]]), col = "red")+ #line of our real slope
    labs(title = paste0("Simulated Number of ILM butternuts in Soil ", deparse(substitute(cat)), " (n=249)"), subtitle = paste0("Quantile of real value: ", real_vs_random_quantiles[[deparse(substitute(cat))]]))+
    theme_classic()
}

plot_out_histogram(BcB)
plot_out_histogram(AaC)
plot_out_histogram(BdC)
plot_out_histogram(EaA)
plot_out_histogram(AaA)
plot_out_histogram(AbA)
plot_out_histogram(CbA)
plot_out_histogram(KbA)
plot_out_histogram(AaB)
plot_out_histogram(BdB)
plot_out_histogram(CdB)
plot_out_histogram(BdD)
plot_out_histogram(KaC)
plot_out_histogram(BcC)
plot_out_histogram(SdA)
plot_out_histogram(KaB)
plot_out_histogram(MaB)
plot_out_histogram(BdE)
plot_out_histogram(KbB)
plot_out_histogram(KcA)
plot_out_histogram(AbC)
#W is water; no butternuts
plot_out_histogram(AbB)


#----
#CPVT
#soil data brought to you by https://websoilsurvey.nrcs.usda.gov/app/WebSoilSurvey.aspx
setwd("~/GitHub/butternut-health-assessment-2025/land_management_analysis/data/wss_aoi_2026_CPVT/spatial")
soil_map_CPVT <- st_read("soilmu_a_aoi.shp", crs = 4326)
soil_map_CPVT_transformed <- st_crop(st_transform(soil_map_CPVT, crs = 26918), CPVT_fixed_field_data_processed_box)
plot(soil_map_CPVT_transformed)
#note that the WSS system labels the soil categories as "MUSYM" and so that column is what denotes soil category for each butternut, determined by where that butternut's point geometry intersects with the soil polygons
CPVT_data_w_soils <- st_join(CPVT_fixed_field_data_processed_sf_transformed, soil_map_CPVT_transformed, join = st_intersects)


#RANDOM POINT ANALYSIS

set.seed(25)

#creating the dataframe that stores how many points are in each soil type for each random distribution
stored_random_distributions <- data.frame(MUSYM = unique(soil_map_CPVT$MUSYM))
#this loop will run 1000 times (so it takes a while to run!!) and do a few things:
#1) create a random 20 point distribution within the bounding box of CPVT points
#2) figure out what soil type each random point is on
#3) count up the hits to each soil type
#4) store those hits on the above dataframe
for (y in 1:1000){ 
  #generating and extracting the randomly distributed population soil categories
  random_163 <- st_sample(CPVT_fixed_field_data_processed_box, 163) #select random 20 points within the cropped CPVT polygon
  random_163 <- random_163 %>%
    st_as_sf() %>% #making sure the random points are stored as simple features
    st_transform(random_163, crs = 26918) #making sure they are in the right CRS
  
  random_163_soil_categories <- st_join(random_163, soil_map_CPVT_transformed, join = st_intersects) #linking each random point to its soil category
  
  random_163_soil_counts <- as.data.frame(random_163_soil_categories %>% count(MUSYM)) #counting up how many random hits we got to each soil category
  random_163_soil_counts <- subset(random_163_soil_counts, select = -x) #removing the geometry column
  
  stored_random_distributions <- stored_random_distributions %>% left_join(random_163_soil_counts, by = join_by(MUSYM)) #adding the number of points in each soil category to a storing dataframe (for each iteration of random points)
  stored_random_distributions[is.na(stored_random_distributions)] <- 0 #if no points were found in that soil type, it's listed as 0 points, not NA
}

stored_random_distributions_cleaned <- data.frame(t(stored_random_distributions))
colnames(stored_random_distributions_cleaned) <- stored_random_distributions_cleaned[1,]
stored_random_distributions_cleaned <- stored_random_distributions_cleaned[-1,]
rownames(stored_random_distributions_cleaned) <- NULL
stored_random_distributions_cleaned <- mutate_all(stored_random_distributions_cleaned, function(x) as.numeric((x)))

#now let's get the data on the real stuff
#creating the dataframe that stores how many points are in each soil type for each random distribution
stored_real_distribution <- data.frame(MUSYM = unique(soil_map_CPVT$MUSYM))

real_163_soil_counts <- as.data.frame(CPVT_data_w_soils %>% count(MUSYM))

real_163_soil_counts <- subset(real_163_soil_counts, select = -geometry) #removing the geometry column
#why the flippity flop is it called geometry here and x above????

stored_real_distribution <- stored_real_distribution %>% left_join(real_163_soil_counts, by = join_by(MUSYM)) #adding the number of points in each soil category to a storing dataframe (for each iteration of random points)

stored_real_distribution[is.na(stored_real_distribution)] <- 0 #if no points were found in that soil type, it's listed as 0 points, not NA

stored_real_distribution_cleaned <- data.frame(t(stored_real_distribution))
colnames(stored_real_distribution_cleaned) <- stored_real_distribution_cleaned[1,]
stored_real_distribution_cleaned <- stored_real_distribution_cleaned[-1,]
rownames(stored_real_distribution_cleaned) <- NULL
# Source - https://stackoverflow.com/a/48419207
# Posted by R. Prost, modified by community. See post 'Timeline' for change history
# Retrieved 2026-08-05, License - CC BY-SA 4.0
stored_real_distribution_cleaned <- mutate_all(stored_real_distribution_cleaned, function(x) as.numeric((x)))


#YOU HAVE TO RUN THIS AS ONE BIG BLOCK
real_vs_random_quantiles <- data.frame(MUSYM = unique(soil_map_CPVT$MUSYM))
real_vs_random_quantiles <- data.frame(t(real_vs_random_quantiles))
colnames(real_vs_random_quantiles) <- real_vs_random_quantiles[1,]
for (i in 1:ncol(stored_random_distributions_cleaned)){
  cat <- real_vs_random_quantiles[[i]]
  real_vs_random_quantiles[1,i] <- ecdf(stored_random_distributions_cleaned[[cat]])(as.numeric(stored_real_distribution_cleaned[[cat]]))
}
#END OF BLOCK

#plotting the histogram of the randomly distributed means and our real mean
plot_out_histogram <- function(cat){
  # Source - https://stackoverflow.com/a/36015931
  ggplot(stored_random_distributions_cleaned, aes(x = {{cat}}))+
    geom_histogram(fill = "dodgerblue1", color = "black")+
    geom_vline(xintercept= as.numeric(stored_real_distribution_cleaned[[deparse(substitute(cat))]]), col = "red")+ #line of our real slope
    labs(title = paste0("Simulated Number of CPVT butternuts in Soil ", deparse(substitute(cat)), " (n=163)"), subtitle = paste0("Quantile of real value: ", real_vs_random_quantiles[[deparse(substitute(cat))]]))+
    theme_classic()
}

plot_out_histogram(SuD)
plot_out_histogram(VeC)
plot_out_histogram(AdD)
plot_out_histogram(PaB)
plot_out_histogram(VeB)
plot_out_histogram(SuB)
plot_out_histogram(SxE)
plot_out_histogram(BlA)
plot_out_histogram(EwA)
plot_out_histogram(GeB)
plot_out_histogram(HnC)
plot_out_histogram(Cv)
plot_out_histogram(FaC)
plot_out_histogram(PaD)
plot_out_histogram(PaC)
plot_out_histogram(SxC)
plot_out_histogram(Lh)
plot_out_histogram(GgC)
plot_out_histogram(SuC)
plot_out_histogram(BlB)
plot_out_histogram(VeD)
#W is water; no butternuts
plot_out_histogram(HnA)
plot_out_histogram(HnB)
plot_out_histogram(FsB)
plot_out_histogram(Le)
plot_out_histogram(MnC)
plot_out_histogram(FaE)
plot_out_histogram(FsE)


#------------
#HEALTH PLOTS
#these next few split up our data by age, and also select all of the living trees (which have more data to their name)
data_ILMadults_w_soil <- CPVT_data_w_soils %>% filter(seedling_or_adult == "Adult")
data_ILMliveadults_w_soil <- data_ILMadults_w_soil %>% filter(adult_dead_or_alive == "Alive")

data_ILMseedlings_w_soil <- ILM_data_w_soils %>% filter(seedling_or_adult == "Seedling")
data_ILMliveseedlings_w_soil <- data_ILMseedlings_w_soil %>% filter(seedling_dead_or_alive == "Alive")

data_CPVTadults_w_soil <- ILM_data_w_soils %>% filter(seedling_or_adult == "Adult")
data_CPVTliveadults_w_soil <- data_CPVTadults_w_soil %>% filter(adult_dead_or_alive == "Alive")

data_CPVTseedlings_w_soil <- CPVT_data_w_soils %>% filter(seedling_or_adult == "Seedling")
data_CPVTliveseedlings_w_soil <- data_CPVTseedlings_w_soil %>% filter(seedling_dead_or_alive == "Alive")

#these plots will show any trends in health status across the two sites based on soil type
barplot(data_CPVT_w_soil, MUSYM, "Soil at CPVT")
barplot(data_CPVTseedlings_w_soil, MUSYM, "Soil at CPVT")

boxplot(data_ILMliveadults_w_soil, MUSYM, adult_percent_live_canopy, "Canopy at ILM by soil")
boxplot(data_CPVTliveadults_w_soil, MUSYM, adult_percent_live_canopy, "Canopy at CPVT by soil")

boxplot(data_ILMliveadults_w_soil, MUSYM, live_adult_girdle, "Girdle at ILM by soil")
boxplot(data_CPVTliveadults_w_soil, MUSYM, live_adult_girdle, "Girdle at CPVT by soil")

boxplot(data_ILMliveadults_w_soil, MUSYM, live_adult_purdue_canker_rating, "Purdue canker ranking at ILM by soil")
boxplot(data_CPVTliveadults_w_soil, MUSYM, live_adult_purdue_canker_rating, "Purdue canker ranking at CPVT by soil")

boxplot(data_ILMliveadults_w_soil, MUSYM, live_adult_purdue_canopy_ranking, "Purdue canopy ranking at ILM by soil")
boxplot(data_CPVTliveadults_w_soil, MUSYM, live_adult_purdue_canopy_ranking, "Purdue canopy ranking at CPVT by soil")

#overall soil maps
map(data_ILM, MUSYM, TRUE, "soils at ILM")
map(data_CPVT, MUSYM, TRUE, "soils at CPVT")