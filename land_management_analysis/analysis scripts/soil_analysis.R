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
library(ggplot2) #to be able to make plots and save them
library(EnvStats) #to show sample size for boxplots and such!
library(tmap) #to make maps

working_directory <- "~/GitHub/butternut-health-assessment-2025/land_management_analysis"
setwd(working_directory)

#from my code ChartBuilder.R, a function that makes easy maps
map <- function(data, colorz, legend_yn, title)
  ggplot(data, aes(gps_w, gps_n)) +
  geom_point(aes(color = {{ colorz }}), show.legend = legend_yn) +
  ggtitle(title)

#and boxplots and barplots!!
boxplot <- function(data, plot_division, metric, title){
  ggplot(data, aes(x=factor({{plot_division}}), y={{metric}})) +
    geom_boxplot(alpha = 0, aes(color = {{plot_division}})) +
    geom_jitter(height = 0, alpha = 0.25, size = 0.5, width = 0.3) +
    #stat_n_text() +
    ggtitle(title)
}

barplot <- function(data, variable, title){
  ggplot(data, aes({{variable}})) +
    geom_bar()+
    #stat_n_text() +
    ggtitle(title)
}

#soil data brought to you by https://websoilsurvey.nrcs.usda.gov/app/WebSoilSurvey.aspx
setwd(working_directory)
setwd("data/wss_aoi_2026_ILM/spatial")
soil_map_ILM <- st_read("soilmu_a_aoi.shp", crs = 4326)
soil_map_ILM_transformed <- st_crop(st_transform(soil_map_ILM, crs = 26918), ILM_fixed_field_data_processed_box)
plot(soil_map_ILM_transformed)
#note that the WSS system labels the soil categories as "MUSYM" and so that column is what denotes soil category for each butternut, determined by where that butternut's point geometry intersects with the soil polygons
ILM_data_w_soils <- st_join(ILM_fixed_field_data_processed_sf_transformed, soil_map_ILM_transformed, join = st_intersects)


#RANDOM POINT ANALYSIS

set.seed(25)

#if you want to save time, you can upload the stored random distribution from here, because it takes a while to run
setwd(working_directory)
setwd("data")
load("ILM_stored_random_distributions_cleaned.Rda")
#if you want to run this by yourself:
# #creating the dataframe that stores how many points are in each soil type for each random distribution
# ILM_stored_random_distributions <- data.frame(MUSYM = unique(soil_map_ILM$MUSYM))
# #this loop will run 1000 times (so it takes a while to run!!) and do a few things:
# #1) create a random 20 point distribution within the bounding box of ILM points
# #2) figure out what soil type each random point is on
# #3) count up the hits to each soil type
# #4) store those hits on the above dataframe
# for (y in 1:1000){ 
#   #generating and extracting the randomly distributed population soil categories
#   random_249 <- st_sample(ILM_fixed_field_data_processed_box, 249) #select random 20 points within the cropped ILM polygon
#   random_249 <- random_249 %>%
#     st_as_sf() %>% #making sure the random points are stored as simple features
#     st_transform(random_249, crs = 26918) #making sure they are in the right CRS
#   
#   random_249_soil_categories <- st_join(random_249, soil_map_ILM_transformed, join = st_intersects) #linking each random point to its soil category
#   
#   random_249_soil_counts <- as.data.frame(random_249_soil_categories %>% count(MUSYM)) #counting up how many random hits we got to each soil category
#   random_249_soil_counts <- subset(random_249_soil_counts, select = -x) #removing the geometry column
# 
#   ILM_stored_random_distributions <- ILM_stored_random_distributions %>% left_join(random_249_soil_counts, by = join_by(MUSYM)) #adding the number of points in each soil category to a storing dataframe (for each iteration of random points)
#   ILM_stored_random_distributions[is.na(ILM_stored_random_distributions)] <- 0 #if no points were found in that soil type, it's listed as 0 points, not NA
# }
# 
# ILM_stored_random_distributions_cleaned <- data.frame(t(ILM_stored_random_distributions))
# colnames(ILM_stored_random_distributions_cleaned) <- ILM_stored_random_distributions_cleaned[1,]
# ILM_stored_random_distributions_cleaned <- ILM_stored_random_distributions_cleaned[-1,]
# rownames(ILM_stored_random_distributions_cleaned) <- NULL
# ILM_stored_random_distributions_cleaned <- mutate_all(ILM_stored_random_distributions_cleaned, function(x) as.numeric((x)))

#now let's get the data on the real stuff
#creating the dataframe that stores how many points are in each soil type for each random distribution
ILM_stored_real_distribution <- data.frame(MUSYM = unique(soil_map_ILM$MUSYM))

real_249_soil_counts <- as.data.frame(ILM_data_w_soils %>% count(MUSYM))

real_249_soil_counts <- subset(real_249_soil_counts, select = -geometry) #removing the geometry column
  
ILM_stored_real_distribution <- ILM_stored_real_distribution %>% left_join(real_249_soil_counts, by = join_by(MUSYM)) #adding the number of points in each soil category to a storing dataframe (for each iteration of random points)

ILM_stored_real_distribution[is.na(ILM_stored_real_distribution)] <- 0 #if no points were found in that soil type, it's listed as 0 points, not NA

ILM_stored_real_distribution_cleaned <- data.frame(t(ILM_stored_real_distribution))
colnames(ILM_stored_real_distribution_cleaned) <- ILM_stored_real_distribution_cleaned[1,]
ILM_stored_real_distribution_cleaned <- ILM_stored_real_distribution_cleaned[-1,]
rownames(ILM_stored_real_distribution_cleaned) <- NULL
# Source - https://stackoverflow.com/a/48419207
# Posted by R. Prost, modified by community. See post 'Timeline' for change history
# Retrieved 2026-08-05, License - CC BY-SA 4.0
ILM_stored_real_distribution_cleaned <- mutate_all(ILM_stored_real_distribution_cleaned, function(x) as.numeric((x)))

#YOU HAVE TO RUN THIS AS ONE BIG BLOCK
ILM_real_vs_random_quantiles <- data.frame(MUSYM = unique(soil_map_ILM$MUSYM))
ILM_real_vs_random_quantiles <- data.frame(t(ILM_real_vs_random_quantiles))
colnames(ILM_real_vs_random_quantiles) <- ILM_real_vs_random_quantiles[1,]
for (i in 1:ncol(ILM_stored_random_distributions_cleaned)){
  cat <- ILM_real_vs_random_quantiles[[i]]
  ILM_real_vs_random_quantiles[1,i] <- ecdf(ILM_stored_random_distributions_cleaned[[cat]])(as.numeric(ILM_stored_real_distribution_cleaned[[cat]]))
}
#END OF BLOCK

#plotting the histogram of the randomly distributed means and our real mean
plot_out_histogram <- function(cat){
  # Source - https://stackoverflow.com/a/36015931
  ggplot(ILM_stored_random_distributions_cleaned, aes(x = {{cat}}))+
    geom_histogram(fill = "dodgerblue1", color = "black")+
    geom_vline(xintercept= as.numeric(ILM_stored_real_distribution_cleaned[[deparse(substitute(cat))]]), col = "red")+ #line of our real slope
    labs(title = paste0("Simulated Number of ILM butternuts in Soil ", deparse(substitute(cat)), " (n=249)"), subtitle = paste0("Quantile of real value: ", ILM_real_vs_random_quantiles[[deparse(substitute(cat))]]))+
    theme_classic()
}

setwd(working_directory)
setwd("saved plots/soil_distribution_plots")
#credit to https://www.geeksforgeeks.org/r-language/save-plot-in-data-object-in-base-r/
ggsave("BcB_ILM.png", plot = plot_out_histogram(BcB), width = 6, height = 4, units = "in")
ggsave("AaC_ILM.png", plot = plot_out_histogram(AaC), width = 6, height = 4, units = "in")
ggsave("BdC_ILM.png", plot = plot_out_histogram(BdC), width = 6, height = 4, units = "in")
ggsave("EaA_ILM.png", plot = plot_out_histogram(EaA), width = 6, height = 4, units = "in")
ggsave("AaA_ILM.png", plot = plot_out_histogram(AaA), width = 6, height = 4, units = "in")
ggsave("AbA_ILM.png", plot = plot_out_histogram(AbA), width = 6, height = 4, units = "in")
ggsave("CbA_ILM.png", plot = plot_out_histogram(CbA), width = 6, height = 4, units = "in")
ggsave("KbA_ILM.png", plot = plot_out_histogram(KbA), width = 6, height = 4, units = "in")
ggsave("AaB_ILM.png", plot = plot_out_histogram(AaB), width = 6, height = 4, units = "in")
ggsave("BdB_ILM.png", plot = plot_out_histogram(BdB), width = 6, height = 4, units = "in")
ggsave("CdB_ILM.png", plot = plot_out_histogram(CbB), width = 6, height = 4, units = "in")
ggsave("BdD_ILM.png", plot = plot_out_histogram(BdD), width = 6, height = 4, units = "in")
ggsave("KaC_ILM.png", plot = plot_out_histogram(KaC), width = 6, height = 4, units = "in")
ggsave("BcC_ILM.png", plot = plot_out_histogram(BcC), width = 6, height = 4, units = "in")
ggsave("SdA_ILM.png", plot = plot_out_histogram(SdA), width = 6, height = 4, units = "in")
ggsave("KaB_ILM.png", plot = plot_out_histogram(KaB), width = 6, height = 4, units = "in")
ggsave("MaB_ILM.png", plot = plot_out_histogram(MaB), width = 6, height = 4, units = "in")
ggsave("BdE_ILM.png", plot = plot_out_histogram(BdE), width = 6, height = 4, units = "in")
ggsave("KbB_ILM.png", plot = plot_out_histogram(KbB), width = 6, height = 4, units = "in")
ggsave("KcA_ILM.png", plot = plot_out_histogram(KcA), width = 6, height = 4, units = "in")
ggsave("AbC_ILM.png", plot = plot_out_histogram(AbC), width = 6, height = 4, units = "in")
#W is water; no butternuts
ggsave("AbB_ILM.png", plot = plot_out_histogram(AbB), width = 6, height = 4, units = "in")


#----
#CPVT
#soil data brought to you by https://websoilsurvey.nrcs.usda.gov/app/WebSoilSurvey.aspx
setwd(working_directory)
setwd("data/wss_aoi_2026_CPVT/spatial")
soil_map_CPVT <- st_read("soilmu_a_aoi.shp", crs = 4326)
soil_map_CPVT_transformed <- st_crop(st_transform(soil_map_CPVT, crs = 26918), CPVT_fixed_field_data_processed_box)
plot(soil_map_CPVT_transformed)
#note that the WSS system labels the soil categories as "MUSYM" and so that column is what denotes soil category for each butternut, determined by where that butternut's point geometry intersects with the soil polygons
CPVT_data_w_soils <- st_join(CPVT_fixed_field_data_processed_sf_transformed, soil_map_CPVT_transformed, join = st_intersects)


#RANDOM POINT ANALYSIS

#I recommend just loading in this data as it takes a while to load
setwd(working_directory)
setwd("data")
load("CPVT_stored_random_distributions_cleaned.Rda")

# #If you want to rerun the dataset
# set.seed(25)
# 
# #creating the dataframe that stores how many points are in each soil type for each random distribution
# CPVT_stored_random_distributions <- data.frame(MUSYM = unique(soil_map_CPVT$MUSYM))
# #this loop will run 1000 times (so it takes a while to run!!) and do a few things:
# #1) create a random 20 point distribution within the bounding box of CPVT points
# #2) figure out what soil type each random point is on
# #3) count up the hits to each soil type
# #4) store those hits on the above dataframe
# for (y in 1:1000){
#   #generating and extracting the randomly distributed population soil categories
#   random_163 <- st_sample(CPVT_fixed_field_data_processed_box, 163) #select random 20 points within the cropped CPVT polygon
#   random_163 <- random_163 %>%
#     st_as_sf() %>% #making sure the random points are stored as simple features
#     st_transform(random_163, crs = 26918) #making sure they are in the right CRS
# 
#   random_163_soil_categories <- st_join(random_163, soil_map_CPVT_transformed, join = st_intersects) #linking each random point to its soil category
# 
#   random_163_soil_counts <- as.data.frame(random_163_soil_categories %>% count(MUSYM)) #counting up how many random hits we got to each soil category
#   random_163_soil_counts <- subset(random_163_soil_counts, select = -x) #removing the geometry column
# 
#   CPVT_stored_random_distributions <- CPVT_stored_random_distributions %>% left_join(random_163_soil_counts, by = join_by(MUSYM)) #adding the number of points in each soil category to a storing dataframe (for each iteration of random points)
#   CPVT_stored_random_distributions[is.na(CPVT_stored_random_distributions)] <- 0 #if no points were found in that soil type, it's listed as 0 points, not NA
# }
# 
# CPVT_stored_random_distributions_cleaned <- data.frame(t(CPVT_stored_random_distributions))
# colnames(CPVT_stored_random_distributions_cleaned) <- CPVT_stored_random_distributions_cleaned[1,]
# CPVT_stored_random_distributions_cleaned <- CPVT_stored_random_distributions_cleaned[-1,]
# rownames(CPVT_stored_random_distributions_cleaned) <- NULL
# CPVT_stored_random_distributions_cleaned <- mutate_all(CPVT_stored_random_distributions_cleaned, function(x) as.numeric((x)))

#now let's get the data on the real stuff
#creating the dataframe that stores how many points are in each soil type for each random distribution
CPVT_stored_real_distribution <- data.frame(MUSYM = unique(soil_map_CPVT$MUSYM))

real_163_soil_counts <- as.data.frame(CPVT_data_w_soils %>% count(MUSYM))

real_163_soil_counts <- subset(real_163_soil_counts, select = -geometry) #removing the geometry column
#why the flippity flop is it called geometry here and x above????

CPVT_stored_real_distribution <- CPVT_stored_real_distribution %>% left_join(real_163_soil_counts, by = join_by(MUSYM)) #adding the number of points in each soil category to a storing dataframe (for each iteration of random points)

CPVT_stored_real_distribution[is.na(CPVT_stored_real_distribution)] <- 0 #if no points were found in that soil type, it's listed as 0 points, not NA

CPVT_stored_real_distribution_cleaned <- data.frame(t(CPVT_stored_real_distribution))
colnames(CPVT_stored_real_distribution_cleaned) <- CPVT_stored_real_distribution_cleaned[1,]
CPVT_stored_real_distribution_cleaned <- CPVT_stored_real_distribution_cleaned[-1,]
rownames(CPVT_stored_real_distribution_cleaned) <- NULL
# Source - https://stackoverflow.com/a/48419207
# Posted by R. Prost, modified by community. See post 'Timeline' for change history
# Retrieved 2026-08-05, License - CC BY-SA 4.0
CPVT_stored_real_distribution_cleaned <- mutate_all(CPVT_stored_real_distribution_cleaned, function(x) as.numeric((x)))


#YOU HAVE TO RUN THIS AS ONE BIG BLOCK
CPVT_real_vs_random_quantiles <- data.frame(MUSYM = unique(soil_map_CPVT$MUSYM))
CPVT_real_vs_random_quantiles <- data.frame(t(CPVT_real_vs_random_quantiles))
colnames(CPVT_real_vs_random_quantiles) <- CPVT_real_vs_random_quantiles[1,]
for (i in 1:ncol(CPVT_stored_random_distributions_cleaned)){
  cat <- CPVT_real_vs_random_quantiles[[i]]
  CPVT_real_vs_random_quantiles[1,i] <- ecdf(CPVT_stored_random_distributions_cleaned[[cat]])(as.numeric(CPVT_stored_real_distribution_cleaned[[cat]]))
}
#END OF BLOCK

#plotting the histogram of the randomly distributed means and our real mean
plot_out_histogram <- function(cat){
  # Source - https://stackoverflow.com/a/36015931
  ggplot(CPVT_stored_random_distributions_cleaned, aes(x = {{cat}}))+
    geom_histogram(fill = "dodgerblue1", color = "black")+
    geom_vline(xintercept= as.numeric(CPVT_stored_real_distribution_cleaned[[deparse(substitute(cat))]]), col = "red")+ #line of our real slope
    labs(title = paste0("Simulated Number of CPVT butternuts in Soil ", deparse(substitute(cat)), " (n=163)"), subtitle = paste0("Quantile of real value: ", CPVT_real_vs_random_quantiles[[deparse(substitute(cat))]]))+
    theme_classic()
}

setwd(working_directory)
setwd("saved plots/soil_distribution_plots")
#credit to https://www.geeksforgeeks.org/r-language/save-plot-in-data-object-in-base-r/
ggsave("SuD_CPVT.png", plot = plot_out_histogram(SuD), width = 6, height = 4, units = "in")
ggsave("VeC_CPVT.png", plot = plot_out_histogram(VeC), width = 6, height = 4, units = "in")
ggsave("AdD_CPVT.png", plot = plot_out_histogram(AdD), width = 6, height = 4, units = "in")
ggsave("PaB_CPVT.png", plot = plot_out_histogram(PaB), width = 6, height = 4, units = "in")
ggsave("VeB_CPVT.png", plot = plot_out_histogram(VeB), width = 6, height = 4, units = "in")
ggsave("SuB_CPVT.png", plot = plot_out_histogram(SuB), width = 6, height = 4, units = "in")
ggsave("SxE_CPVT.png", plot = plot_out_histogram(SxE), width = 6, height = 4, units = "in")
ggsave("BlA_CPVT.png", plot = plot_out_histogram(BlA), width = 6, height = 4, units = "in")
ggsave("EwA_CPVT.png", plot = plot_out_histogram(EwA), width = 6, height = 4, units = "in")
ggsave("GeB_CPVT.png", plot = plot_out_histogram(GeB), width = 6, height = 4, units = "in")
ggsave("HnC_CPVT.png", plot = plot_out_histogram(HnC), width = 6, height = 4, units = "in")
ggsave("Cv_CPVT.png", plot = plot_out_histogram(Cv), width = 6, height = 4, units = "in")
ggsave("FaC_CPVT.png", plot = plot_out_histogram(FaC), width = 6, height = 4, units = "in")
ggsave("PaD_CPVT.png", plot = plot_out_histogram(PaD), width = 6, height = 4, units = "in")
ggsave("PaC_CPVT.png", plot = plot_out_histogram(PaC), width = 6, height = 4, units = "in")
ggsave("SxC_CPVT.png", plot = plot_out_histogram(SxC), width = 6, height = 4, units = "in")
ggsave("Lh_CPVT.png", plot = plot_out_histogram(Lh), width = 6, height = 4, units = "in")
ggsave("GgC_CPVT.png", plot = plot_out_histogram(GgC), width = 6, height = 4, units = "in")
ggsave("SuC_CPVT.png", plot = plot_out_histogram(SuC), width = 6, height = 4, units = "in")
ggsave("BlB_CPVT.png", plot = plot_out_histogram(BlB), width = 6, height = 4, units = "in")
ggsave("VeD_CPVT.png", plot = plot_out_histogram(VeD), width = 6, height = 4, units = "in")
#W is water; no butternuts
ggsave("HnA_CPVT.png", plot = plot_out_histogram(HnA), width = 6, height = 4, units = "in")
ggsave("HnB_CPVT.png", plot = plot_out_histogram(HnB), width = 6, height = 4, units = "in")
ggsave("FsB_CPVT.png", plot = plot_out_histogram(FsB), width = 6, height = 4, units = "in")
ggsave("Le_CPVT.png", plot = plot_out_histogram(Le), width = 6, height = 4, units = "in")
ggsave("MnC_CPVT.png", plot = plot_out_histogram(MnC), width = 6, height = 4, units = "in")
ggsave("FaE_CPVT.png", plot = plot_out_histogram(FaE), width = 6, height = 4, units = "in")
ggsave("FsE_CPVT.png", plot = plot_out_histogram(FsE), width = 6, height = 4, units = "in")
write.csv(CPVT_real_vs_random_quantiles, "CPVT_butternut_proportion_quantiles_by_soil_type.csv")

#------------
#Plots
#overall soil maps
tm_shape(ILM_polygon) +
  tm_polygons() +
  tm_shape(ILM_data_w_soils) +
  tm_dots("MUSYM")

tm_shape(CPVT_polygon) +
  tm_polygons() +
  tm_shape(CPVT_data_w_soils) +
  tm_dots("MUSYM")

#HEALTH PLOTS
#these next few split up our data by age, and also select all of the living trees (which have more data to their name)
data_ILMadults_w_soil <- ILM_data_w_soils %>% filter(seedling_or_adult == "Adult")
data_ILMliveadults_w_soil <- data_ILMadults_w_soil %>% filter(adult_dead_or_alive == "Alive")

data_ILMseedlings_w_soil <- ILM_data_w_soils %>% filter(seedling_or_adult == "Seedling")
data_ILMliveseedlings_w_soil <- data_ILMseedlings_w_soil %>% filter(seedling_dead_or_alive == "Alive")

data_CPVTadults_w_soil <- CPVT_data_w_soils %>% filter(seedling_or_adult == "Adult")
data_CPVTliveadults_w_soil <- data_CPVTadults_w_soil %>% filter(adult_dead_or_alive == "Alive")

data_CPVTseedlings_w_soil <- CPVT_data_w_soils %>% filter(seedling_or_adult == "Seedling")
data_CPVTliveseedlings_w_soil <- data_CPVTseedlings_w_soil %>% filter(seedling_dead_or_alive == "Alive")

#let's just see how many butternuts we have at each soil type
barplot(ILM_data_w_soils, MUSYM, "Soil at ILM")
barplot(data_ILMseedlings_w_soil, MUSYM, "Soil for ILM seedlings")

barplot(CPVT_data_w_soils, MUSYM, "Soil at CPVT")
barplot(data_CPVTseedlings_w_soil, MUSYM, "Soil for CPVT seedlings")

#these plots will show any trends in health status across the two sites based on soil type
boxplot(data_ILMliveadults_w_soil, MUSYM, adult_percent_live_canopy, "Canopy at ILM by soil") + stat_n_text()
boxplot(data_CPVTliveadults_w_soil, MUSYM, adult_percent_live_canopy, "Canopy at CPVT by soil") + stat_n_text()

boxplot(data_ILMliveadults_w_soil, MUSYM, live_adult_girdle, "Girdle at ILM by soil") + stat_n_text()
boxplot(data_CPVTliveadults_w_soil, MUSYM, live_adult_girdle, "Girdle at CPVT by soil") + stat_n_text()

boxplot(data_ILMliveadults_w_soil, MUSYM, live_adult_purdue_canker_rating, "Purdue canker ranking at ILM by soil") + stat_n_text()
boxplot(data_CPVTliveadults_w_soil, MUSYM, live_adult_purdue_canker_rating, "Purdue canker ranking at CPVT by soil") + stat_n_text()

boxplot(data_ILMliveadults_w_soil, MUSYM, live_adult_purdue_canopy_ranking, "Purdue canopy ranking at ILM by soil") + stat_n_text()
boxplot(data_CPVTliveadults_w_soil, MUSYM, live_adult_purdue_canopy_ranking, "Purdue canopy ranking at CPVT by soil") + stat_n_text()