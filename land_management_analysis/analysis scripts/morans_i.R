#Catherine dell'Olio
#Moran's I testing of spatial autocorrelation in trees' canopy and canker metrics
#based on Rebecca Wanger's script https://github.com/HobanLab/QUBR_GenGeoEcoDemoCorr/blob/b85cbd70de1591ef81ef08558da2d392a1040a42/analyses/hypothesis_2_moran_CLEANED.R

# The purpose of this script is to evaluated whether the canker and canopy metrics (health status) of butternut trees is impacted by the distance to other individuals of the same species, either due to competition or facilitation.
# To test this, we used Global and Local Moran's I to determine whether values of adult_percent_live_canopy and live_adult_girdle that were closer together were more similar in value or not. 
# The global Moran's I looked for general spatial autocorrelation
# The local Moran's I looked for areas were values were more similar than other areas. 

# The script is broken into sections of 
# 1) loading and processing the packages and spatial/health data for the trees/sites, 
# 2) creating a dataframe with the coordinates and average distance between each tree and their 5 nearest neighbors
# 3) graphing descriptive summary histograms and calculating summary statistics (mean, median, sd, etc.) for each response variable, and
# 4) Running the global and local Moran's I analyses, 

# NOTE: Uncomment and run line 30, sourcing DataProcessing.R, if the line has not yet to be run across any of the scripts/the environment has been cleared 

#### Loading libraries and relevant data ####

library(tidyverse)
library(broom) #for function tidy() to save outputs of statistical tests
library(moments) # for casdulating the moments of each variable
library(sf) # for plotting spatial objects
library(spatstat) # to run the nndist function
library(spdep) # to use Moran's I functions like lag.listw
library(ape) # for computing the Moran's I stat
library(raster) #to use point distance
library(nlme) # linear mixed effect models
library(MuMIn) #to be able to use model.sel for fitting linear models with spatial autocorrelation
library(geoR) # to be able to use variograms with the lme, requires XQuartz from 
library(Kendall)# to use the Kendall's Tau test to look for non-parametric correlations in the data

#set your working directory
#working_directory <- "~/GitHub/butternut-health-assessment-2025/land_management_analysis/"
setwd(working_directory)
setwd("saved plots/morans_i")

# loading in the processed tree data 
# NOTE: Uncomment and run line 41, sourcing Data_Processing_Script.R, if the line has not yet to be run across any of the scripts/the environment has been cleared 
#source("./analysis scripts/DataProcessing.R")

# Make a function that is the opposite of the %in% function
`%notin%` <- Negate(`%in%`) 

#### Moran's I ####

# For each response variable, we ran these analyses:
# 1) Global Moran's I
# a) We first created a matrix of the inverse distances between each tree 
# b) computed the average response variable of the neighboring trees for each tree (lagged response variable)
# c) plot the regression line of the lagged response variable vs. the actual response values for each tree
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
# d) Calculate the Moran's I statistic and then use a Monte Carlo test and plot to see if it is significant spatial autocorrelation 
# 2) Local Moran's I
# a) calculate and plot the expected local Moran's I vs. local Moran's I for each tree 
# b) calculate the number of trees with significant local Moran's I (p-values < 0.05)


###Test for ILM###

#creating a shapefile for plotting later
plot_box <- ILM_box
dataframe <- ILM_fixed_field_data_processed

## Global Moran's I

#removing NAs in our health metric (needed to compute nearest neighbor and the metric of focus)
dataframe <- dataframe %>%
  filter(!is.na(adult_percent_live_canopy)) %>%
  filter(!is.na(.data$adult_percent_live_canopy))

#creating a matrix of the tree locations
tree.coord.matrix <- as.matrix(cbind(dataframe$X, 
                                     dataframe$Y))

#creates nearest neighbor matrix of the tree coordinates within 40 meters of the mean of our health metric of the population
knn.dist <- dnearneigh(tree.coord.matrix, d1 = 0, d2 = (40*mean(dataframe$adult_percent_live_canopy)))

#inverse distance weighting with raw distance-based weights without applying any normalization
lw.dist <- nb2listwdist(knn.dist, dataframe, type="idw", style="W", #fixed_field_data_processed_sf_trans_coordinates
                        alpha = 1, dmax = NULL, longlat = NULL, zero.policy=T) # had to set zero.policy to true because of empty neighbor sets

#creating lags for each tree, which computes the average neighboring health metric for each tree
dataframe$lag.metric <- lag.listw(lw.dist, dataframe$adult_percent_live_canopy)

# Create a regression model of the lagged response variable (average amongst closest trees) vs. the known response variable 
lm <- lm(lag.metric ~ dataframe$adult_percent_live_canopy, dataframe)

#computing the Moran's I statistic
global.moran.I <- moran(dataframe$adult_percent_live_canopy, listw = lw.dist, n = length(lw.dist$neighbours), S0 = Szero(lw.dist))
global.moran.I

#assessing statistical significance with a Monte-Carlo simulation
MC.LM.metric <- moran.mc(dataframe$adult_percent_live_canopy, lw.dist, nsim = 999)
MC.LM.metric

#plot of simulated Moran's I values against our value
plot(MC.LM.metric, main=paste0("Global Moran's I MC values ILM PLC, p value = ", MC.LM.metric$p.value), las=1, xlab = "Percent live canopy")

cat(paste0("Global Moran's I: ", global.moran.I$I, "    Global Moran's I Monte Carlo P-value: ", MC.LM.metric$p.value), file = "Global Morans I and MC P value ILM PLC.txt")

## Local Morcapture.output()## Local Moran's I

#setting a seed for all of the results from using the morans_I() function because localmoran_perm() uses random permutations
set.seed(25)

#using the weighted neighbors to simulate health metric values at random
MC_local <- localmoran_perm(dataframe$adult_percent_live_canopy, lw.dist, nsim = 9999, alternative = "greater")
MC_local.df <- as.data.frame(MC_local)

#calculating the p-values for each individual tree Moran's I, observed vs. expected
dataframe$p.metric  <- MC_local.df$`Pr(folded) Sim`
#adjusting the p-vlaues to take into account multiple tests
dataframe$p.metric.adjusted <- p.adjust(dataframe$p.metric, 
                                        method = "fdr", n=length(dataframe$p.metric))

#Number of trees with significant adjusted Moran's I
sig.tree <- length(which(dataframe$p.metric.adjusted < 0.05))

#filtering out significant p-values
dataframe_sign <- dataframe %>%
  mutate(pval_sig = p.metric.adjusted <= .05) %>%
  filter(pval_sig == T)

results <- list(dataframe$lag.metric, lm, MC.LM.metric,
            MC_local.df, sig.tree, dataframe_sign, dataframe$p.metric.adjusted, global.moran.I, dataframe)


#regression for ANN canopy metric vs. tree metric 
results[[2]]

#Monte Carlo Simulation for Global Moran's I
results[[3]]

#storing the dataframe used for the Moran's Is with the NAs removed
ILM_results_dataframe <- results[[9]]


#creating a column for the lagged canopy metric
ILM_results_dataframe$lag_adult_percent_live_canopy <- results[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, healthier trees are closer together and less healthy trees are closer together
# negative slope, negative spatial autocorrelation, variation in health(canopy) of trees close together
P <- ggplot(data=ILM_results_dataframe, aes(x=adult_percent_live_canopy, y=lag_adult_percent_live_canopy))+
  geom_point()+
  geom_smooth(method = "lm", col="blue")+
  xlab("Percent live canopy")+
  ylab("Lagged percent live canopy")
P
ggsave("Values vs Neighboring Values ILM PLC.png", plot = P, width = 6, height = 4, units = "in")

#Local Moran's I 

#assigning a Monta Carlo dataframe for plotting
MC_local_ILM_results_df <- results[[4]]

#number of trees with Significant Local Moran's I
results[[5]]

#assigning the trees with the significant local Moran's I to a dataframe
ILM_fixed_field_data_processed_sign <- results[[6]]

#assigning the p-values of the adjusted local Moran's I to a dataframe
ILM_results_dataframe$p_live_canopy <- results[[7]]

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local_ILM_results_df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic")+
  ylab("Expected Moran's I Statistic")+
  theme_gray()

#plotting the local Moran's I 
ggplot() +
  geom_sf(data =trails_ILM_trans) +
  geom_sf(data =ILM_results_dataframe, aes(color = p_live_canopy)) +
  geom_sf(data = ILM_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  coord_sf(xlim = c(ILM_box[1], ILM_box[3]), ylim = c(ILM_box[2], ILM_box[4]))+
  labs(color = "Adjusted P Value for PLC")


#attempting to zoom on the health metrics of the significant points
P <- ggplot() +
  geom_sf(data =trails_ILM_trans) +
  geom_sf(data =ILM_results_dataframe, aes(size = adult_percent_live_canopy)) +
  geom_sf(data =ILM_results_dataframe, aes(color = p_live_canopy))+
  geom_sf(data = ILM_fixed_field_data_processed_sign, color = "red", aes(fill = "red"))
P
ggsave("Local Morans I for PLC ILM.png", plot = P, width = 6,height = 4, units = "in")
#it's the less healthy ones that seem to be clustered??

###Live adult girdle

#creating a shapefile for plotting later
plot_box <- ILM_box
dataframe <- ILM_fixed_field_data_processed

## Global Moran's I

#removing NAs in our health metric (needed to compute nearest neighbor and the metric of focus)
dataframe <- dataframe %>%
  filter(!is.na(live_adult_girdle)) %>%
  filter(!is.na(.data$live_adult_girdle))

#creating a matrix of the tree locations
tree.coord.matrix <- as.matrix(cbind(dataframe$X, 
                                     dataframe$Y))

#creates nearest neighbor matrix of the tree coordinates within 40 meters of the mean of our health metric of the population
knn.dist <- dnearneigh(tree.coord.matrix, d1 = 0, d2 = (40*mean(dataframe$live_adult_girdle)))

#inverse distance weighting with raw distance-based weights without applying any normalization
lw.dist <- nb2listwdist(knn.dist, dataframe, type="idw", style="W", #fixed_field_data_processed_sf_trans_coordinates
                        alpha = 1, dmax = NULL, longlat = NULL, zero.policy=T) # had to set zero.policy to true because of empty neighbor sets

#creating lags for each tree, which computes the average neighboring health metric for each tree
dataframe$lag.metric <- lag.listw(lw.dist, dataframe$live_adult_girdle)

# Create a regression model of the lagged response variable (average amongst closest trees) vs. the known response variable 
lm <- lm(lag.metric ~ dataframe$live_adult_girdle, dataframe)

#computing the Moran's I statistic
global.moran.I <- moran(dataframe$live_adult_girdle, listw = lw.dist, n = length(lw.dist$neighbours), S0 = Szero(lw.dist))
global.moran.I

#assessing statistical significance with a Monte-Carlo simulation
MC.LM.metric <- moran.mc(dataframe$live_adult_girdle, lw.dist, nsim = 999)
MC.LM.metric

#plot of simulated Moran's I values against our value
plot(MC.LM.metric, main=paste0("Global Morans I MC Values PG ILM, P value = ", MC.LM.metric$p.value), las=1, xlab = "Percent girdled")
#saved manually

cat(paste0("Global Moran's I: ", global.moran.I$I, "    Global Moran's I Monte Carlo P-value: ", MC.LM.metric$p.value), file = "Global Morans I and MC P value ILM PG.txt")

## Local Moran's I

#setting a seed for all of the results from using the morans_I() function because localmoran_perm() uses random permutations
set.seed(25)

#using the weighted neighbors to simulate health values at random
MC_local <- localmoran_perm(dataframe$live_adult_girdle, lw.dist, nsim = 9999, alternative = "greater")
MC_local.df <- as.data.frame(MC_local)

#calculating the p-values for each individual tree Moran's I, observed vs. expected
dataframe$p.metric  <- MC_local.df$`Pr(folded) Sim`
#adjusting the p-vlaues to take into account multiple tests
dataframe$p.metric.adjusted <- p.adjust(dataframe$p.metric, 
                                        method = "fdr", n=length(dataframe$p.metric))

#Number of trees with significant adjusted Moran's I
sig.tree <- length(which(dataframe$p.metric.adjusted < 0.05))

#filtering out significant p-values
dataframe_sign <- dataframe %>%
  mutate(pval_sig = p.metric.adjusted <= .05) %>%
  filter(pval_sig == T)

results <- list(dataframe$lag.metric, lm, MC.LM.metric,
                MC_local.df, sig.tree, dataframe_sign, dataframe$p.metric.adjusted, global.moran.I, dataframe)


#regression for ANN canopy metric vs. tree metric 
results[[2]]

#Monte Carlo Simulation for Global Moran's I
results[[3]]

#storing the dataframe used for the Moran's Is with the NAs removed
ILM_results_dataframe <- results[[9]]


#creating a column for the lagged canker metric
ILM_results_dataframe$lag_live_adult_girdle <- results[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, healthier trees are closer together and less healthy trees are closer together
# negative slope, negative spatial autocorrelation, variation in health(canopy) of trees close together
P <- ggplot(data=ILM_results_dataframe, aes(x=live_adult_girdle, y=lag_live_adult_girdle))+
  geom_point()+
  geom_smooth(method = "lm", col="blue")+
  xlab("Percent girdled")+
  ylab("Lagged percent girdled")
P
ggsave("Values vs Neighboring Values ILM PG.png", plot = P, width = 6, height = 4, units = "in")

#Local Moran's I 

#assigning a Monta Carlo dataframe for plotting
MC_local_ILM_results_df <- results[[4]]

#number of trees with Significant Local Moran's I
results[[5]]

#assigning the trees with the significant local Moran's I to a dataframe
ILM_fixed_field_data_processed_sign <- results[[6]]

#assigning the p-values of the adjusted local Moran's I to a dataframe
ILM_results_dataframe$p_live_girdle <- results[[7]]

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local_ILM_results_df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic")+
  ylab("Expected Moran's I Statistic")+
  theme_gray()

#plotting the local Moran's I 
ggplot() +
  geom_sf(data =trails_ILM_trans) +
  geom_sf(data =ILM_results_dataframe, aes(color = p_live_girdle)) +
  geom_sf(data = ILM_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  coord_sf(xlim = c(ILM_box[1], ILM_box[3]), ylim = c(ILM_box[2], ILM_box[4]))+
  labs(color = "Adjusted P Value for PLC")

#attempting to zoom on the healths of the significant point
P <- ggplot() +
  geom_sf(data =trails_ILM_trans) +
  geom_sf(data =ILM_results_dataframe, aes(size = live_adult_girdle)) +
  geom_sf(data =ILM_results_dataframe, aes(color = p_live_girdle))+
  geom_sf(data = ILM_fixed_field_data_processed_sign, color = "red", aes(fill = "red"))
P
ggsave("Local Morans I PG ILM.png", plot = P, width = 6, height = 4, units = "in")


###Test for CPVT###


#creating a shapefile for plotting later
plot_box <- CPVT_box
dataframe <- CPVT_fixed_field_data_processed

## Global Moran's I

#removing NAs in our health metric (needed to compute nearest neighbor and the metric of focus)
dataframe <- dataframe %>%
  filter(!is.na(adult_percent_live_canopy)) %>%
  filter(!is.na(.data$adult_percent_live_canopy))

#creating a matrix of the tree locations
tree.coord.matrix <- as.matrix(cbind(dataframe$X, 
                                     dataframe$Y))

#creates nearest neighbor matrix of the tree coordinates within 40 meters of the mean of our health metric of the population
knn.dist <- dnearneigh(tree.coord.matrix, d1 = 0, d2 = (40*mean(dataframe$adult_percent_live_canopy)))

#inverse distance weighting with raw distance-based weights without applying any normalization
lw.dist <- nb2listwdist(knn.dist, dataframe, type="idw", style="W", #fixed_field_data_processed_sf_trans_coordinates
                        alpha = 1, dmax = NULL, longlat = NULL, zero.policy=T) # had to set zero.policy to true because of empty neighbor sets

#creating lags for each tree, which computes the average neighboring size metric for each tree
dataframe$lag.metric <- lag.listw(lw.dist, dataframe$adult_percent_live_canopy)

# Create a regression model of the lagged response variable (average amongst closest trees) vs. the known response variable 
lm <- lm(lag.metric ~ dataframe$adult_percent_live_canopy, dataframe)

#computing the Moran's I statistic
global.moran.I <- moran(dataframe$adult_percent_live_canopy, listw = lw.dist, n = length(lw.dist$neighbours), S0 = Szero(lw.dist))
global.moran.I

#assessing statistical significance with a Monte-Carlo simulation
MC.LM.metric <- moran.mc(dataframe$adult_percent_live_canopy, lw.dist, nsim = 999)
MC.LM.metric

#plot of simulated Moran's I values against our value
plot(MC.LM.metric, main=paste0("Global Morans I MC Values PLC CPVT, p value = ", MC.LM.metric$p.value), las=1, xlab = "Percent live canopy")
#save manually

cat(paste0("Global Moran's I: ", global.moran.I$I, "    Global Moran's I Monte Carlo P-value: ", MC.LM.metric$p.value), "Global Morans I and MC P value CPVT PLC.txt")

## Local Moran's I

#setting a seed for all of the results from using the morans_I() function because localmoran_perm() uses random permutations
set.seed(25)

#using the weighted neighbors to simulate size values at random
MC_local <- localmoran_perm(dataframe$adult_percent_live_canopy, lw.dist, nsim = 9999, alternative = "greater")
MC_local.df <- as.data.frame(MC_local)

#calculating the p-values for each individual tree Moran's I, observed vs. expected
dataframe$p.metric  <- MC_local.df$`Pr(folded) Sim`
#adjusting the p-vlaues to take into account multiple tests
dataframe$p.metric.adjusted <- p.adjust(dataframe$p.metric, 
                                        method = "fdr", n=length(dataframe$p.metric))

#Number of trees with significant adjusted Moran's I
sig.tree <- length(which(dataframe$p.metric.adjusted < 0.05))

#filtering out significant p-values
dataframe_sign <- dataframe %>%
  mutate(pval_sig = p.metric.adjusted <= .05) %>%
  filter(pval_sig == T)

results <- list(dataframe$lag.metric, lm, MC.LM.metric,
                MC_local.df, sig.tree, dataframe_sign, dataframe$p.metric.adjusted, global.moran.I, dataframe)


#regression for ANN canopy metric vs. tree metric 
results[[2]]

#Monte Carlo Simulation for Global Moran's I
results[[3]]

#storing the dataframe used for the Moran's Is with the NAs removed
CPVT_results_dataframe <- results[[9]]


#creating a column for the lagged canopy metric
CPVT_results_dataframe$lag_adult_percent_live_canopy <- results[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, healthier trees are closer together and less healthy trees are closer together
# negative slope, negative spatial autocorrelation, variation in health(canopy) of trees close together
P <- ggplot(data=CPVT_results_dataframe, aes(x=adult_percent_live_canopy, y=lag_adult_percent_live_canopy))+
  geom_point()+
  geom_smooth(method = "lm", col="blue")+
  xlab("Percent live canopy")+
  ylab("Lagged percent live canopy")
P
ggsave("CPVT PLC values vs. neighboring (lagged) values.png", plot = P, width = 6, height = 4, units = "in")

#Local Moran's I 

#assigning a Monta Carlo dataframe for plotting
MC_local_CPVT_results_df <- results[[4]]

#number of trees with Significant Local Moran's I
results[[5]]

#assigning the trees with the significant local Moran's I to a dataframe
CPVT_fixed_field_data_processed_sign <- results[[6]]

#assigning the p-values of the adjusted local Moran's I to a dataframe
CPVT_results_dataframe$p_live_canopy <- results[[7]]

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local_CPVT_results_df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic")+
  ylab("Expected Moran's I Statistic")+
  theme_gray()

#plotting the local Moran's I 
ggplot() +
  geom_sf(data =trails_CPVT_trans) +
  geom_sf(data =CPVT_results_dataframe, aes(color = p_live_canopy)) +
  geom_sf(data = CPVT_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  coord_sf(xlim = c(CPVT_box[1], CPVT_box[3]), ylim = c(CPVT_box[2], CPVT_box[4]))+
  labs(color = "Adjusted P Value for PLC")

#attempting to zoom on the sizes of the significant point
P <- ggplot() +
  geom_sf(data =trails_CPVT_trans) +
  geom_sf(data =CPVT_results_dataframe, aes(size = adult_percent_live_canopy)) +
  geom_sf(data =CPVT_results_dataframe, aes(color = p_live_canopy))+
  geom_sf(data = CPVT_fixed_field_data_processed_sign, color = "red", aes(fill = "red"))
P
ggsave("Local Morans I Values CPVT PLC.png", plot = P, width = 6, height = 4, units = "in")

###Live adult girdle

#creating a shapefile for plotting later
plot_box <- CPVT_box
dataframe <- CPVT_fixed_field_data_processed

## Global Moran's I

#removing NAs in our health metric (needed to compute nearest neighbor and the metric of focus)
dataframe <- dataframe %>%
  filter(!is.na(live_adult_girdle)) %>%
  filter(!is.na(.data$live_adult_girdle))

#creating a matrix of the tree locations
tree.coord.matrix <- as.matrix(cbind(dataframe$X, 
                                     dataframe$Y))

#creates nearest neighbor matrix of the tree coordinates within 40 meters of the mean of our health metric of the population
knn.dist <- dnearneigh(tree.coord.matrix, d1 = 0, d2 = (40*mean(dataframe$live_adult_girdle)))

#inverse distance weighting with raw distance-based weights without applying any normalization
lw.dist <- nb2listwdist(knn.dist, dataframe, type="idw", style="W", #fixed_field_data_processed_sf_trans_coordinates
                        alpha = 1, dmax = NULL, longlat = NULL, zero.policy=T) # had to set zero.policy to true because of empty neighbor sets

#creating lags for each tree, which computes the average neighboring health metric for each tree
dataframe$lag.metric <- lag.listw(lw.dist, dataframe$live_adult_girdle)

# Create a regression model of the lagged response variable (average amongst closest trees) vs. the known response variable 
lm <- lm(lag.metric ~ dataframe$live_adult_girdle, dataframe)

#computing the Moran's I statistic
global.moran.I <- moran(dataframe$live_adult_girdle, listw = lw.dist, n = length(lw.dist$neighbours), S0 = Szero(lw.dist))
global.moran.I

#assessing statistical significance with a Monte-Carlo simulation
MC.LM.metric <- moran.mc(dataframe$live_adult_girdle, lw.dist, nsim = 999)
MC.LM.metric

#plot of simulated Moran's I values against our value
plot(MC.LM.metric, main=paste0("Global Morans I Values CPVT PG, p value = ", MC.LM.metric$p.value), las=1, xlab = "Percent girdled")
#save manually


cat(paste0("Global Moran's I: ", global.moran.I$I, "     Global Moran's I Monte Carlo P-value: ", MC.LM.metric$p.value), "Global Morans I and MC P Value CPVT PG.txt")

## Local Moran's I

#setting a seed for all of the results from using the morans_I() function because localmoran_perm() uses random permutations
set.seed(25)

#using the weighted neighbors to simulate health values at random
MC_local <- localmoran_perm(dataframe$live_adult_girdle, lw.dist, nsim = 9999, alternative = "greater")
MC_local.df <- as.data.frame(MC_local)

#calculating the p-values for each individual tree Moran's I, observed vs. expected
dataframe$p.metric  <- MC_local.df$`Pr(folded) Sim`
#adjusting the p-vlaues to take into account multiple tests
dataframe$p.metric.adjusted <- p.adjust(dataframe$p.metric, 
                                        method = "fdr", n=length(dataframe$p.metric))

#Number of trees with significant adjusted Moran's I
sig.tree <- length(which(dataframe$p.metric.adjusted < 0.05))

#filtering out significant p-values
dataframe_sign <- dataframe %>%
  mutate(pval_sig = p.metric.adjusted <= .05) %>%
  filter(pval_sig == T)

results <- list(dataframe$lag.metric, lm, MC.LM.metric,
                MC_local.df, sig.tree, dataframe_sign, dataframe$p.metric.adjusted, global.moran.I, dataframe)


#regression for ANN canopy metric vs. tree metric 
results[[2]]

#Monte Carlo Simulation for Global Moran's I
results[[3]]

#storing the dataframe used for the Moran's Is with the NAs removed
CPVT_results_dataframe <- results[[9]]


#creating a column for the lagged canker metric
CPVT_results_dataframe$lag_live_adult_girdle <- results[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, healthier trees are closer together and less healthy trees are closer together
# negative slope, negative spatial autocorrelation, variation in health(canopy) of trees close together
P <- ggplot(data=CPVT_results_dataframe, aes(x=live_adult_girdle, y=lag_live_adult_girdle))+
  geom_point()+
  geom_smooth(method = "lm", col="blue")+
  xlab("Percent girdled")+
  ylab("Lagged percent girdled")
P
ggsave("CPVT PG values vs. neighboring (lagged) values.png", plot = P, width = 6, height = 4, units = "in")
#Local Moran's I 

#assigning a Monta Carlo dataframe for plotting
MC_local_CPVT_results_df <- results[[4]]

#number of trees with Significant Local Moran's I
results[[5]]

#assigning the trees with the significant local Moran's I to a dataframe
CPVT_fixed_field_data_processed_sign <- results[[6]]

#assigning the p-values of the adjusted local Moran's I to a dataframe
CPVT_results_dataframe$p_live_girdle <- results[[7]]

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local_CPVT_results_df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic")+
  ylab("Expected Moran's I Statistic")+
  theme_gray()

#plotting the local Moran's I 
ggplot() +
  geom_sf(data =trails_CPVT_trans) +
  geom_sf(data =CPVT_results_dataframe, aes(color = p_live_girdle)) +
  geom_sf(data = CPVT_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  coord_sf(xlim = c(CPVT_box[1], CPVT_box[3]), ylim = c(CPVT_box[2], CPVT_box[4]))+
  labs(color = "Adjusted P Value for PLC")

#attempting to zoom on the healths of the significant point
P <- ggplot() +
  geom_sf(data =trails_CPVT_trans) +
  geom_sf(data =CPVT_results_dataframe, aes(size = live_adult_girdle)) +
  geom_sf(data =CPVT_results_dataframe, aes(color = p_live_girdle))+
  geom_sf(data = CPVT_fixed_field_data_processed_sign, color = "red", aes(fill = "red"))
P
ggsave("Local Morans I CPVT PG.png", plot = P, width = 6, height= 4, units = "in")