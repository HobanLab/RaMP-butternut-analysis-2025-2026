#RUN_ME
#Catherine dell'Olio
#dataAnalysis, associate trees experiment

#This project will help me analyze my associate trees protocol data, collected in March/April fieldwork trip.  In this trip, I am collecting information about the associate trees around a subsample of butternuts at a site in northwestern Vermont.  These trees were selected to represent the wide range of health statuses of butternuts at this field site.  I think of these as belonging to four distinct groups: those that are healthy with little evidence of canker, those that are very healthy despite a lot of canker, those dying but little canker, and those dying with a lot of canker.  I am recording all of the trees that occur within a 10 meter radius around each butternut selected---their height, species, distance to the focal tree, DBH, and canopy width.

#I can use the DBH and distance to the focal tree for all individuals of a given species to calculate an impact coefficient of each species around the butternut, and determine the impacts (both scale and direction) for species around the butternut.

#The basic steps of this analysis process are:
#0 (because it is done BEFORE any script is run): Clean up data in google sheets
#1: load data (and all data subsets), cleaning the data
#2: looking for species clustering
#3: initial data familiarization
#4: determining a model (and testing the model-fitting process)
#5: testing the chosen model

#Starting off strong with step 1: the initial data cleaning in R.  There are more extensive notes in this script that I recommend following along with, but essentially, this script makes the datasheets formed in Google Sheets easier to work with in R (setting column names properly, setting variables as numerics, converting percentages, etc.)
setwd("~/assoc_trees_analysis_cfd/scripts")
source("STEP_1_data_upload_and_cleaning.R")

#Step 2: I looked for initial species clustering with this script to get a sense of what types of communities were seen.  I include both dead and alive trees in these communities to get a sense of what they look like over a longer term, so this is only run with abundance data including both Dead and Alive (DA) trees. 
#this will only work with abundance_DA!
df <- abundance_DA
setwd("~/assoc_trees_analysis_cfd/scripts")
source("STEP_2_species_clustering.R")
#In this, I found that the clustering of trees doesn't have a super strong signal; the most important species driving community differences is also the species just by far and away the most common (red cedar).
#There is borderline too much stress to really trust the 2d visualization of community differences, but it seems to be somewhat differentiated by the species clustering given in the dendrogram.  Again, most of this is driven by differences in abundance of red cedar.

#This is when we create three brazillion plots to demonstrate the way the data looks
df <- impact_TOTAL
#Brace yourself; this will produce, conservatively, 9.6 billion plots
source("STEP_3_initial_data_familiarization.R")
#It's worth going through all of this pretty carefully to really understand the distribution and behavior of the data.
#Note that the correlations between some species categories are significant, but none of those are between two species with enough data that we will be testing for both, and thus have to worry about correlations between predictor variables (a big no-no in GLMs)

#Now it is time to determine a model
source("STEP_4_determining_a_model.R")
#IMPORTANT NOTE: this creates a version of our data, where the response variables fit the (0,1) range specified by the beta distribution so we can use beta regression, by transforming 0 and 1 values to 0.1 and 0.99.  So note that df is now df_without_boundaries (those boundary values being 0 and 1).  I have this specified below as "response", which we will need in step 5 anyways
response <- df_wo_boundaries$adult_percent_live_canopy
modelSelection(df_wo_boundaries$adult_percent_live_canopy)
#for our models based on impact_TOTAL (total impact of all trees, alive or dead, we find that alternate-leaf dogwood is the best explainer of percent live canopy (K=3, AICc = -26.49).

source("STEP_5_testing_the_model.R")
chosen_model <- glmmTMB(adult_percent_live_canopy ~ `alternate-leaf_dogwood`, data = df_wo_boundaries, family = beta_family(link = "logit"))
modelTest(chosen_model)
#model coefficients
# [[1]]
# Family: beta  ( logit )
# Formula:          adult_percent_live_canopy ~ `alternate-leaf_dogwood`
# Data: df_wo_boundaries
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# -27.0     -20.9      16.5     -33.0        52 
# 
# 
# Dispersion parameter for beta family (): 1.53 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)               0.70028    0.17970   3.897 9.74e-05 ***
#   `alternate-leaf_dogwood` -0.22259    0.07727  -2.881  0.00397 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# residuals of model (should be close to 0:
# [1] 0.03073619
modelOutputPlots(chosen_model)
#these look pretty good!

#now it is time to test for spatial autocorrelation.  This script also relies on chosen_model, df, and response (all predefined!)
source("STEP_6_morans_i.R")
morans_i_test(chosen_model, df)
# Moran I test under randomisation
# 
# data:  residuals_regression  
# weights: lw    
# 
# Moran I statistic standard deviate = 2.1731, p-value = 0.01489
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# 0.14102901       -0.01851852        0.00539048 
moran_plots(chosen_model, response)
#these all look good as well.  Moran's I globally shows some positive autocorrelation in percent live canopy (so trees with similar values cluster together) and don't portray much of a local pattern

#now I will shift into looking a different response live_adult_girdle (back to step 4)
#just a reminder, we are still looking at impact_TOTAL
response <- df_wo_boundaries$live_adult_girdle
modelSelection(df_wo_boundaries$live_adult_girdle)
#for our models based on impact_TOTAL (total impact of all trees, alive or dead, we find that bitternut_hickory + american elm + alternate leaf dogwood is the best explainer of percent girdle (K=5, AICc = -6.42), but that butternut is a large player as well.
chosen_model <- glmmTMB(adult_percent_live_canopy ~ bitternut_hickory+american_elm+`alternate-leaf_dogwood`, data = df_wo_boundaries, family = beta_family(link = "logit"))
source("STEP_5_testing_the_model.R")
modelTest(chosen_model)
# [[1]]
# Family: beta  ( logit )
# Formula:          adult_percent_live_canopy ~ bitternut_hickory + american_elm +  
#   `alternate-leaf_dogwood`
# Data: df_wo_boundaries
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# -23.4     -13.3      16.7     -33.4        50 
# 
# 
# Dispersion parameter for beta family (): 1.54 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)               0.78753    0.23241   3.389 0.000703 ***
#   bitternut_hickory        -0.01684    0.03522  -0.478 0.632493    
# american_elm             -0.01169    0.02627  -0.445 0.656260    
# `alternate-leaf_dogwood` -0.21100    0.08248  -2.558 0.010524 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# [[2]]
# [1] 0.0302792
modelOutputPlots(chosen_model)
#looks even better!

source("STEP_6_morans_i.R")
morans_i_test(chosen_model, df)
#spatial autocorrelation found!
# Moran I test under randomisation
# 
# data:  residuals_regression  
# weights: lw    
# 
# Moran I statistic standard deviate = 1.9225, p-value = 0.02727
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# 0.122498407      -0.018518519       0.005380433 

#It's worth noting that we measured many different metrics of canker (percent girdle, canker to the first live branch, canker at the base, Purdue canker rating).  I use percent girdle here because I think it gives a good sense of the amount of stress canker is creating for the tree's nutrient flow up the trunk.

#Now I'd also like to see if the models chosen differ meaningfully under three different conditions; if we use a more complex version of impact coefficient that takes into account the height and directionality of the associate trees (impact_height_direction_TOTAL), a relativized version of impact coefficient (impact_PERCENT), and finally, if we calculate impact coefficients only from still-living trees.  The response variables tested (percent live canopy and percent girdling) remain the same.

#--------------------------------------------------------
#ACCOUNTING FOR THE SIZE AND DIRECTION TO ASSOCIATE TREES
df <- impact_height_direction_TOTAL
source("STEP_3_initial_data_familiarization.R")
#similar spread of data and correlations found, which isn't surprising

#Now it is time to determine a model
response <- df_wo_boundaries$adult_percent_live_canopy
source("STEP_4_determining_a_model.R")
modelSelection(df_wo_boundaries$adult_percent_live_canopy)
#Model chosen is the same as the one chosen when using plain impact_TOTAL, based only on `alternate-leaf_dogwood` impact (K = 3 AICc = -22.65)

source("STEP_5_testing_the_model.R")
chosen_model <- glmmTMB(adult_percent_live_canopy ~ `alternate-leaf_dogwood`, data = df_wo_boundaries, family = beta_family(link = "logit"))
modelTest(chosen_model)
# [[1]]
# Family: beta  ( logit )
# Formula:          adult_percent_live_canopy ~ `alternate-leaf_dogwood`
# Data: df_wo_boundaries
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# -23.1     -17.1      14.6     -29.1        52 
# 
# 
# Dispersion parameter for beta family (): 1.42 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)                0.5798     0.1784   3.250  0.00116 **
#   `alternate-leaf_dogwood`  -0.3609     0.1396  -2.585  0.00973 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# [[2]]
# [1] 0.03555328
#the intercept, coefficient, and dispersion parameter are only slightly different, but within the standard deviation of its estimate, when you compare to the model based solely on the simplified impact_TOTAL.
#This model also is nearly identical in mean residuals
modelOutputPlots(chosen_model)
#there is a bit of an issue in the combined adjusted quantile test here for this model, which didn't rise to the level of an issue when just using the impact_TOTAL dataset.

source("STEP_6_morans_i.R")
morans_i_test(chosen_model, df)

#let's see what happens with live adult girdle, our other response of interest, when we keep using this more complex way of calculating impact coefficient
response <- df_wo_boundaries$live_adult_girdle
modelSelection(df_wo_boundaries$live_adult_girdle)
#same model as well: bitternut_hickory+american_elm+`alternate-leaf_dogwood` fits the data best (K = 5, AICc = -8.93)
source("STEP_5_testing_the_model.R")
chosen_model <- glmmTMB(adult_percent_live_canopy ~ bitternut_hickory+american_elm+`alternate-leaf_dogwood`, data = df_wo_boundaries, family = beta_family(link = "logit"))
modelTest(chosen_model)
# [[1]]
# Family: beta  ( logit )
# Formula:          adult_percent_live_canopy ~ bitternut_hickory + american_elm +      `alternate-leaf_dogwood`
# Data: df_wo_boundaries
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# -19.6      -9.6      14.8     -29.6        50 
# 
# 
# Dispersion parameter for beta family (): 1.43 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)               0.68003    0.23190   2.933  0.00336 **
#   bitternut_hickory        -0.01247    0.03456  -0.361  0.71837   
# american_elm             -0.01433    0.02392  -0.599  0.54897   
# `alternate-leaf_dogwood` -0.34803    0.14559  -2.390  0.01683 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# [[2]]
# [1] 0.0350625
#very similar again! it seems like using impact_height_direction_TOTAL instead of impact_TOTAL has no real effect on the fit of our models
modelOutputPlots(chosen_model)
#looking good!
#because of these results, it is probably not even worth it to use the more complex impact coefficient, as it doesn't change our findings in any meaningful way

source("STEP_6_morans_i.R")
morans_i_test(chosen_model, df)

#--------------------------------
#RELATIVIZING IMPACT COEFFICIENTS
df <- impact_PERCENT
source("STEP_3_initial_data_familiarization.R")
#similar spread of data and correlations found, which, again, isn't surprising

#Now it is time to determine a model
response <- df_wo_boundaries$adult_percent_live_canopy
source("STEP_4_determining_a_model.R")
modelSelection(df_wo_boundaries$adult_percent_live_canopy)
#Grand_Total is meaningless here, because it will always add up to 100%. so every model with Grand_Total in it will be dropped, and that's okay
#Model chosen is the same as the one chosen when using plain impact_TOTAL, based only on `alternate-leaf_dogwood` impact (K = 3 AICc = -24.88)

source("STEP_5_testing_the_model.R")
chosen_model <- glmmTMB(adult_percent_live_canopy ~ `alternate-leaf_dogwood`, data = df_wo_boundaries, family = beta_family(link = "logit"))
modelTest(chosen_model)
# [[1]]
# Family: beta  ( logit )
# Formula:          adult_percent_live_canopy ~ `alternate-leaf_dogwood`
# Data: df_wo_boundaries
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# -25.4     -19.3      15.7     -31.4        52 
# 
# 
# Dispersion parameter for beta family (): 1.48 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                0.6834     0.1819   3.758 0.000172 ***
#   `alternate-leaf_dogwood` -20.7937     7.9450  -2.617 0.008865 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# [[2]]
# [1] 0.0346243
#Again, what a similar performing model!!!
modelOutputPlots(chosen_model)
#again, quantile deviations detected, but overall, not an issue
source("STEP_6_morans_i.R")
morans_i_test(chosen_model, df)
# Moran I test under randomisation
# 
# data:  residuals_regression  
# weights: lw    
# 
# Moran I statistic standard deviate = 2.0438, p-value = 0.02049
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# 0.131685747      -0.018518519       0.005401133 

response <- df_wo_boundaries$live_adult_girdle
modelSelection(df_wo_boundaries$live_adult_girdle)
#same model as well: bitternut_hickory+american_elm+`alternate-leaf_dogwood` fits the data best (K = 5, AICc = -8.55)
source("STEP_5_testing_the_model.R")
chosen_model <- glmmTMB(live_adult_girdle ~ bitternut_hickory+american_elm+`alternate-leaf_dogwood`, data = df_wo_boundaries, family = beta_family(link = "logit"))
modelTest(chosen_model)
# [[1]]
# Family: beta  ( logit )
# Formula:          live_adult_girdle ~ bitternut_hickory + american_elm + `alternate-leaf_dogwood`
# Data: df_wo_boundaries
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# -9.8       0.3       9.9     -19.8        50 
# 
# 
# Dispersion parameter for beta family (): 3.36 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)               -0.5113     0.1840  -2.779  0.00546 **
#   bitternut_hickory          7.7777     3.3514   2.321  0.02030 * 
#   american_elm               8.0255     2.6673   3.009  0.00262 **
#   `alternate-leaf_dogwood`  14.7923     7.8229   1.891  0.05864 . 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# [[2]]
# [1] 0.003103197
modelOutputPlots(chosen_model)
#again, there are issues in quantile deviations

source("STEP_6_morans_i.R")
morans_i_test(chosen_model, df)
#NONSIGNIFICANT in autocorrelation
# Moran I test under randomisation
# 
# data:  residuals_regression  
# weights: lw    
# 
# Moran I statistic standard deviate = -0.50236, p-value = 0.6923
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# -0.055434486      -0.018518519       0.005400133

#--------------------------------
#ACCOUNTING ONLY FOR LIVING TREES
df <- impact_alive_TOTAL
source("STEP_3_initial_data_familiarization.R")

response <- df_wo_boundaries$adult_percent_live_canopy
source("STEP_4_determining_a_model.R")
modelSelection(df_wo_boundaries$adult_percent_live_canopy)
#`alternate-leaf_dogwood` remains victorious! (K = 3, AICc = -26.49)
source("STEP_5_testing_the_model.R")
chosen_model <- glmmTMB(adult_percent_live_canopy ~ `alternate-leaf_dogwood`, data = df_wo_boundaries, family = beta_family(link = "logit"))
modelTest(chosen_model)
# [[1]]
# Family: beta  ( logit )
# Formula:          adult_percent_live_canopy ~ `alternate-leaf_dogwood`
# Data: df_wo_boundaries
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# -27.0     -20.9      16.5     -33.0        52 
# 
# 
# Dispersion parameter for beta family (): 1.53 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)               0.70028    0.17970   3.897 9.74e-05 ***
#   `alternate-leaf_dogwood` -0.22259    0.07727  -2.881  0.00397 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# [[2]]
# [1] 0.03073619
modelOutputPlots(chosen_model)

source("STEP_6_morans_i.R")
morans_i_test(chosen_model, df)
# Moran I test under randomisation
# 
# data:  residuals_regression  
# weights: lw    
# 
# Moran I statistic standard deviate = 2.1731, p-value = 0.01489
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# 0.14102901       -0.01851852        0.00539048 

response <- df_wo_boundaries$live_adult_girdle
modelSelection(df_wo_boundaries$live_adult_girdle)
#you guessed it, bitternut_hickory+american_elm+`alternate-leaf_dogwood` wins, with K = 5 and AICc = -7.69
source("STEP_5_testing_the_model.R")
chosen_model <- glmmTMB(live_adult_girdle ~ bitternut_hickory+american_elm+`alternate-leaf_dogwood`, data = df_wo_boundaries, family = beta_family(link = "logit"))
modelTest(chosen_model)
# [[1]]
# Family: beta  ( logit )
# Formula:          live_adult_girdle ~ bitternut_hickory + american_elm + `alternate-leaf_dogwood`
# Data: df_wo_boundaries
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# -9.8       0.3       9.9     -19.8        50 
# 
# 
# Dispersion parameter for beta family (): 3.36 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)               -0.5113     0.1840  -2.779  0.00546 **
#   bitternut_hickory          7.7777     3.3514   2.321  0.02030 * 
#   american_elm               8.0255     2.6673   3.009  0.00262 **
#   `alternate-leaf_dogwood`  14.7923     7.8229   1.891  0.05864 . 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# [[2]]
# [1] 0.003103197
modelOutputPlots(chosen_model)

source("STEP_6_morans_i.R")
morans_i_test(chosen_model, df)
# Moran I test under randomisation
# 
# data:  residuals_regression  
# weights: lw    
# 
# Moran I statistic standard deviate = -0.50236, p-value = 0.6923
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# -0.055434486      -0.018518519       0.005400133 

#------------------------------
#ONLY LIVING TREES, RELATIVIZED
df <- impact_alive_PERCENT
source("STEP_3_initial_data_familiarization.R")

response <- df_wo_boundaries$adult_percent_live_canopy
source("STEP_4_determining_a_model.R")
modelSelection(df_wo_boundaries$adult_percent_live_canopy)
#`alternate-leaf_dogwood`*white_cedar (K = 5, AICc = -25.57) is now the best model, different from each other dataset, but by a very slight margin, compared to just `alternate-leaf_dogwood` (K = 3, AICc = -25.13)
source("STEP_5_testing_the_model.R")
chosen_model <- glmmTMB(adult_percent_live_canopy ~ `alternate-leaf_dogwood`*white_cedar, data = df_wo_boundaries, family = beta_family(link = "logit"))
modelTest(chosen_model)
# [[1]]
# Family: beta  ( logit )
# Formula:          adult_percent_live_canopy ~ `alternate-leaf_dogwood` * white_cedar
# Data: df_wo_boundaries
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# -26.8     -16.8      18.4     -36.8        50 
# 
# 
# Dispersion parameter for beta family (): 1.65 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                            0.7928     0.2116   3.746 0.000179 ***
#   `alternate-leaf_dogwood`             -34.8049    10.4782  -3.322 0.000895 ***
#   white_cedar                           -1.6134     1.6455  -0.980 0.326842    
# `alternate-leaf_dogwood`:white_cedar 117.8321    58.4573   2.016 0.043832 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# [[2]]
# [1] 0.02901244
modelOutputPlots(chosen_model)
#all good here though, thank goodness.

source("STEP_6_morans_i.R")
# morans_i_test(chosen_model, df)
# Moran I test under randomisation
# 
# data:  residuals_regression  
# weights: lw    
# 
# Moran I statistic standard deviate = -0.50236, p-value = 0.6923
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# -0.055434486      -0.018518519       0.005400133 

response <- df_wo_boundaries$live_adult_girdle
modelSelection(df_wo_boundaries$live_adult_girdle)
#still, bitternut_hickory+american_elm+`alternate-leaf_dogwood` wins, with K = 5 and AICc = -7.33
source("STEP_5_testing_the_model.R")
chosen_model <- glmmTMB(live_adult_girdle ~ bitternut_hickory+american_elm+`alternate-leaf_dogwood`, data = df_wo_boundaries, family = beta_family(link = "logit"))
modelTest(chosen_model)
# [[1]]
# Family: beta  ( logit )
# Formula:          live_adult_girdle ~ bitternut_hickory + american_elm + `alternate-leaf_dogwood`
# Data: df_wo_boundaries
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# -8.6       1.5       9.3     -18.6        50 
# 
# 
# Dispersion parameter for beta family (): 3.29 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)               -0.4551     0.1785  -2.550  0.01077 * 
#   bitternut_hickory          6.0602     3.1053   1.952  0.05098 . 
# american_elm               7.1538     2.4066   2.973  0.00295 **
#   `alternate-leaf_dogwood`  13.1402     6.8354   1.922  0.05456 . 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# [[2]]
# [1] 0.003399762
modelOutputPlots(chosen_model)
#all good!!

source("STEP_6_morans_i.R")
morans_i_test(chosen_model, df)
# Moran I test under randomisation
# 
# data:  residuals_regression  
# weights: lw    
# 
# Moran I statistic standard deviate = -0.56224, p-value = 0.713
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# -0.060038287      -0.018518519       0.005453391 