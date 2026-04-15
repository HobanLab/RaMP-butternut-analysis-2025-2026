#Catherine dell'Olio
#dataAnalysis, associate trees experiment

#This code will help me analyze my associate trees protocol data, collected in March/April fieldwork trip.  In this trip, I am collecting information about the associate trees around a subsample of butternuts.  These trees were selected to represent the wide range of health statuses of butternuts at this fieldsite.  I think of these as belonging to four distinct groups: those that are healthty with little evidence of canker, those that are very healthy despite a lot of canker, those dying but little canker, and those dying with a lot of canker.  I am recording all of the trees that occur within a 10 meter radius around each butternut selected---their height, species, distance to the focal tree, DBH, and canopy width.

#I can use the DBH and distance to the focal tree for all individuals of a given species to calculate an impact coefficient of each species around the butternut, and determine the impacts (both scale and direction) for species around the butternut.

#The basic steps of this analysis process are:
#0 (because it is done BEFORE this script is run): Clean up data in google sheets
#1: load packages and data (and all data subsets)
#+++++++++++++++++++++++++++++++++++++++++++++


#----------------------------
#Step zero
#DATA CLEANING OFF R (SHEETS)
#----------------------------
#data should look like: https://docs.google.com/spreadsheets/d/1uJCOvdRocM-r9MQwYGmdiycVxVrh9uFi29-qGX8h0sk/edit?gid=900838604#gid=900838604
#you will need to collect focal tree	species	distance (m)	direction	height (Ft)	DBH (cm)	canopy max (m)	canopy min (m) in the field
#calculate impact coefficient (=DBH/distance) #BONUS then calculate basal area (=(DBH^2)*0.00007854)
#sum these by species in a pivot table
#copy over pivot table values to a new spreadsheet that has each focal tree's health status
#upload data into "~/dataAnalysisCfd/april_assoc_trees_analysis"

#create a SECOND datasheet that doesn't include dead trees, but is otherwise exactly the same.  The easiest way to do this is to copy the original datasheet and filter out dead trees, then repeat all the steps above.

assoc_data <- read.csv("26_04_associate_tree_data.csv")
assoc_data_larger_trees <- assoc_data %>% filter(DBH..cm. >= 5)

assoc_data_larger_trees_only_living <- assoc_data_larger_trees %>% filter(dead.alive. == "alive")

write.csv(assoc_data_larger_trees, "Assoc_data_larger_trees.csv")

#----------------------------------------
#Step one
#DATA + PACKAGE UPLOAD, DATA VERIFICATION
#----------------------------------------
library(tidyverse)
library(dplyr)
library(spdep)
library(gam)
library(AICcmodavg)
library(lme4)
library(RVAideMemoire) # for function "overdisp.glmer" to test overdispersion
library(MuMIn) # for the r.squaredGLMM
library(DHARMa) # for a simulation based approach to inspection of residuals for model diagnostics 
library(pscl) #zero inflated
library(tidyverse)
library(lsmeans)
library(multcompView)
library(multcomp)
library(plotrix) #std err
library(glmmTMB) #zeroinfl poisson (necessary?[])

#upload data
setwd("~/dataAnalysisCfd/april_assoc_trees_analysis")
dummy_data <- read.csv("dummy_data.csv")
#IMPORTANT: I am changing the datasheet name from a slightly more informative one I just defined in the R environment so that I can see I'm working with the dummy data right now, but don't have to change every single instance of its name when I do this with real data. This will also come in handy when I rerun the model taking out all dead trees---I will only need to add in
#dummy_data_wo_dead_trees <- read.csv("fake_data_wo_dead_trees.csv")
#df <- dummy_data_wo_dead_trees
#instead of changing every instance of <data> to <data_without_dead_trees>
df <- dummy_data

#I lied: there is a small amount of data cleaning that is easier to do in R because it is about making the data nice to play with in R (Ensuring header rows are recognized as such, that data are categorized as numeric, etc.)
#R DATA CHECKING
#remove duplicate labels in second row
df <- df[-c(1),]
#convert percentages into actual proportions
for (i in 1:nrow(df)){
  df$adult_percent_live_canopy[i] <- df$adult_percent_live_canopy[i]/100
  df$live_adult_girdle[i] <- df$live_adult_girdle[i]/100
  df$live_adult_canker_1stlivebranch[i] <- df$live_adult_canker_1stlivebranch[i]/100
}

#coerce all impact coefficients to be numeric, because clearly they are
df$Grand.Total.impact.coeff <- as.numeric(df$Grand.Total.impact.coeff)
df$basswood.impact.coeff <- as.numeric(df$basswood.impact.coeff)
df$bitternut.hickory.impact.coeff <- as.numeric(df$bitternut.hickory.impact.coeff)
df$black.cherry.impact.coeff <- as.numeric(df$black.cherry.impact.coeff)
df$butternut.impact.coeff <- as.numeric(df$butternut.impact.coeff)
df$european.beech.impact.coeff <- as.numeric(df$european.beech.impact.coeff)
df$red.hickory.impact.coeff <- as.numeric(df$red.hickory.impact.coeff)
df$sugar.maple.impact.coeff <- as.numeric(df$sugar.maple.impact.coeff)
df$sycamore.impact.coeff <- as.numeric(df$sycamore.impact.coeff)
df$white.cedar.impact.coeff <- as.numeric(df$white.cedar.impact.coeff)
#stop: is this all of your tree species of interest?
str(df)
#CHECKPOINT ONE: do all the columns look to be the right kind of data?

#----------------------------
#Step two
#INITIAL DATA FAMILIARIZATION
#----------------------------

#let's visualize the distribution of our response variables (metrics of health)
hist(df$adult_percent_live_canopy)#notes
hist(df$live_adult_girdle)#notes
hist(df$live_adult_canker_1stlivebranch)#notes
hist(df$adult_purdue_canker_rating)#notes
hist(df$adult_purdue_canopy_ranking)#notes

#now we can start visualizing the response as related to our predictors
#fixed effects: crown class
#random effects: potentially many---each species can be an effect, as well as the grand total impact coefficient

#crown class is the only factor here with discrete levels, so it needs to be called a factor
df$adult_crown_class <- factor(df$adult_crown_class)

plot(adult_percent_live_canopy ~ adult_crown_class, data = df) #what does it look like?

plot(adult_percent_live_canopy ~ Grand.Total.impact.coeff,
     data = df,
     xlab = 'Summed impact coefficient',
     ylab = 'Percent live canopy')
#what does it look like?

#with continuous predictor variables comes testing of correlations
#[]+++++++++++++++++++Emma's code?

#testing out different models, and comparing them with AICc (correction due to small sample size)
modList <- list()

modList[["grand.total*crown_class"]] <- glmmPQL1 <- 
  glmmPQL(adult_percent_live_canopy ~ adult_crown_class, random = ~ 1 | Grand.Total.impact.coeff, 
          data = df, 
          family = quasibinomial(link = "logit"))

testmodel <- glmer(adult_percent_live_canopy ~ Grand.Total.impact.coeff | adult_crown_class, family = binomial(link = "logit"), data = df)

glm.comple.pooling <- glm(adult_percent_live_canopy ~ Grand.Total.impact.coeff, 
                          data = df, 
                          family = quasibinomial)
#*this is for a quasibinomial distribution.  It seems we have two paths according to https://stats.stackexchange.com/questions/333410/comparing-quasibinomial-glms-in-r
#*1) beta distribution using the betareg package, and then use stats::AIC to compare
#*2) use binomial/quasibinomial and then use the quasilikelihood AIC (which is a bit more precarious)

#recommended for determining if beta is what I want, but I don't really understand why as it doesn't seem you can compare it to other distributions easily
#util_beta_aic(df$adult_percent_live_canopy)

#-----------------------
#ORDERED BETA REGRESSION involves ordbetareg instead of this one
#-----------------------







#---------
#MORAN'S I
#---------
#Testing moran's i, which is our way of testing the impact of site on the residuals etc.
regression <- glm(adult_percent_live_canopy ~ Grand.Total.impact.coeff, data = df)
summary(regression)
structural_diff = glm(adult_percent_live_canopy ~ Grand.Total.impact.coeff, data = df)

#gam_model <- gam(adult_percent_live_canopy ~ s(Grand.Total.impact.coeff),
#                 family=quasibinomial(link="logit"), data = dummy_data, na.action = na.omit)
#summary(gam_model)

plot(regression, shade = TRUE, rug = TRUE) #this makes a few diagnostic plots that show you how well your model satisfies assumptions
#first should be a starry sky, straight line across 0.  if not you can transform the predictors but not response with the quasibinomial model I am using
# second should be a 1 to 1 relationship
#third should also be nonstructured (starry sky) but okay if the red line is weird as long as there is no structure to it
#last one tests for influential outliers, and if there is a point on or in the cook's distance lines it is a potential influential outlier

residuals_regression <- residuals(regression, type = "pearson")

coords <- df[, c("gps_w", "gps_n")]
nb <- knn2nb(knearneigh(coords, k = 5))  # Use 5 nearest neighbors
lw <- nb2listw(nb, style = "W")

#Compute Moran's I
moran_test <- moran.test(residuals_regression, lw)
print(moran_test)