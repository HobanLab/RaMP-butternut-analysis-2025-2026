#Catherine dell'Olio
# This code allows me to select a stratified sample from the ILM adult butternuts based on health (both % girdle and % live canopy, in a 3x3 subgroup stratification for high, medium, and low health in both metrics independently.)
#After I generated a sample from this code, I then tested its representativeness of the pool of eligible adult butternuts in terms of spatial extent and crown class.
#I call the package dpylr to filter my data with its syntax
#I call one datasheet: ILM_assoc_trees, which is just a datasheet of all the ILM trees.  Alternatively, this can be run with the full datasheet.csv filtered by site == "ILM".

#setting working directory
set.seed(5005)
setwd("~/dataviz_cfd")
df <- read.csv("ILM_assoc_trees.csv")
library(dplyr)

#filter to eligible butternuts: adults with DBH greater than 4.5 cm and less than 50 cm
adults_df <- df %>% filter(seedling_or_adult == "Adult")
adults_df <- adults_df %>% filter(number != "212" & number != "314" & number != "309" & number != "428" & number != "387" & number != "041" & number != "558" & number != "001" & number != "544")
adults_df <- adults_df %>% filter(dbh_cm >= 4.5)
adults_df <- adults_df %>% filter(dbh_cm <=50)
#excluding trees in agricultural lands
adults_df <- adults_df %>% filter(number != "193" & number != "424" & number != "425" & number != 119 & number != 120)

#creating the subgroups:
#health is represented by two metrics: girdle percentage and percent live canopy
# girdle is considered high if it is over 60%, low if it is less than or equal to 30%, and medium if in between
# percent live canopy is considered high if it is over 70%, low if it is less than or equal to 30%, and medium if in between
# these are independent metrics so they create a 3x3 grid of possible subgroups
cat1 <- adults_df %>%  filter(adult_percent_live_canopy <= 70 & live_adult_girdle <= 50) #low canopy low girdle
cat2 <- adults_df %>%  filter(adult_percent_live_canopy <= 70 & live_adult_girdle > 50) #low canopy high girdle
cat3 <- adults_df %>%  filter(adult_percent_live_canopy > 70 & live_adult_girdle > 50) #high canopy low girdle
cat4 <- adults_df %>%  filter(adult_percent_live_canopy > 70 & live_adult_girdle <= 50) #high canopy high girdle
scatterplot(adults_df, adult_percent_live_canopy, live_adult_girdle, "df")

ggplot() +
  geom_point(data = cat1, aes(gps_w, gps_n), colour="green2") +
  geom_point(data = cat2, aes(gps_w, gps_n), colour="red") +
  geom_point(data = cat3, aes(gps_w, gps_n), colour="yellow") +
  geom_point(data = cat4, aes(gps_w, gps_n), colour="blue") +
  ggtitle("Map of all eligible adults by category")
#blue (high canopy and low girdle) conains a lot of the extremities and also many in middle

#a function to sample a certain size from each category.  the reason I don't have it sample 6 from every category is because category 1 is too small--only 2 trees.  So those will both be sampled.
randomsample <- function(cat, size){
  s <-  sample(1:nrow({{cat}}), size, replace = FALSE)
  return(cat[s,])
}

#we need to figure out a way to have the computer only select trees that are more than 20 meters apart (2 radii)
#shout out rodrigo on stack overflow for this function: https://stackoverflow.com/questions/639695/how-to-convert-latitude-or-longitude-to-meters
measure <- function(lon1,lat1,lon2,lat2) {
  R <- 6378.137                                # radius of earth in Km
  dLat <- (lat2-lat1)*pi/180
  dLon <- (lon2-lon1)*pi/180
  a <- (sin(dLat/2))^2 + cos(lat1*pi/180)*cos(lat2*pi/180)*(sin(dLon/2))^2
  c <- 2 * atan2(sqrt(a), sqrt(1-a))
  d <- R * c
  return (d * 1000)                            # distance in meters
}

#RUN THIS ONCE. It resets the sample back to nothing.
totalsample <- data.frame()
adults_remaining <- adults_df

#RUN FROM HERE TO DASHED LINE AS MANY TIMES AS NEEDED TO CREATE A SAMPLE OF THE DESIRED SIZE (in my case, 13 times to get a sample of 52).  This builds the sample by selecting one tree from each category, checking to exclude any very nearby trees to the selected one that now would be spatially autocorrelated with it, and filtering the remaining trees to only the ones that are still eligible (that is, not too close to any of the previously selected trees).
set.seed(5005)
t1 <- randomsample(cat3, 1)
adults_cycle <- data.frame()
for (i in 1:nrow(adults_remaining)){
  if (measure(t1$gps_w, t1$gps_n, adults_remaining[i,]$gps_w, adults_remaining[i,]$gps_n) > 20){
    adults_cycle <- rbind(adults_cycle, adults_remaining[i,])
  }
}
adults_remaining <- adults_cycle
totalsample <- rbind(totalsample, t1)
set.seed(5005)
t1 <- randomsample(cat1, 1)
adults_cycle <- data.frame()
for (i in 1:nrow(adults_remaining)){
  if (measure(t1$gps_w, t1$gps_n, adults_remaining[i,]$gps_w, adults_remaining[i,]$gps_n) > 20){
    adults_cycle <- rbind(adults_cycle, adults_remaining[i,])
  }
}
adults_remaining <- adults_cycle
totalsample <- rbind(totalsample, t1)
set.seed(5005)
t1 <- randomsample(cat2, 1)
adults_cycle <- data.frame()
for (i in 1:nrow(adults_remaining)){
  if (measure(t1$gps_w, t1$gps_n, adults_remaining[i,]$gps_w, adults_remaining[i,]$gps_n) > 20){
    adults_cycle <- rbind(adults_cycle, adults_remaining[i,])
  }
}
adults_remaining <- adults_cycle
totalsample <- rbind(totalsample, t1)
set.seed(5005)
t1 <- randomsample(cat4, 1)
adults_cycle <- data.frame()
for (i in 1:nrow(adults_remaining)){
  if (measure(t1$gps_w, t1$gps_n, adults_remaining[i,]$gps_w, adults_remaining[i,]$gps_n) > 20){
    adults_cycle <- rbind(adults_cycle, adults_remaining[i,])
  }
}
adults_remaining <- adults_cycle
totalsample <- rbind(totalsample, t1)
ggplot() +
  geom_point(data = adults_df, aes(gps_w, gps_n), colour="#FDE000") +
  geom_point(data = adults_remaining, aes(gps_w, gps_n), colour="#3B528B") +
  geom_point(data = totalsample, aes(gps_w, gps_n), colour="#5EC962") +
  ggtitle("Remaining eligible adults (yellow) with sample (green) and excluded autocorrelated trees (blue)")

cat1 <- adults_remaining %>%  filter(adult_percent_live_canopy <= 70 & live_adult_girdle <= 50) #low canopy low girdle (dying of other causes)
cat2 <- adults_remaining %>%  filter(adult_percent_live_canopy <= 70 & live_adult_girdle > 50) #low canopy high girdle (dying of canker)
cat3 <- adults_remaining %>%  filter(adult_percent_live_canopy > 70 & live_adult_girdle > 50) #high canopy high girdle (surviving despite canker)
cat4 <- adults_remaining %>%  filter(adult_percent_live_canopy > 70 & live_adult_girdle <= 50) # high canopy low girdle (spared from canker)
#---------------------------------------------------------------------------
#END OF CHUNK TO RUN MULTIPLE TIMES. Run as one big block 13 times to recreate my sample (the set seed will ensure you create the same sample every time!)

#retrieve tree IDs for the proposed sample
totalsample$number
#This gives me the list of trees to sample at ILM for the associate trees analysis!

#this creates a matrix of the distances between trees to ensure we aren't getting an autocorrelated sample.
distances <- data.frame()
for (i in 1:52){
  for (j in 1:52){
    distances[i,j] <- measure(totalsample[i,4], totalsample[i,3], totalsample[j,4], totalsample[j,3])
  }
}

#------------------------------
#TESTING FOR REPRESENTATIVENESS
#------------------------------
#does this sample consist of mostly similar proportions of trees in each crown class as the eligible population?  Trees of different crown classes receive different amounts of light which might really impact our analysis.
barplot <- function(data, variable, title){
  ggplot(data, aes({{variable}})) +
    geom_bar()+
    ggtitle(title)
}
barplot(totalsample, adult_crown_class, "crown class of proposed sample")
barplot(adults_df, adult_crown_class, "crown class of total pool of eligible adults")
#the answer is definitely yes
#does the sample appear to cover the spatial extent of eligible trees?
ggplot() +
  geom_point(data = adults_df, aes(gps_w, gps_n), colour="green2") +
  geom_point(data = totalsample, aes(gps_w, gps_n), colour="red") +
  ggtitle("Map of all eligible adults (green) with proposed sample (red) overlaid")

#produce a csv of which trees to sample
write.csv(totalsample, "proposed_sample_ILM_assoc_trees.csv")
