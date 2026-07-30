#Catherine dell'Olio
#This code is intended to take the whole of our health assessment data and split it into more workable pieces, so that we can make meaningful comparisons.

#installing dplyr for ease of use working with large datasets
install.packages("dplyr")
library(dplyr)
#setting working directory
setwd("~/dataAnalysisCfd")
#data will be the whole of our data
data <- read.csv("datasheet.csv")

#removing CPVT29, where DBH is two numbers listed for two trunks.  This is an issue because R cannot see that value as a number, and all values in the same column must be of the same type, so including this coerces every other DBH to not be an number too.
data <- data[c(1:6, 8:429),]
#now that that datapoint is gone, I can convert all those numbers back to numbers.
data$dbh_cm <- as.numeric(data$dbh_cm)

#separating all the seedlings out
data_seedlings <- data %>% filter(seedling_or_adult == "Seedling")
#all the seedlings have specific variables to create, listed below
average_densi <- numeric(nrow(data_seedlings))
mowing <- logical(nrow(data_seedlings))
vole_nibbling <- logical(nrow(data_seedlings))
deer_browse <- logical(nrow(data_seedlings))
deer_rub <- logical(nrow(data_seedlings))
rabbit_browse <- logical(nrow(data_seedlings))
data_seedlings <- cbind(data_seedlings, average_densi, mowing, vole_nibbling, deer_browse, deer_rub, rabbit_browse)
#we will now fill in these extra variables
#for average densiometer, we will average the four densiometer readings and multiply by 1.04
for (i in 1:nrow(data_seedlings)) {
  data_seedlings$average_densi[i] <- 1.04 * mean(c(data_seedlings$densi_n[i], data_seedlings$densi_s[i], data_seedlings$densi_e[i], data_seedlings$densi_w[i]))
}
#all of the types of damage are combined in the seedling_damage category, and so each combination is treated as a different thing.  I am going to extract out each type of damage into different variables.
data_seedlings[grep('Having been mowed in the past', data_seedlings$seedling_damage), "mowing"] <- TRUE
data_seedlings[grep("Vole (nibbling on the bark)", data_seedlings$seedling_damage), "vole_nibbling"] <- TRUE
data_seedlings[grep("Deer browse", data_seedlings$seedling_damage), "deer_browse"] <- TRUE
data_seedlings[grep("Deer rub", data_seedlings$seedling_damage), "deer_rub"] <- TRUE
data_seedlings[grep("Rabbit browse", data_seedlings$seedling_damage), "rabbit_browse"] <- TRUE

#now separating out all the adults
data_adults <- data %>% filter(seedling_or_adult == "Adult")

#separating out all of the living seedlings out of the seedlings data, so that they retain all of the addtional variables added in
data_liveseedlings <- data_seedlings %>% filter(seedling_dead_or_alive == "Alive")

#creating different datasheets for each year's set of seedlings.  In all honesty, I thought this would be more useful then it turned out, so I am going to leave it commented out for now.
#for (i in 2017:2025){
#  nam <- paste("data_seedlings", i, sep = "")
#  assign(nam, data_seedlings %>% filter(seedling_germ_yr == {{i}}))
#}

#separating out the live adults from just the adults
data_liveadults <- data_adults %>% filter(adult_dead_or_alive == "Alive")

#for the living adults or seedlings, if we found no canker on them, we would not be directed to the health assessment questions about extent of canker, so those would be left with NA.  NA values make sense for dead adults or seedlings, but living adults can safely be considered a 0 on canker metrics if they do not have canker reported.

#this function converts na to zero for any variable and return it to the same data frame
na.zero <- function (x) {
  x[is.na(x)] <- 0
  x
}

#I am converting all NAs to 0 for metrics of canker on live adults and seedlings
data_liveadults$live_adult_canker_base <- na.zero(data_liveadults$live_adult_canker_base)
data_liveadults$live_adult_canker_1stlivebranch <- na.zero(data_liveadults$live_adult_canker_1stlivebranch)
data_liveadults$live_adult_girdle <- na.zero(data_liveadults$live_adult_girdle)
data_liveseedlings$seedling_canker_base <- na.zero(data_liveseedlings$seedling_canker_base)
data_liveseedlings$seedling_canker_stem <- na.zero(data_liveseedlings$seedling_canker_stem)
data_liveseedlings$seedling_girdle <- na.zero(data_liveseedlings$seedling_girdle)
for (i in 1:nrow(data_liveadults)){
  if (data_liveadults$live_adult_canker_presence[i] == "N" && is.na(data_liveadults$live_adult_canker_1stlivebranch[i]))
    data_liveadults$live_adult_canker_1stlivebranch[i] <- 0
}

#for live adults, I am also going to use the retention strategy based on Farley 2010 as a metric of health, in which butternuts with at least 70% live canopy and less than 20% canker were considered healthy.  Alternatively, tree with at least 50% live canopy but no canker are also considered healthy.
#creating a variable to record retention strategy
retention_strategy <- numeric(nrow(data_liveadults))
data_liveadults <- cbind(data_liveadults, retention_strategy)

for(i in 1:nrow(data_liveadults)){
  if ((data_liveadults$adult_percent_live_canopy[i] > 70  &
      data_liveadults$live_adult_canker_base[i] < 20) |
      (data_liveadults$adult_percent_live_canopy[i] >= 50 &
      data_liveadults$live_adult_canker_base[i] == 0))
      data_liveadults$retention_strategy[i] <- "healthy"
      else data_liveadults$retention_strategy[i] <- "unhealthy"
}

data_healthyadults <- data_liveadults %>% filter(retention_strategy == "healthy")


#these next few split up our data by site, so that we can see if patterns across the two sites remain consistent when looking at either site or the other.
data_ILM <- data %>% filter(site == "ILM")
data_CPVT <- data %>% filter(site == "CPVT")

data_CPVTadults <- data_adults %>% filter(site == "CPVT")
data_ILMadults <- data_adults %>% filter(site == "ILM")

data_CPVTseedlings <- data_seedlings %>% filter(site == "CPVT")
data_CPVTliveseedlings <- data_liveseedlings %>% filter(site == "CPVT")

data_ILMseedlings <- data_seedlings %>% filter(site == "ILM")
data_ILMliveseedlings <- data_liveseedlings %>% filter(site == "ILM")

data_CPVTliveadults <- data_liveadults %>% filter(site == "CPVT")
data_ILMliveadults <- data_liveadults %>% filter(site == "ILM")