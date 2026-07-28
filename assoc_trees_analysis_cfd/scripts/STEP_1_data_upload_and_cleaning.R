#Catherine dell'Olio
#datacleaning for associate trees analysis
#step 1: data upload and cleaning

#You will need tidyverse and dplyr to run this, 

library(tidyverse)
library(dplyr)

#--------------------------------------------
#Step zero
#DATA CLEANING BEFORE MOVING TO R (IN SHEETS)
#all of the data collected: distance to focal tree (m),	direction (degree 0-360),	height (Ft),	and DBH (cm)	
#for two areas, we also collected canopy max (m) and	canopy min (m) in the field
#calculate impact coefficient (=DBH/distance)
#BONUS then calculate basal area (=(DBH^2)*0.00007854)
#sum up the data collected by species in a pivot table
#copy over pivot table values to a new spreadsheet that has each focal tree's health status

#then I create a SECOND datasheet that doesn't include dead trees, but is otherwise exactly the same.  The easiest way to do this is to copy the original datasheet and filter out dead trees, then repeat all the steps above.

#once this is done in Sheets, all of the data sheets I use are downloaded onto this repository as CSVs, and it's time to get moving in R!
#---------------------
#Step one
#DATA UPLOAD
#make sure to set the working directory first!

setwd("~/assoc_trees_analysis_cfd/datasheets")
#abundance_DA (dead or alive) simply counts the number of each species around each focal butternut, and includes counts of both dead and alive trees.
abundance_DA <- read.csv("26_04_associate_tree_data - abundance_data_DA.csv")
#impact_TOTAL sums up the impact coefficient (DBH divided by distance) for each species, around each of the focal butternuts, combined with the 2025 health assessment data from that tree
impact_TOTAL <- read.csv("26_04_associate_tree_data - impact_TOTAL.csv")
#impact_PERCENT relativizes the impacts of all tree species
impact_PERCENT <- read.csv("26_04_associate_tree_data - impact_PERCENT.csv")
#I also experimented with a slightly more complex version of impact coefficient that accounts for the relative size of the associate (multiplying the impact coefficient by the ratio of the height of the associate to the height of the focal butternut---in essence this boosts the impact if the associate is taller and lowers it if the associate is shorter than the focal butternut) and also accounts for the direction of the associate tree (where southernly associates are boosted in impact coefficient by 1.5 to account for the extra light they block).4
impact_height_direction_TOTAL <- read.csv("26_04_associate_tree_data - impact_height_direction_TOTAL.csv")

#I ALSO ran both the impact_TOTAL and impact_PERCENT based models again with only currently living trees included in the impact coefficients of each species and overall, in order to see if the dead trees were causing the models to change significantly.
impact_alive_TOTAL <- read.csv("26_04_associate_tree_data - impact_alive_TOTAL.csv")
impact_alive_PERCENT <- read.csv("26_04_associate_tree_data - impact_alive_PERCENT.csv")

#Abundance_DA is a smaller dataframe than each other so it gets a separate cleaning process!
#deleting spaces from what will become the column names
# Source - https://stackoverflow.com/a/46421817
# Posted by ZWL, modified by community. See post 'Timeline' for change history
# Retrieved 2026-04-13, License - CC BY-SA 4.0
for (i in 1:ncol(abundance_DA)){
  abundance_DA[1,i] <- gsub(" ", "_", abundance_DA[1, i]) 
}
#remove duplicate labels in second row
# Source - https://stackoverflow.com/a/49634498
# Posted by Shalini Baranwal
# Retrieved 2026-04-13, License - CC BY-SA 3.0
names(abundance_DA) <- abundance_DA[1,]
abundance_DA <- abundance_DA[-1,]
#trim off the grand total column/row
abundance_DA <- abundance_DA[-57, -35]

#replace all blank cells with a 0
abundance_DA[abundance_DA == ''] <- 0
abundance_DA <- abundance_DA %>% mutate_at(-c(1), as.numeric)

#CLEANING EACH OTHER DATA SHEET (full health data plus impact coefficients!)
clean_data <- function(uncleaned_df){
  #R DATA CLEANING (run this on again whenever you reset the current 'uncleaned_df' as above)
  #deleting spaces from what will become the column names
  # Source - https://stackoverflow.com/a/46421817
  # Posted by ZWL, modified by community. See post 'Timeline' for change history
  # Retrieved 2026-04-13, License - CC BY-SA 4.0
  for (i in 1:ncol(uncleaned_df)){
    uncleaned_df[1,i] <- gsub(" ", "_", uncleaned_df[1, i]) 
  }
  #remove duplicate labels in second row
  # Source - https://stackoverflow.com/a/49634498
  # Posted by Shalini Baranwal
  # Retrieved 2026-04-13, License - CC BY-SA 3.0
  names(uncleaned_df) <- uncleaned_df[1,]
  uncleaned_df <- uncleaned_df[-1,]
  #trim off column labels in last row, as well as associate tree data at the end of the health assessment (which has duplicate names with the new data)
  uncleaned_df <- uncleaned_df[1:55,c(1:95)]
  
  #replace all blank cells with a 0
  # Source - https://stackoverflow.com/a/67404924
  # no longer works so I need to make a new version
  uncleaned_df[uncleaned_df == ''] <- 0
  uncleaned_df <- uncleaned_df %>% mutate_at(c(2:37, 39:41, 43, 45, 56, 61:66), as.numeric)
  #convert percentages into actual proportions
  for (i in 1:nrow(uncleaned_df)){
    uncleaned_df$adult_percent_live_canopy[i] <- uncleaned_df$adult_percent_live_canopy[i]/100
    uncleaned_df$live_adult_girdle[i] <- uncleaned_df$live_adult_girdle[i]/100
    uncleaned_df$live_adult_canker_1stlivebranch[i] <- uncleaned_df$live_adult_canker_1stlivebranch[i]/100
    uncleaned_df$live_adult_canker_base[i] <- uncleaned_df$live_adult_canker_base[i]/100 #%%%FIX%%%for some inexplicable reason, this didn't work for impact_alive_TOTAL but works for everything else.
    #I am also going to add the values of each individual ash species to our ash_spp sum, as we will be summing the impact of all ashes together
    uncleaned_df$ash_spp[i] <- (uncleaned_df$ash_spp[i] + uncleaned_df$black_ash[i] + uncleaned_df$green_ash[i] + uncleaned_df$white_ash[i])
  }
  return(uncleaned_df)
}

impact_TOTAL <- clean_data(impact_TOTAL)
impact_PERCENT <- clean_data(impact_PERCENT)
impact_height_direction_TOTAL <- clean_data(impact_height_direction_TOTAL)
impact_alive_TOTAL <- clean_data(impact_alive_TOTAL)
impact_alive_PERCENT <- clean_data(impact_alive_PERCENT)
