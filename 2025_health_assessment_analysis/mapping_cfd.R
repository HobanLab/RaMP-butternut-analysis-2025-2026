#Catherine dell'Olio, cdellolio@mortonarb.org
#This code relies on the script dataSplitter.R, which in turn relies on the data which is collected in datasheet.csv, and the package dplyr for its syntax.  The purpose of dataSplitter.R is to filter and process my datasheet in a consistent way across different analyses (with the benefit of streamlining those analysis scripts).  Please see notes in that script for the explanation of what's going on within.
#I also call the script chartBuilder.R, which is again consistent across my analyses and creates a set of handy functions so that I can standardize plots (maps, boxplots, barplots, histograms, and scatterplots) and create them in a hurry.
#As implied by the name, for this analysis, I heavily rely on a function I call map.  Map uses ggplots to create a scatterplot of latitude and longitude of butternuts---essentially mapping their distribution by GPS points collected in the field.  If I am using a visualization technique that requires it, I have built in the option to include a legend for clarity.

#The purpose of this code is to map our butternuts and in doing so, show spatial patterns of age and health.

#setting working directory and running dataSplitter to get our data into readable and processable units.
setwd("~/dataAnalysisCfd")
source('dataSplitter.R')
source('chartBuilder.R')
#I can now work with relevant subsets of data.  Because this data spans two sites quite far from each other, mapping sites separately just makes sense for this.

#Here I have the map of all trees surveyed at each site, categorized by their age class (seedling or adult; adults are defined as having a DBH).
map(data_CPVT, seedling_or_adult, TRUE, "CPVT butternuts")
map(data_ILM, seedling_or_adult, TRUE, "ILM butternuts")


map(data_CPVTliveadults, retention_strategy, TRUE, "CPVT live adults")
map(data_ILMliveadults, retention_strategy, TRUE, "ILM live adults")

#We collected a lot of data on the health of these butternuts, and so we have a few different metrics that measure the similar markers of canker in different ways.  As such, it can be a bit confusing to determine how these disparate markers all relate to and match up to ambiguous terms like butternut "health" or "risk".  These two show a map of ILM live adults based on two different metrics of canker---the canker severity rating (adapted from Purdue's methodology), which ranks canker severity over the whole tree on a scale of 1-5 (5 being the tree basically dead of canker), and the percentage of canker area at the base of the trunk.
map(data_ILMliveadults, live_adult_purdue_canker_rating, TRUE, "ILM adults canker severity")
map(data_ILMliveadults, live_adult_canker_base, TRUE, "ILM adults canker severity")
#We can see at this site that the two metrics, owing to their high correlationm, reveal similar spatial patterns when mapped.

#I will now map deer browse and other damage spatially to check for patterns.
map(data_CPVTseedlings, deer_browse, TRUE, "Deer browse patterns at CPVT")
map(data_ILMseedlings, deer_browse, TRUE, "Deer browse patterns at ILM")
map(data_CPVTseedlings, mowing, TRUE, "Mowing patterns at CPVT")
map(data_ILMseedlings, mowing, TRUE, "Mowing patterns at ILM")
map(data_CPVTseedlings, rabbit_browse, TRUE, "Rabbit browse patterns at CPVT")
map(data_ILMseedlings, rabbit_browse, TRUE, "Rabbit browse patterns at ILM")
map(data_CPVTseedlings, vole_nibbling, TRUE, "Vole nibbling patterns at CPVT")
map(data_ILMseedlings, vole_nibbling, TRUE, "Vole nibbling patterns at ILM")
map(data_CPVTseedlings, deer_rub, TRUE, "Deer rub patterns at CPVT")
map(data_ILMseedlings, deer_rub, TRUE, "Deer rub patterns at ILM")
#Mowing is the only type of damage measured with enough data to see a pattern, and again, only at CPVT, where it is concentrated in the eastern half of the park.
map(data_CPVTseedlings, mowing, TRUE, "Mowing damage across ages in CPVT seedlings") + facet_wrap( ~ seedling_germ_yr)
#Mowing damage is more common in the older trees, and seems to have been largely not performed in the most recent years.  New seedlings from 2025 are almost entirely in previously mowed areas.


