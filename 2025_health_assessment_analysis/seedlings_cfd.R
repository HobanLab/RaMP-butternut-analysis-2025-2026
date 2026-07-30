#Catherine dell'Olio, cdellolio@mortonarb.org
#This code relies on the script dataSplitter.R, which in turn relies on the data which is collected in datasheet.csv, and the package dplyr for its syntax.  The purpose of dataSplitter.R is to filter and process my datasheet in a consistent way across different analyses (with the benefit of streamlining those analysis scripts).  Please see notes in that script for the explanation of what's going on within.
#I also call the script chartBuilder.R, which is again consistent across my analyses and creates a set of handy functions so that I can standardize plots (maps, boxplots, barplots, histograms, and scatterplots) and create them in a hurry.
#The purpose of this code is to relate different metrics of seedling population structure and health to each other.

#setting working directory and running dataSplitter to get our data into readable and processable units.
setwd("~/dataAnalysisCfd")
source('dataSplitter.R')
source('chartBuilder.R')

#Let's look at the seedlings spatially and get an idea of their age and size
map(data_CPVTseedlings, seedling_germ_yr, TRUE, "CPVT seedlings by germination year")
map(data_ILMseedlings, seedling_germ_yr, TRUE, "ILM seedlings by germination year")
map(data_CPVTseedlings, height_ft, TRUE, "CPVT seedlings by height")
map(data_ILMseedlings, height_ft, TRUE, "ILM seedlings by height")
#Generally, we're seeing that the area is dominated with younger seedlings, and with shorter seedlings (makes sense).

#This graphs the seedling's germination year and its size (measured by height and basal diameter (which is a measure of stem thickness))
scatterplot(data_liveseedlings, seedling_germ_yr, height_ft, "Seedling germination year compared to height")
scatterplot(data_liveseedlings, seedling_germ_yr, basal_diameter_mm, "Seedling germination year compared to basal diameter")
#There is a clear trend downwards, as expected (younger seedlings are smaller and shorter than seedlings that have had more time to grow), but what I really find interesting is that there is much more variation within a year than between years

#These maps show the distribution of canker on live seedlings at both sites.  Again, we have multiple metrics of canker for seedlings: percent area of canker at the seedling's base, percent area along the rest of the seedling's stem, and percentage of the stem girdled (encircled) by canker).  The former two show the extent and thus severity of canker infection, whereas girdling is the canker's most direct threat to the seedling's survival.
map(data_ILMliveseedlings, seedling_canker_base, TRUE, "ILM live seedlings by canker AT BASE") + scale_color_gradient(low = "blue", high = "red")
map(data_CPVTliveseedlings, seedling_canker_base, TRUE, "CPVT live seedlings by canker AT BASE") + scale_color_gradient(low = "blue", high = "red")
map(data_ILMliveseedlings, seedling_canker_stem, TRUE, "ILM live seedlings by canker ON STEM") + scale_color_gradient(low = "blue", high = "red")
map(data_CPVTliveseedlings, seedling_canker_stem, TRUE, "CPVT live seedlings by canker ON STEM") + scale_color_gradient(low = "blue", high = "red")
map(data_ILMliveseedlings, seedling_girdle, TRUE, "ILM live seedlings by canker GIRDLING") + scale_color_gradient(low = "blue", high = "red")
map(data_CPVTliveseedlings, seedling_girdle, TRUE, "CPVT live seedlings by canker GIRDLING") + scale_color_gradient(low = "blue", high = "red")
#Seedlings are tiny, and girdling can happen really easily on them compared to mature trees---a single canker can more easily encircle 50% or more of the stem.  So we see that girdle values are higher than either base or stem, but in general, the tree are quite correlated with each other (i.e., seedlings with high canker at base tend to also have a lot of canker on the stem or girdling them, consistent with their infection.)
