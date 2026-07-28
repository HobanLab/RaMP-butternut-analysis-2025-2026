library(ggplot2)
library(EnvStats)
library(corrplot)

#Catherine dell'Olio

#----------------------------
#Step three
#INITIAL DATA FAMILIARIZATION
#let's visualize the distribution of our response variables (metrics of health)
hist(df$adult_percent_live_canopy)#clusters at the ends of the data
hist(df$live_adult_girdle)#mostly even across the board
hist(df$live_adult_canker_1stlivebranch)#kind of linearly dropping off from 0 to 100
hist(df$live_adult_purdue_canker_rating)#this looks almost normal, with a skew towards the right
hist(df$live_adult_purdue_canopy_ranking)#this is skewed towards the left

hist(sqrt(df$adult_percent_live_canopy))#clusters at the ends of the data
hist(sqrt(df$live_adult_girdle))#increases towards higher values
hist(sqrt(df$live_adult_canker_1stlivebranch))#kind of linearly dropping off from 0 to 100
hist(sqrt(df$live_adult_purdue_canopy_ranking))#this is now normal for the most part

#now we can start visualizing the response as related to our predictors
#fixed effects: crown class
#random effects: potentially many---each species can be an effect, as well as the grand total impact coefficient

#crown class is the only factor here with discrete levels, so it needs to be called a factor
df$adult_crown_class <- factor(df$adult_crown_class)

plot(adult_percent_live_canopy ~ adult_crown_class, data = df) #overall, no real pattern, except that the suppressed trees have higher percent live canopy
ggplot(data = df, mapping = aes(adult_crown_class, adult_percent_live_canopy)) + geom_boxplot() + stat_n_text()
kruskal.test(adult_percent_live_canopy ~ adult_crown_class, data = df)
#so there is no real difference in our sample between the canopy of each crown class, according to Kruskal and Wallis
plot(live_adult_girdle ~ adult_crown_class, data = df) #medians look different, but again, this is a small sample. Suppressed trees tend to have lower canker girdling in this sample, but again, small sample size 
ggplot(data = df, mapping = aes(adult_crown_class, live_adult_girdle)) + geom_boxplot() + stat_n_text()
kruskal.test(live_adult_girdle ~ adult_crown_class, data = df)
#so there is no real difference in our sample between the girdlement of each crown class

#NOTE: Grand_Total is ONLY meaningful for impact_TOTAl, not impact_PERCENT.  Obviously, the grand_total of percentages will sum to 100%, every time.  SO if your df is currently set to a relative calculation of impact, ignore this section!
plot(adult_percent_live_canopy ~ Grand_Total,
     data = df,
     xlab = 'Summed impact coefficient',
     ylab = 'Percent live canopy')
#what does it look like?
plot(live_adult_girdle ~ Grand_Total,
     data = df,
     xlab = 'Summed impact coefficient',
     ylab = 'Girdle % of canker')
plot(live_adult_canker_1stlivebranch ~ Grand_Total,
     data = df,
     xlab = 'Summed impact coefficient',
     ylab = 'Canker to 1st live branch')
plot(live_adult_purdue_canopy_ranking ~ Grand_Total,
     data = df,
     xlab = 'Summed impact coefficient',
     ylab = 'Purdue CANOPY ranking')
plot(live_adult_purdue_canker_rating ~ Grand_Total,
     data = df,
     xlab = 'Summed impact coefficient',
     ylab = 'Purdue CANKER ranking')

#with continuous predictor variables comes testing of correlations.  code taken from https://www.geeksforgeeks.org/machine-learning/how-to-test-for-multicollinearity-in-r/
#0.7 or above would be highly correlated
cor_matrix <- cor(df[, 2:35]) #NOTE: change to 2:35 for impact_alive ()
corrplot::corrplot(cor_matrix, method = "square", tl.pos = "tl", tl.cex = 0.8, tl.col = "black")
#lilac correlated with paper birch, elm_spp correlated with white_ash, black_walnut correlated with white_ash, black_walnut correlated with elm_spp, black_ash correlated with hickory_spp, malus correlated with paper_birch.  None of these pairs will be used because they all include at least one species with insufficient data, so I think we are good here
#I am going to run a ton of plots to see which species of interest can even be used in the model because they have enough sample size.  If the data in these plots clusters at zero, that means that for most of our focal butternuts, the impact from that species was 0.  Therefore, there is not a large enough sample of butternuts that had any impact from the species at all to include it in our models.

plot(adult_percent_live_canopy ~ `alternate-leaf_dogwood`, data = df) 
plot(adult_percent_live_canopy ~ american_elm, data = df) 
plot(adult_percent_live_canopy ~ ash_spp, data = df) 
plot(adult_percent_live_canopy ~ basswood, data = df) 
plot(adult_percent_live_canopy ~ big_tooth_aspen, data = df) 
plot(adult_percent_live_canopy ~ birch_spp, data = df) 
plot(adult_percent_live_canopy ~ bitternut_hickory, data = df) 
plot(adult_percent_live_canopy ~ black_ash, data = df) 
plot(adult_percent_live_canopy ~ black_cherry, data = df) 
plot(adult_percent_live_canopy ~ black_walnut, data = df) 
plot(adult_percent_live_canopy ~ box_elder, data = df) 
plot(adult_percent_live_canopy ~ buckthorn, data = df) 
plot(adult_percent_live_canopy ~ butternut, data = df) 
plot(adult_percent_live_canopy ~ cedar_spp, data = df) 
plot(adult_percent_live_canopy ~ elderberry, data = df) 
if(any(colnames(df) %in% "elm_spp")){
  plot(adult_percent_live_canopy ~ elm_spp, data = df)
  } #this is only a dead tree, so the if statement just checks if there is anything to plot so we don't get errors. credit to https://stackoverflow.com/questions/1169248/test-if-a-vector-contains-a-given-element
plot(adult_percent_live_canopy ~ green_ash, data = df) 
plot(adult_percent_live_canopy ~ hickory_spp, data = df) 
plot(adult_percent_live_canopy ~ honeysuckle, data = df) 
plot(adult_percent_live_canopy ~ ironwood, data = df) 
plot(adult_percent_live_canopy ~ lilac, data = df) 
plot(adult_percent_live_canopy ~ malus, data = df) 
plot(adult_percent_live_canopy ~ paper_birch, data = df)
plot(adult_percent_live_canopy ~ pin_cherry, data = df) 
plot(adult_percent_live_canopy ~ red_cedar, data = df) 
plot(adult_percent_live_canopy ~ red_maple, data = df) 
plot(adult_percent_live_canopy ~ serviceberry, data = df) 
plot(adult_percent_live_canopy ~ shagbark_hickory, data = df) 
plot(adult_percent_live_canopy ~ slippery_elm, data = df) 
plot(adult_percent_live_canopy ~ sugar_maple, data = df) 
plot(adult_percent_live_canopy ~ white_ash, data = df) 
plot(adult_percent_live_canopy ~ white_cedar, data = df) 
plot(adult_percent_live_canopy ~ white_oak, data = df) 
plot(adult_percent_live_canopy ~ yellow_birch, data = df)
#tree species that have some amount of data to examine are: white_cedar,sugar_maple,red_cedar,malus,buckthorn, butternut,bitternut_hickory,basswood, ash_spp,  shagbark_hickory,`alternate-leaf_dogwood`,american_elm

#now let's look at the impact of each of these in turn
plot(adult_percent_live_canopy ~ white_cedar,
     data = df,
     xlab = 'White cedar impact coefficient',
     ylab = 'Percent live canopy')
plot(live_adult_girdle ~ white_cedar,
     data = df,
     xlab = 'White cedar impact coefficient',
     ylab = 'Girdle % of canker')

plot(adult_percent_live_canopy ~ red_cedar,
     data = df,
     xlab = 'Red cedar impact coefficient',
     ylab = 'Percent live canopy')
plot(live_adult_girdle ~ red_cedar,
     data = df,
     xlab = 'Red cedar impact coefficient',
     ylab = 'Girdle % of canker')

plot(adult_percent_live_canopy ~ basswood,
     data = df,
     xlab = 'Basswood impact coefficient',
     ylab = 'Percent live canopy')
#this actually looks like a downward trend
plot(live_adult_girdle ~ basswood,
     data = df,
     xlab = 'Basswood impact coefficient',
     ylab = 'Girdle % of canker')

plot(adult_percent_live_canopy ~ american_elm,
     data = df,
     xlab = 'Elm impact coefficient',
     ylab = 'Percent live canopy')
#nothing much
plot(live_adult_girdle ~ american_elm,
     data = df,
     xlab = 'Elm impact coefficient',
     ylab = 'Girdle % of canker')

plot(adult_percent_live_canopy ~ buckthorn,
     data = df,
     xlab = 'Buckthorn impact coefficient',
     ylab = 'Percent live canopy')
#nothing much
plot(live_adult_girdle ~ buckthorn,
     data = df,
     xlab = 'Buckthorn impact coefficient',
     ylab = 'Girdle % of canker')
#potentially an upward trend, but also looks like an outlier is doing heavy lifting

plot(adult_percent_live_canopy ~ ash_spp,
     data = df,
     xlab = 'Ash impact coefficient',
     ylab = 'Percent live canopy')
#this actually looks like a upward trend slightly
plot(live_adult_girdle ~ ash_spp,
     data = df,
     xlab = 'Ash impact coefficient',
     ylab = 'Girdle % of canker')

plot(adult_percent_live_canopy ~ sugar_maple,
     data = df)
plot(live_adult_girdle ~ sugar_maple,
     data = df)
#again, looks like an upward trend, but there are two outliers that might be weighting.  Strange because we would have expected sugar maples to create litter that decreases girdling, though they do grow thick enough that they could be influencing local humidity

plot(adult_percent_live_canopy ~ bitternut_hickory,
     data = df)
#could be a downward trend but very slight
plot(live_adult_girdle ~ bitternut_hickory,
     data = df)
#could be a downward trend but very slight

plot(adult_percent_live_canopy ~ butternut,
     data = df)
plot(live_adult_girdle ~ butternut,
     data = df)
#this is definitely only an outlier weighting it