#Catherine dell'Olio
#Mowing plots
#This script creates the plots comparing seedlings with and without mowing damage.
#This was intended to help land managers understand the impacts of mowing on butternut seedlings.  Do the benefits (increased light availability for butternut regeneration) outweigh the costs (mowing down new growth in unprotected butternut seedlings?)
library(dplyr)
library(ggplot2)
library(EnvStats)

#You will need the dataframe liveseedlings, which is a set of all living seedlings across two VT sites (ILM and CPVT)
#this should be populated from DataSplitter; alternatively, it can be uploaded from the CSV saved to this directory
#set this to be YOUR working directory path before starting!!!
working_directory <- "~/GitHub/butternut-health-assessment-2025/land_management_analysis/"
setwd(working_directory)
setwd("data")
data_liveseedlings <- read.csv("data_liveseedlings.csv")

#from my other code ChartBuilder.R, a function to make simple boxplots
boxplot <- function(data, plot_division, metric, title){
  ggplot(data, aes(x=factor({{plot_division}}), y={{metric}})) +
    geom_boxplot(alpha = 0, aes(color = {{plot_division}})) +
    geom_jitter(height = 0, alpha = 0.25, size = 0.5, width = 0.3) +
    #stat_n_text() +
    ggtitle(title)
}

#plotting the height of seedlings that have sustained mowing damage (mowing = TRUE) against those that have not sustained any mowing damage (FALSE)
ggplot(data = data_liveseedlings, aes(x = mowing, y= height_ft)) +
  geom_boxplot(alpha = 0, aes(color = mowing)) +
  geom_jitter(height = 0, alpha = 0.25, size = 0.5, width = 0.3) +
  stat_n_text() +
  labs(x = "Mowing damage seen", y ="Height of seedling (ft)", title = "Height of mowed and unmowed seedlings across two sites")

#of course, a big predictor of height for a seedling is age, so I separate the seedlings in our dataset by germination year so that I can make fair comparisons within age classes
seedlings2020 <- data_liveseedlings %>% filter(seedling_germ_yr == 2020)
seedlings2021 <- data_liveseedlings %>% filter(seedling_germ_yr == 2021)
seedlings2022 <- data_liveseedlings %>% filter(seedling_germ_yr == 2022)
seedlings2023 <- data_liveseedlings %>% filter(seedling_germ_yr == 2023)
seedlings2024 <- data_liveseedlings %>% filter(seedling_germ_yr == 2024)
seedlings2025 <- data_liveseedlings %>% filter(seedling_germ_yr == 2025)

#let's visualize the height data we have in each of these age-classes
boxplot(seedlings2021, mowing, height_ft, "Height of 2021 seedings, mowed or unmowed")
boxplot(seedlings2022, mowing, height_ft, "Height of 2022 seedings, mowed or unmowed")
boxplot(seedlings2023, mowing, height_ft, "Height of 2023 seedings, mowed or unmowed")
boxplot(seedlings2024, mowing, height_ft, "Height of 2024 seedings, mowed or unmowed")
boxplot(seedlings2025, mowing, height_ft, "Height of 2025 seedings, mowed or unmowed")

#statistically testing between these treatments (mowing damage and no mowing damage) for each individual age class using a KS test
ks.test(seedlings2020$height_ft[data_liveseedlings$mowing == TRUE], seedlings2020$height_ft[data_liveseedlings$mowing == FALSE])
ks.test(seedlings2021$height_ft[data_liveseedlings$mowing == TRUE], seedlings2021$height_ft[data_liveseedlings$mowing == FALSE])
ks.test(seedlings2022$height_ft[data_liveseedlings$mowing == TRUE], seedlings2022$height_ft[data_liveseedlings$mowing == FALSE])
ks.test(seedlings2023$height_ft[data_liveseedlings$mowing == TRUE], seedlings2023$height_ft[data_liveseedlings$mowing == FALSE])
ks.test(seedlings2024$height_ft[data_liveseedlings$mowing == TRUE], seedlings2024$height_ft[data_liveseedlings$mowing == FALSE])
ks.test(seedlings2025$height_ft[data_liveseedlings$mowing == TRUE], seedlings2025$height_ft[data_liveseedlings$mowing == FALSE])

#plotting all of these age-classes on the same plot to demonstrate the finding that there is no pattern in height measurement between mowed and unmowed seedlings
data_mowing <- data_liveseedlings %>% filter(seedling_germ_yr>2020)
boxplot(data_mowing, mowing, height_ft, "Height of seedlings when mowed (blue) and unmowed (red)") + facet_wrap( ~ seedling_germ_yr) + labs( x = "Mowing damage", y = "Height of seedling (ft)") + theme_classic() + stat_n_text()