#Catherine dell'Olio
#Mowing plots
#This script creates the plots comparing seedlings with and without mowing damage.
#This was intended to help land managers understand the impacts of mowing on butternut seedlings.  Do the benefits (increased light availability for butternut regeneration) outweigh the costs (mowing down new growth in unprotected butternut seedlings?)

#You will need the dataframe liveseedlings, which is a set of all living seedlings across two VT sites (ILM and CPVT)
write.csv(data_liveseedlings, "data_liveseedlings.csv")
setwd()
data_liveseedlings <- read.csv("data_liveseedlings.csv")

ggplot(data = data_liveseedlings, aes(x = height_ft, y= mowing)) +
  geom_boxplot(alpha = 0, aes(color = mowing)) +
  geom_jitter(height = 0, alpha = 0.25, size = 0.5, width = 0.3) +
  stat_n_text() +
  labs(x = "Mowing damage seen", y ="Height of seedling (ft)", title = "Height of mowed and unmowed seedlings across two sites")
library(EnvStats)


s1 <- seedlings2022$height_ft[data_liveseedlings$mowing == TRUE]
s2 <- seedlings2022$height_ft[data_liveseedlings$mowing == FALSE]
#ks(s1, s2, "Mowing and height of seedlings")
ks.test(s1, s2)
boxplot(data_liveseedlings, mowing, seedling_germ_yr, "title")

seedlings2020 <- data_liveseedlings %>% filter(seedling_germ_yr == 2020)
seedlings2021 <- data_liveseedlings %>% filter(seedling_germ_yr == 2021)
seedlings2022 <- data_liveseedlings %>% filter(seedling_germ_yr == 2022)
seedlings2023 <- data_liveseedlings %>% filter(seedling_germ_yr == 2023)
seedlings2024 <- data_liveseedlings %>% filter(seedling_germ_yr == 2024)
seedlings2025 <- data_liveseedlings %>% filter(seedling_germ_yr == 2025)

data_mowing <- data_liveseedlings %>% filter(seedling_germ_yr>2020)
boxplot(data_mowing, mowing, height_ft, "Height of seedlings when mowed (blue) and unmowed (red)") + facet_wrap( ~ seedling_germ_yr) + labs( x = "Mowing damage", y = "Height of seedling (ft)") + theme_classic() + stat_n_text()

boxplot(seedlings2021, mowing, height_ft, "title")
boxplot(seedlings2022, mowing, height_ft, "title")
boxplot(seedlings2023, mowing, height_ft, "title")
boxplot(seedlings2024, mowing, height_ft, "title")
boxplot(seedlings2025, mowing, height_ft, "title")