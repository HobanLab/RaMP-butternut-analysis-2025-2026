#Catherine dell'Olio
library(spdep)
library(ggplot2)
library(patchwork)

map <- function(data, colorz, legend_yn, title)
  ggplot(data, aes(gps_w, gps_n)) +
  geom_point(aes(color = {{ colorz }}), show.legend = legend_yn) +
  ggtitle(title)

#---------------------
#Step six
#CALCULATING MORAN'S I
residuals_regression <- residuals(chosen_model, type = "pearson")

coords <- df[, c("gps_w", "gps_n")]
nb <- knn2nb(knearneigh(coords, k = 5))# Use 5 nearest neighbors
nbw <- nbw <- nb2listw(nb, style = "W")
lw <- nb2listw(nb, style = "W")


#Testing moran's i, which is our way of testing the impact of site on the residuals etc.regression <- chosen_model
morans_i_test <- function(chosen_model, df){
  moran_test <- moran.test(residuals_regression, lw)
  return(moran_test)
}

moran_plots <- function(chosen_model, response){
  #credit to https://www.paulamoraga.com/book-spatial/spatial-autocorrelation.html
  moran_testMC <- moran.mc(df$live_adult_girdle, nbw, nsim=999)
  moran_testMC
  hist(moran_testMC$res)
  abline(v = moran_testMC$statistic, col = "red")
  
  moran.plot(df$live_adult_girdle, nbw)
  
  lmoran <- localmoran(df$adult_percent_live_canopy, nbw, alternative = "greater")
  head(lmoran)
  df$lmI <- lmoran[, "Ii"] # local Moran's I
  df$lmZ <- lmoran[, "Z.Ii"] # z-scores
  # p-values corresponding to alternative greater
  df$lmp <- lmoran[, "Pr(z > E(Ii))"]
  map(df, lmI, TRUE, "Map of local moran's I")
  map(df, lmZ, TRUE, "Map of local moran's Z scores")
  map(df, lmp, TRUE, "Map of local moran's pvalues")
}

print("Function morans_i_test(chosen_model, df) and moran_plots(chosen_model, response) are now ready for use!")