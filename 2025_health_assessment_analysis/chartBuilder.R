#Catherine dell'Olio
#This code creates functions to make the most common plots I use (maps, histograms, boxplots, barplots, and scatterplots) much more streamlined.  Instead of making tweaks each time in ggplots, here I have a standardized way to make plots of different butternuts and their health attributes.
#I use two packages here---ggplots2 and EnvStats (which allows me to add in stats like n= values to bar and boxplots).
#I do not call any datasheets in this code.

#Making charts
install.packages("ggplot2")
install.packages("EnvStats")
library(EnvStats)
library(ggplot2)

#spatial visualization
map <- function(data, colorz, legend_yn, title)
  ggplot(data, aes(gps_w, gps_n)) +
  geom_point(aes(color = {{ colorz }}), show.legend = legend_yn) +
  ggtitle(title)

#data spread visualization
histogram <- function(data, metric, binzwidth, title)
  ggplot(data, aes({{ metric }})) +
  geom_histogram(bins = binzwidth) +
  ggtitle(title)

histogramwidth <- function(data, metric, binzwidth, title)
  ggplot(data, aes({{ metric }})) +
  geom_histogram(binwidth = binzwidth) +
  ggtitle(title)

boxplot <- function(data, plot_division, metric, title){
  ggplot(data, aes(x=factor({{plot_division}}), y={{metric}})) +
    geom_boxplot(alpha = 0, aes(color = {{plot_division}})) +
    geom_jitter(height = 0, alpha = 0.25, size = 0.5, width = 0.3) +
    #stat_n_text() +
    ggtitle(title)
}

barplot <- function(data, variable, title){
  ggplot(data, aes({{variable}})) +
    geom_bar()+
    #stat_n_text() +
    ggtitle(title)
}

#associations
scatterplot <- function(data, x, y, title)
  ggplot(data, aes({{ x }}, {{ y }})) +
  geom_point() +
  geom_smooth(method='lm') +
  ggtitle(title)

#kolmogonov-smirnov test plots
ks <- function(s1, s2, title){
  ks_results <- ks.test(s1, s2)
  plot(ecdf(s1), xlim = range(c(s1, s2)), col = "blue", main = title)
  plot(ecdf(s2), add = TRUE, lty = "dashed", col = "red")
  title(sub= paste0("K-S Statistic = ", round(ks_results$statistic, 3)), line=2)
}
