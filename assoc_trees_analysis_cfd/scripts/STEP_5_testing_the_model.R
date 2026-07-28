library(patchwork)
library(DHARMa)

#------------------------
#Step five
#TESTING THE CHOSEN MODEL
#This will help us visualize the residuals and fitted values of our model, by storing them in a dataframe.
#i found this code online but can no longer find the link to it, so thanks to the unknown author
modelOutputPlots <- function(chosen_model){
  ModelOutputs<-data.frame(Fitted=fitted(chosen_model),
                           Residuals=resid(chosen_model))
  
  p3<-ggplot(ModelOutputs)+
    geom_point(aes(x=Fitted,y=Residuals))+
    theme_classic()+
    labs(y="Residuals",x="Fitted Values")
  
  p4<-ggplot(ModelOutputs) +
    stat_qq(aes(sample=Residuals))+
    stat_qq_line(aes(sample=Residuals))+
    theme_classic()+
    labs(y="Sample Quartiles",x="Theoretical Quartiles")
  
  #credit to https://www.datanovia.com/learn/data-visualization/ggplot2/combine-plots
  p3 + p4
  
  simulationOutput <- simulateResiduals(fittedModel = chosen_model)
  plot(simulationOutput)
  testDispersion(simulationOutput) #tests overdispersion
}

modelTest <- function(chosen_model){
  #plot(chosen_model, shade = TRUE, rug = TRUE) #this makes a few diagnostic plots that show you how well your model satisfies assumptions
  #first should be a starry sky, straight line across 0.  if not you can transform the predictors but not response with the quasibinomial model I am using
  # second should be a 1 to 1 relationship
  #third should also be nonstructured (starry sky) but okay if the red line is weird as long as there is no structure to it
  #last one tests for influential outliers, and if there is a point on or in the cook's distance lines it is a potential influential outlier
  list(summary(chosen_model),  mean(residuals(chosen_model)))
}

print("Functions modelTest(model) and modelOutputPlots(model) are ready for use!")