library(glmmTMB)
library(AICcmodavg)

#-------------------
#Step four
#DETERMINING A MODEL
#Goal here is to test out different models, and compare them with AICc (correction due to small sample size)
#It seems we have two paths according to https://stats.stackexchange.com/questions/333410/comparing-quasibinomial-glms-in-r
#*1) beta distribution using the betareg package (or glmmTMB, does the same thing with beta regression if you specify family = beta_family), and then use stats::AIC to compare
#*2) use binomial/quasibinomial and then use the quasilikelihood AIC (which is a bit more precarious)

#beta requires a predictor variable within (0, 1), without the two boundary values included.  Because we have trees with 0% or 100% canopy or girdlement, I am altering them slightly to be 1% or 99%
df_wo_boundaries <- df

for (i in 1:nrow(df_wo_boundaries)){
  if (df_wo_boundaries$adult_percent_live_canopy[i] == 0) {
    df_wo_boundaries$adult_percent_live_canopy[i] <- 0.01 
  }
  if (df_wo_boundaries$adult_percent_live_canopy[i] == 1) {
    df_wo_boundaries$adult_percent_live_canopy[i] <- 0.99
  }
}
for (i in 1:nrow(df_wo_boundaries)){
  if (df_wo_boundaries$live_adult_canker_1stlivebranch[i] == 0) {
    df_wo_boundaries$live_adult_canker_1stlivebranch[i] <- 0.01 
  }
  if (df_wo_boundaries$live_adult_canker_base[i] == 1) {
    df_wo_boundaries$live_adult_canker_base[i] <- 0.99
  }
}
for (i in 1:nrow(df_wo_boundaries)){
  if (df_wo_boundaries$live_adult_girdle[i] == 0) {
    df_wo_boundaries$live_adult_girdle[i] <- 0.01 
  }
  if (df_wo_boundaries$live_adult_girdle[i] == 1) {
    df_wo_boundaries$live_adult_girdle[i] <- 0.99
  }
}

#TESTING OUR MODELLING AND MODEL-FITTING CODE (this is good for debugging)
#we are starting off small with a simple model, linking adult_percent_live_canopy with the total impact of all associate trees (which can be thought of as the overall competition around the tree)
# glm_grand_total <-glmmTMB(df_wo_boundaries[[response_var]] ~ Grand_Total, data = df_wo_boundaries, family = beta_family(link = "logit"))
# 
# #testing our model by looking at the residuals
# ModelOutputs<-data.frame(Fitted=fitted(glm_grand_total),
#                          Residuals=resid(glm_grand_total))
# 
# p3<-ggplot(ModelOutputs)+
#   geom_point(aes(x=Fitted,y=Residuals))+
#   theme_classic()+
#   labs(y="Residuals",x="Fitted Values")
# 
# p4<-ggplot(ModelOutputs) +
#   stat_qq(aes(sample=Residuals))+
#   stat_qq_line(aes(sample=Residuals))+
#   theme_classic()+
#   labs(y="Sample Quartiles",x="Theoretical Quartiles")
# 
# p3
# p4
# summary(glm_grand_total)
# 
# #same model with different packages
# glm_grand_total_betareg <- betareg(adult_percent_live_canopy ~ Grand_Total, data = df_wo_boundaries)
# summary(glm_grand_total_betareg)
# 
# plot(residuals(glm_grand_total) ~ Grand_Total,
#      data = df) 
# abline(h = 0, lty = 2)
# 
# fittedModel  <- glm_grand_total
# simulationOutput <- simulateResiduals(fittedModel = fittedModel)
# plot(simulationOutput)
# testDispersion(simulationOutput) #not overdispersed


#testModel #if you want to run the test model, just uncomment this line!
#So this test model is bad at predicting things, because the data doesn't track well to the predictor.  It explains 4.3% of the deviance. none of the p values are significant except for the estimate of phi.
#NOTE TO SELF on the use of different packages: It's worth running the final model with both betareg and glmmTMB with family = beta_family.  They give the same coefficients, but betareg will give you the p value of phi the dispersion parameter and a pseudo R squared.  Gam puts an estimate of phi that is similar but not exact in family: beta regression(phi), and has a different deviance explained though their R squared is similar to betareg; probably not worth it.  glmmTMb plays well with DHARMA and aictab; betareg gives the same results in DHARMA with a warning message and also aictab works. Betareg and aictab give the same results if you strip the random effects structure and don't include a starting precision parameter in betareg. The jury is out if betareg can handle mixed models, so I'm choosing to continue with glmmTMB for now, but will compare with betareg.


#CHOOSING A MODEL
#if the model with only Grand_total impact is bad, what model might be better?  Here I will test a long list of models using AIC to see what combination of predictors will best model the data without overfitting.  Buckle up because it is in fact a very long list of combinations (and includes some models with interactions between certain impacts, denoted by *)
modelSelection <- function(response_var){
  modList <- list()
  
  modList[["1"]] <- glmmTMB1 <- glmmTMB(response_var ~ 1, data = df_wo_boundaries, family = beta_family)
  modList[["Grand_Total"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["basswood"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ basswood, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["white_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ white_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["red_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ red_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["malus"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ malus, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["butternut"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ butternut, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["bitternut_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ bitternut_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["buckthorn"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ buckthorn, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["ash_spp"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ ash_spp, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["shagbark_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ shagbark_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["`alternate-leaf_dogwood`"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood`, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["Grand_Total+basswood"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total+basswood, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total+white_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total+white_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total+red_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total+red_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total+malus"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total+malus, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total+butternut"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total+butternut, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total+bitternut_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total+bitternut_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total+buckthorn"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total+buckthorn, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total+ash_spp"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total+ash_spp, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total+shagbark_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total+shagbark_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total+`alternate-leaf_dogwood`"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total+`alternate-leaf_dogwood`, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total+american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total+american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["`alternate-leaf_dogwood`+basswood"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood`+basswood, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["`alternate-leaf_dogwood`+white_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood`+white_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["`alternate-leaf_dogwood`+red_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood`+red_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["`alternate-leaf_dogwood`+malus"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood`+malus, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["`alternate-leaf_dogwood`+butternut"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood`+butternut, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["`alternate-leaf_dogwood`+bitternut_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood`+bitternut_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["`alternate-leaf_dogwood`+buckthorn"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood`+buckthorn, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["`alternate-leaf_dogwood`+ash_spp"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood`+ash_spp, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["`alternate-leaf_dogwood`+shagbark_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood`+shagbark_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["`alternate-leaf_dogwood`+american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood`+american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["buckthorn+basswood"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ buckthorn+basswood, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["buckthorn+white_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ buckthorn+white_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["buckthorn+red_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ buckthorn+red_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["buckthorn+malus"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ buckthorn+malus, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["buckthorn+butternut"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ buckthorn+butternut, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["buckthorn+bitternut_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ buckthorn+bitternut_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["buckthorn+ash_spp"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ buckthorn+ash_spp, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["buckthorn+shagbark_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ buckthorn+shagbark_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["buckthorn+american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ buckthorn+american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["basswood+white_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ basswood+white_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["basswood+red_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ basswood+red_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["basswood+malus"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ basswood+malus, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["basswood+butternut"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ basswood+butternut, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["basswood+bitternut_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ basswood+bitternut_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["basswood+ash_spp"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ basswood+ash_spp, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["basswood+shagbark_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ basswood+shagbark_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["basswood+american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ basswood+american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["white_cedar+red_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ white_cedar+red_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["white_cedar+malus"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ white_cedar+malus, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["white_cedar+butternut"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ white_cedar+butternut, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["white_cedar+bitternut_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ white_cedar+bitternut_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["white_cedar+ash_spp"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ white_cedar+ash_spp, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["white_cedar+shagbark_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ white_cedar+shagbark_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["white_cedar+american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ white_cedar+american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["red_cedar+malus"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ red_cedar+malus, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["red_cedar+butternut"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ red_cedar+butternut, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["red_cedar+bitternut_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ red_cedar+bitternut_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["red_cedar+ash_spp"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ red_cedar+ash_spp, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["red_cedar+shagbark_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ red_cedar+shagbark_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["red_cedar+american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ red_cedar+american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["malus+butternut"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ malus+butternut, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["malus+bitternut_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ malus+bitternut_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["malus+ash_spp"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ malus+ash_spp, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["malus+shagbark_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ malus+shagbark_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["malus+american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ malus+american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["butternut+bitternut_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ butternut+bitternut_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["butternut+ash_spp"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ butternut+ash_spp, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["butternut+shagbark_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ butternut+shagbark_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["butternut+american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ butternut+american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["bitternut_hickory+ash_spp"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ bitternut_hickory+ash_spp, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["bitternut_hickory+shagbark_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ bitternut_hickory+shagbark_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["bitternut_hickory+american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ bitternut_hickory+american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["ash_spp+shagbark_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ ash_spp+shagbark_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["ash_spp+american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ ash_spp+american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["shagbark_hickory+american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ shagbark_hickory+american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["Grand_Total*basswood"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total*basswood, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total*white_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total*white_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total*red_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total*red_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total*malus"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total*malus, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total*butternut"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total*butternut, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total*bitternut_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total*bitternut_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total*buckthorn"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total*buckthorn, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total*ash_spp"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total*ash_spp, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total*shagbark_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total*shagbark_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total*`alternate-leaf_dogwood`"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total*`alternate-leaf_dogwood`, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total*american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total*american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["`alternate-leaf_dogwood`*basswood"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood`*basswood, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["`alternate-leaf_dogwood`*white_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood`*white_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["`alternate-leaf_dogwood`*red_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood`*red_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["`alternate-leaf_dogwood`*malus"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood`*malus, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["`alternate-leaf_dogwood`*butternut"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood`*butternut, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["`alternate-leaf_dogwood`*bitternut_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood`*bitternut_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["`alternate-leaf_dogwood`*buckthorn"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood`*buckthorn, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["`alternate-leaf_dogwood`*ash_spp"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood`*ash_spp, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["`alternate-leaf_dogwood`*shagbark_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood`*shagbark_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["`alternate-leaf_dogwood`*american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood`*american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["buckthorn*basswood"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ buckthorn*basswood, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["buckthorn*white_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ buckthorn*white_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["buckthorn*red_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ buckthorn*red_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["buckthorn*malus"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ buckthorn*malus, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["buckthorn*butternut"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ buckthorn*butternut, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["buckthorn*bitternut_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ buckthorn*bitternut_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["buckthorn*ash_spp"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ buckthorn*ash_spp, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["buckthorn*shagbark_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ buckthorn*shagbark_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["buckthorn*american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ buckthorn*american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["basswood*white_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ basswood*white_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["basswood*red_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ basswood*red_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["basswood*malus"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ basswood*malus, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["basswood*butternut"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ basswood*butternut, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["basswood*bitternut_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ basswood*bitternut_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["basswood*ash_spp"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ basswood*ash_spp, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["basswood*shagbark_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ basswood*shagbark_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["basswood*american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ basswood*american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["white_cedar*red_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ white_cedar*red_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["white_cedar*malus"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ white_cedar*malus, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["white_cedar*butternut"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ white_cedar*butternut, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["white_cedar*bitternut_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ white_cedar*bitternut_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["white_cedar*ash_spp"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ white_cedar*ash_spp, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["white_cedar*shagbark_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ white_cedar*shagbark_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["white_cedar*american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ white_cedar*american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["red_cedar*malus"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ red_cedar*malus, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["red_cedar*butternut"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ red_cedar*butternut, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["red_cedar*bitternut_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ red_cedar*bitternut_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["red_cedar*ash_spp"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ red_cedar*ash_spp, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["red_cedar*shagbark_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ red_cedar*shagbark_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["red_cedar*american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ red_cedar*american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["malus*butternut"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ malus*butternut, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["malus*bitternut_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ malus*bitternut_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["malus*ash_spp"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ malus*ash_spp, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["malus*shagbark_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ malus*shagbark_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["malus*american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ malus*american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["butternut*bitternut_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ butternut*bitternut_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["butternut*ash_spp"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ butternut*ash_spp, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["butternut*shagbark_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ butternut*shagbark_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["butternut*american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ butternut*american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["bitternut_hickory*ash_spp"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ bitternut_hickory*ash_spp, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["bitternut_hickory*shagbark_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ bitternut_hickory*shagbark_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["bitternut_hickory*american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ bitternut_hickory*american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["ash_spp*shagbark_hickory"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ ash_spp*shagbark_hickory, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["ash_spp*american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ ash_spp*american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["shagbark_hickory*american_elm"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ shagbark_hickory*american_elm, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  modList[["Grand_Total*`alternate-leaf_dogwood`*buckthorn"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total * `alternate-leaf_dogwood` * buckthorn, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["`alternate-leaf_dogwood`*malus*white_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood` * malus * white_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["`alternate-leaf_dogwood`+malus+white_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ `alternate-leaf_dogwood` + malus + white_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total*`alternate-leaf_dogwood`*malus*white_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total * `alternate-leaf_dogwood` * malus * white_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total+`alternate-leaf_dogwood`+malus+white_cedar"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total + `alternate-leaf_dogwood` + malus + white_cedar, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  modList[["Grand_Total+`alternate-leaf_dogwood`+buckthorn"]] <- glmmTMB1 <- 
    glmmTMB(response_var ~ Grand_Total + `alternate-leaf_dogwood` + buckthorn, 
            data = df_wo_boundaries, family = beta_family(link = "logit"))
  
  #I added these in because they were unexplored combinations of models, similar to some models that were performing really well, for live_adult_girdle ~.  There are too many combinations to test literally everything, so I just added some in as I went on my first go around of coding this---this is the final list.
  modList[["bitternut_hickory*red_cedar+`alternate-leaf_dogwood`"]] <- glmmTMB1 <- glmmTMB(
    response_var ~ bitternut_hickory * red_cedar + `alternate-leaf_dogwood`,
    data = df_wo_boundaries, family = beta_family
  )
  modList[["bitternut_hickory+red_cedar+`alternate-leaf_dogwood`"]] <- glmmTMB1 <-
    glmmTMB(
      response_var ~ bitternut_hickory + red_cedar + `alternate-leaf_dogwood`,
      data = df_wo_boundaries, family = beta_family(link = "logit")
    )
  modList[["bitternut_hickory+butternut+`alternate-leaf_dogwood`"]] <- glmmTMB1 <-
    glmmTMB(
      response_var ~ bitternut_hickory + butternut + `alternate-leaf_dogwood`,
      data = df_wo_boundaries, family = beta_family(link = "logit")
    )
  modList[["bitternut_hickory+butternut+american_elm+`alternate-leaf_dogwood`"]] <- glmmTMB1 <-
    glmmTMB(
      response_var ~ bitternut_hickory + butternut + american_elm + `alternate-leaf_dogwood`,
      data = df_wo_boundaries, family = beta_family(link = "logit")
    )
  modList[["bitternut_hickory+american_elm+`alternate-leaf_dogwood`"]] <- glmmTMB1 <-
    glmmTMB(
      response_var ~ bitternut_hickory + american_elm + `alternate-leaf_dogwood`,
      data = df_wo_boundaries, family = beta_family(link = "logit")
    )
  modList[["butternut+american_elm+`alternate-leaf_dogwood`"]] <- glmmTMB1 <-
    glmmTMB(
      response_var ~ butternut + american_elm + `alternate-leaf_dogwood`,
      data = df_wo_boundaries, family = beta_family(link = "logit")
    )
  return(aictab(modList)) 
}

print("function modelSelection(df_wo_boundaries$'response_var') now ready for use with the dataframe selected!")