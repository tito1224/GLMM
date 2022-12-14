library(tidyverse)
library(lme4)
library(data.table)
source("./R/estimation_functions.R")
source("./R/simulateRandomEffect.R")

bootstrapSE = function(n,df,dfname){
  lstSE = numeric(0)
  for(item in 1:n){
    # simulate zi,mu_i,y_ij and add simulated seizures to dataset
    dfSimulate = copy(epilepsy)
    dfSimulate$seizures = simulateRandomEffect(df,dfname)$seizuresSimulate

    # fit the model
    model = run_model(dfSimulate, dfname)

    # estimated beta parameters
    estBetaSimulate = model$beta
    B_expind_treat = as.numeric(unname(estBetaSimulate[4]))
    lstSE = c(lstSE,B_expind_treat)
  }

  # calculate standard error of B_expind_treat
  estSE = sd(lstSE)

  return(estSE)
}

set.seed(1)
