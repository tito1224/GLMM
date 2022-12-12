library(tidyverse)
library(lme4)
library(data.table)
source("./R/estimation_functions.R")

# simulate our zi and mu_i values
simulateRandomEffect = function(df,dfname){
  # fit model to the epilepsy dataset
  model = run_model(df, dfname)

  # retrieve components
  # estimates are for:
  ##  (Intercept) -- the intercept (B0)
  ##  age -- age as a continuous predictor (B1)
  ##  expind -- a categorical variable with two levels (0 for before and 1 for after) (B2)
  ##  expind:treat -- a categorical variable with two levels (1 for observations from individuals on the drug in the after period and 0 otherwise) (B3)

  ## estimated beta parameters
  estBeta = model$beta

  ## estimated test statistic for beta parameters
  ts = model$test_stat

  ## estimated variance
  estVar= model$sigmasq

  ## estimated random effect (Z_i)
  estRe = model$re

  # simulate zi values using computed estVar MLE
  dfZi = tibble(id = 1:length(estRe),
                zi =  rnorm(max(id),0,sqrt(estVar)))

  # add simulated random effects to epilepsy dataframe and compute mu_i
  epilepsySimulate= left_join(epilepsy,dfZi,by = "id") %>%
    mutate(mu_i = exp(estBeta[1] + estBeta[2]*age + estBeta[3]*expind + estBeta[4]*expind*treat+zi))

  # generate y_ij (ie seizures)
  epilepsySimulate$seizuresSimulate = sapply(epilepsySimulate$mu_i,function(x) rpois(1,x) )

  # return updated dataframe with random errors + mu_i
  return(epilepsySimulate)
}

