library(tidyverse)
library(lme4)
library(data.table)
library(logr)
source("./R/btsp.R")

ci <- function(data, example="epilepsy", B,seed=NULL,ignore.warning=TRUE) {

  # retrieve estimates from data
  estimates <- run_model(data, example)
  n <- length(estimates$re)
  btsp <- btsp(data, example, B,seed=seed,ignore.warning=ignore.warning)

  # create ci using bootstrap estimates
  interceptCI <- c(estimates$beta[1] - 1.96*btsp$interceptSE, estimates$beta[1] + 1.96*btsp$interceptSE)
  ageCI <- c(estimates$beta[2] - 1.96*btsp$ageSE, estimates$beta[2] + 1.96*btsp$ageSE)
  expindCI <- c(estimates$beta[3] - 1.96*btsp$expindSE, estimates$beta[3] + 1.96*btsp$expindSE)
  expindtreatCI <- c(estimates$beta[4] - 1.96*btsp$`expind:treatSE`, estimates$beta[4] + 1.96*btsp$`expind:treatSE`)
  sigmasqCI <- c((n-1)*estimates$sigmasq/qchisq(.025, df=n-1), (n-1)*estimates$sigmasq/qchisq(.025, df=n-1, lower.tail=F))
  list(interceptCI = unname(interceptCI),
       ageCI = unname(ageCI),
       expindCI = unname(expindCI),
       expindtreatCI = unname(expindtreatCI),
       sigmasqCI = unname(sigmasqCI))
}
