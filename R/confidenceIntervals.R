ci <- function(data, example="epilepsy", B) {

  # retrieve estimates from data
  estimates <- run_model(data, example)
  n <- length(estimates$re)
  btsp <- btsp(data, example, B,1)
  interceptCI <- c(estimates$beta[1] - 1.96*btsp$interceptSE, estimates$beta[1] + 1.96*btsp$interceptSE)
  ageCI <- c(estimates$beta[2] - 1.96*btsp$ageSE, estimates$beta[2] + 1.96*btsp$ageSE)
  expindCI <- c(estimates$beta[3] - 1.96*btsp$expindSE, estimates$beta[3] + 1.96*btsp$expindSE)
  expindtreatCI <- c(estimates$beta[4] - 1.96*btsp$expindtreatSE, estimates$beta[4] + 1.96*btsp$expindtreatSE)
  sigmasqCI <- c((n-1)*estimates$sigmasq/qchisq(.025, df=n-1), (n-1)*estimates$sigmasq/qchisq(.025, df=n-1, lower.tail=F))
  list(interceptCI = interceptCI,
       ageCI = ageCI,
       expindCI = expindCI,
       expindtreatCI = expindtreatCI,
       sigmasqCI = sigmasqCI)
}
