#' Compute Confidence Interval
#'
#' @param data The epilepsy data set with id, age, expind, treat and seizure columns
#' @param example The name of the data set (called "epilepsy" in this case)
#' @param B Numeric value to represent the number of times bootstrapping should be done
#' @param seed Optional argument to set a seed before bootstrapping so that results can be saved
#'
#' @importFrom dplyr %>% left_join mutate n select
#' @importFrom lme4 glmer ranef
#' @importFrom rlang .data
#' @importFrom stats binomial coefficients poisson qchisq rnorm rpois sd
#' @importFrom tibble tibble
#'
#' @return A list of confidence intervals for each variable in the data set
#' @export
#'
#' @examples
#' ci(data = epilepsy,example = "epilepsy",B = 20,seed = 1)
#' ci(data = epilepsy,example = "epilepsy",B = 20)
ci <- function(data, example="epilepsy", B,seed=NULL) {

  #retrieve estimates from data
  estimates <- run_model(data, example)
  n <- length(estimates$re)
  resultsBtsp <- btsp(data, example, B,seed)

  #create ci using bootstrap estimates
  interceptCI <- c(estimates$beta[1] - 1.96*resultsBtsp$interceptSE, estimates$beta[1] + 1.96*resultsBtsp$interceptSE)
  ageCI <- c(estimates$beta[2] - 1.96*resultsBtsp$ageSE, estimates$beta[2] + 1.96*resultsBtsp$ageSE)
  expindCI <- c(estimates$beta[3] - 1.96*resultsBtsp$expindSE, estimates$beta[3] + 1.96*resultsBtsp$expindSE)
  expindtreatCI <- c(estimates$beta[4] - 1.96*resultsBtsp$`expind:treatSE`, estimates$beta[4] + 1.96*resultsBtsp$`expind:treatSE`)
  sigmasqCI <- c((n-1)*estimates$sigmasq/qchisq(.025, df=n-1,lower.tail = FALSE), (n-1)*estimates$sigmasq/qchisq(.975, df=n-1,lower.tail=FALSE))
  list(interceptCI = unname(interceptCI),
       ageCI = unname(ageCI),
       expindCI = unname(expindCI),
       expindtreatCI = unname(expindtreatCI),
       sigmasqCI = unname(sigmasqCI))
}
