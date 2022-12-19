
#' Conduct bootstrapping to find standard error using point estimates of variables
#'
#' @param data The epilepsy dataframe with id, age, expind, treat and seizure columns
#' @param example The name of the dataframe - in this scenario it is "epilepsy"
#' @param B Numeric value to represent the number of times boostrapping should be done
#' @param seed optional argument to set a seed before bootstrapping so that results can be replicated
#'
#' @return A list of standard errors for each variable in the dataset
#' @export
#'
#' @import dplyr
#' @import lme4
#' @import tibble
#' @rawNamespace import(stats, except=c(filter,lag))
#'
#' @examples
#' btsp(data = epilepsy,example = "epilepsy",B = 20,seed = 1)
#' btsp(data = epilepsy,example = "epilepsy",B = 20)
btsp <- function(data, example="epilepsy", B, seed=NULL) {
  btsp_replicate <- function(fit, data, n) {
    sigmasq <- fit$sigmasq
    betas <- fit$beta
    x1 <- data$age
    x2 <- data$treat
    x3 <- data$expind
    z_i <- rnorm(n, mean=0, sd=sqrt(sigmasq))
    mu_i <- exp(betas[1] + betas[2]*x1 + betas[3]*x3 + betas[4]*x2*x3 + z_i)
    y_i <- c()
    for (i in 1:n) {
      y_i[i] <- rpois(1, mu_i[i])
    }
    sim_data <- data %>% mutate(seizures=y_i)
    epilepsy_fit <- run_model(sim_data, example)
    return (c(epilepsy_fit$beta, epilepsy_fit$sigmasq))
  }
  # set seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  epilepsy_fit <- run_model(data, example)
  n <- nrow(data)
  r <- replicate(B, btsp_replicate(epilepsy_fit, data, n))
  list(interceptSE = sd(r[1,1:B]),
       ageSE = sd(r[2,1:B]),
       expindSE = sd(r[3,1:B]),
       `expind:treatSE` = sd(r[4,1:B]),
       sigmasqSE = sd(r[5,1:B]))
}



