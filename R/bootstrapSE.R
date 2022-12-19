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
btsp <- function(data, example="epilepsy", B,seed=NULL) {

  btsp_replicate <- function(fit, data) {

    # retrieve model output
    # estimates are for:
    ##  (Intercept) -- the intercept (B0)
    ##  age -- age as a continuous predictor (B1)
    ##  expind -- a categorical variable with two levels (0 for before and 1 for after) (B2)
    ##  expind:treat -- a categorical variable with two levels
    ### 1 for observations from individuals on the drug in the after period and 0 otherwise) (B3)

    sigmasq <- fit$sigmasq
    betas <- fit$beta
    re <- fit$re

    #simulate zi values using computed estVar MLE
    dfZi <- tibble(id = 1:length(re),
                  zi =  rnorm(max(id),0,sqrt(sigmasq)))

    #add simulated random effects to epilepsy dataframe and compute mu_i
    sim_data <- left_join(data,dfZi,by = "id") %>%
      mutate(mu_i = exp(betas[1] + betas[2]*(.data$age) + betas[3]*(.data$expind) + betas[4]*(.data$expind*.data$treat)+.data$zi))

    #generate y_ij (ie seizures)
    sim_data <- sim_data %>%
      mutate(seizures=rpois(n(),.data$mu_i)) %>%
      select(-.data$mu_i,-.data$zi)

    #fit the model on simulated data and return estimates
    epilepsy_fit <- run_model(sim_data, example)
    flagWarning <- epilepsy_fit$convergence

    #use optional information to remove instances of no convergence
    #we can't trust estimates produced if the model did not converge
    if(length(flagWarning)>1){
      return(rep(NA,5))
    } else {
      return (c(unname(epilepsy_fit$beta), unname(epilepsy_fit$sigmasq)))
    }

  }

  #set seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  #run the bootstrap function to return estimates
  epilepsy_fit <- run_model(data, example)

  #supress warning about model not fitting because we use the flag to filter out bad results
  r <- suppressWarnings(replicate(B, btsp_replicate(epilepsy_fit, data)))
  finalValues <- list(interceptSE = sd(r[1,1:B],na.rm=TRUE),
             ageSE = sd(r[2,1:B],na.rm=TRUE),
             expindSE = sd(r[3,1:B],na.rm=TRUE),
             `expind:treatSE` = sd(r[4,1:B],na.rm=TRUE),
             sigmasqSE = sd(r[5,1:B],na.rm=TRUE))

  return(finalValues)
}
