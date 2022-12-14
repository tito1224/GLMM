btsp <- function(data, example="epilepsy", B,seed=NULL,ignore.error=TRUE) {
  btsp_replicate <- function(fit, data) {

    # retrieve model output
    # estimates are for:
    ##  (Intercept) -- the intercept (B0)
    ##  age -- age as a continuous predictor (B1)
    ##  expind -- a categorical variable with two levels (0 for before and 1 for after) (B2)
    ##  expind:treat -- a categorical variable with two levels (1 for observations from individuals on the drug in the after period and 0 otherwise) (B3)

    sigmasq <- fit$sigmasq
    betas <- fit$beta
    re <- fit$re

    # simulate zi values using computed estVar MLE
    dfZi <- tibble(id = 1:length(re),
                  zi =  rnorm(max(id),0,sqrt(sigmasq)))

    # add simulated random effects to epilepsy dataframe and compute mu_i
    sim_data <- left_join(epilepsy,dfZi,by = "id") %>%
      mutate(mu_i = exp(betas[1] + betas[2]*age + betas[3]*expind + betas[4]*expind*treat+zi))

    # generate y_ij (ie seizures)
    simYij <- sapply(sim_data$mu_i,function(x) rpois(1,x) )
    sim_data$seizures <- simYij

    # fit the model on simulated data and return estimates
    # add some error handling to figure out which values are causing
    # lack of convergence
    if(ignore.error==TRUE){
      epilepsy_fit <- run_model(sim_data, example)
      return (c(epilepsy_fit$beta, epilepsy_fit$sigmasq))
    } else {
      tryCatch({
        epilepsy_fit <- run_model(sim_data, example)
        return (c(epilepsy_fit$beta, epilepsy_fit$sigmasq))
      }, error = function(e) {
        return(sim_data)
      })
    }
  }

  # set seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  # run the bootstrap function to return estimates
  epilepsy_fit <- run_model(data, example)
  r <- replicate(B, btsp_replicate(epilepsy_fit, data))
  list(interceptSE = sd(r[1,1:B],na.rm = TRUE),
             ageSE = sd(r[2,1:B],na.rm = TRUE),
             expindSE = sd(r[3,1:B],na.rm = TRUE),
             `expind:treatSE` = sd(r[4,1:B],na.rm = TRUE),
             sigmasqSE = sd(r[5,1:B]),na.rm = TRUE)
}

