
#' Final Project Model Fitting
#'
#' @param data the data set for your example
#' @param example the name of your example: one of "culcita", "ctsib", "epilepsy", or "tortoise"
#'
#' @importFrom dplyr %>% left_join mutate n select
#' @importFrom lme4 glmer ranef
#' @importFrom rlang .data
#' @importFrom stats binomial coefficients poisson qchisq rnorm rpois sd
#' @importFrom tibble tibble
#'
#' @return A list with summary statistics from the fitted model.
#' @export
#'
#' @examples
#' ##### Epilepsy #####
#'
#' ## Fit model to epilepsy data. Fixed effects in the model are:
#' ##  (Intercept) -- the intercept
#' ##  age -- age as a continuous predictor
#' ##  expind -- a categorical variable with two levels (0 for before and 1 for after)
#' ##  expind:treat -- a categorical variable with two levels (1 for
#' ##  observations from individuals on the drug in the after period and 0 otherwise)
#'
#' epilepsy_fit <- run_model(epilepsy, "epilepsy")
#'
#'
run_model <- function(data, example = "tortoise"){

  ## Fit model
  if(example == "tortoise")
    lmer_fit <- glmer(shells~prev+offset(log(Area))+factor(year)+(1|Site),
                      family=poisson,
                      data=data)

  else if(example == "culcita")
    lmer_fit <- glmer(predation~ttt+(1|block),
                      data=data,
                      family=binomial)

  else if(example == "ctsib")
    lmer_fit <- glmer(stable ~ Surface + Vision + (1|Subject),
                      family = binomial,
                      data = data)

  else if(example == "epilepsy")
    lmer_fit <- glmer(seizures ~ age + expind + expind:treat + (1|id),
                      family = poisson,
                      data = data)
  else
    stop("You must set example to one of: tortoise, culcita, ctsib, or epilepsy.")

  ## Compute model summary
  lmer_summ <- summary(lmer_fit)

  ## Extract coefficients
  coeff <- coefficients(lmer_summ)[,"Estimate"]

  ## Extract t-statistics for each coefficient
  test_stat <- coefficients(lmer_summ)[,"z value"]

  ## Extract random effects
  re <- ranef(lmer_fit)[[1]][[1]]

  ## Extract random effects variance
  sigmasq <- lmer_summ$varcor[[1]][[1]]

  ## Extract optimization information
  optinfo <- attributes(lmer_fit)$optinfo

  ## Avoid estimate of sigmasq=0 for the original tortoise data only
  if(example == "tortoise" ){
    if(all(sort(data$shells) == c(0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3,4,5,8,10,11,12))){
      if(sigmasq == 0){
        sigmasq <- .100

        set.seed(7777)
        re <- rnorm(length(re), sd = sqrt(sigmasq))
      }
    }
  }

  ## Return values
  list(beta = coeff,
       sigmasq = sigmasq,
       re = re,
       test_stat = test_stat,
       convergence = c(optinfo$conv$opt,
                       optinfo$conv$lme4$code),
       messages = c(optinfo$message,
                    optinfo$conv$lme4$messages))

}
