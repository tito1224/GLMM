#' Epilepsy Trial Data
#'
#' Data on 59 subjects who participated in a clinical trial to test a new drug
#' called Progabide to treat epilepsy. Patients were randomly assigned into two
#' groups: treatment (31 patients) or placebo (28 patients). The number of
#' seizures each individual had was counted before the trial period and after the
#' trial period.
#'
#' @format ## `epilepsy`
#' A data frame with 118 rows and 5 columns:
#' \describe{
#'   \item{id}{ID to represent the individual}
#'   \item{age}{Contionus variable to represent age of the individual}
#'   \item{expind}{A categorical variable with two levels (0 for baseline
#'   measurement and 1 for trial measurement)}
#'   \item{treat}{A categorical variable with two levels (0 for placebo
#'   group or 1 for drug group)}
#'   ...
#' }
#' @source Leppik, I.E., et al. (1985) A double-blind crossover evaluation of
#' progabide in partial seizures. Neurology 35(285).
"epilepsy"
