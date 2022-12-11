## code to prepare `cleanData` dataset goes here


## load dataset- no need to clean for now because
# estimation_functions.R does it for us
epilepsy = read.csv("./data-raw/epilepsy.csv")


usethis::use_data(epilepsy, overwrite = TRUE)
