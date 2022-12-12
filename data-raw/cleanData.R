## code to prepare `cleanData` dataset goes here


## load dataset
epilepsy = read.csv("./data-raw/epilepsy.csv")


## clean data by adding separate observations for patients in period2
epilepsy = epilepsy %>%
  group_by(id,treat,expind,age) %>%
  summarize(seizures = sum(seizures),
            .groups = "drop")


usethis::use_data(epilepsy, overwrite = TRUE)
