## Load packages
library(tidyverse)
library(lme4)

## Load helper functions
source("estimation_functions.R")

## Load data
epilepsy <- read_csv("epilepsy.csv")

## Sum separate observations for each patient in the after period
epilepsy <- epilepsy %>%
  group_by(id,treat,expind,age) %>%
  summarize(seizures = sum(seizures),
            .groups = "drop")

## Fit model to epilepsy data. Fixed effects in the model are:
##  (Intercept) -- the intercept
##  age -- age as a continuous predictor
##  expind -- a categorical variable with two levels (0 for before and 1 for after)
##  expind:treat -- a categorical variable with two levels (1 for observations from individuals on the drug in the after period and 0 otherwise)

epilepsy_fit <- run_model(epilepsy, "epilepsy")

## Components

# 1) Estimated beta parameters
epilepsy_fit$beta

# 2) Test statistics for each of the beta parameters
epilepsy_fit$test_stat

# 3) Estimated variance parameter
epilepsy_fit$sigmasq

# 4) Estimated random effects
epilepsy_fit$re

