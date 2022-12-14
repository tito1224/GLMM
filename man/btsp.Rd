% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrapSE.R
\name{btsp}
\alias{btsp}
\title{Conduct bootstrapping to find standard error using point estimates of variables}
\usage{
btsp(data, example = "epilepsy", B, seed = NULL)
}
\arguments{
\item{data}{The epilepsy dataframe with id, age, expind, treat and seizure columns}

\item{example}{The name of the dataframe - in this scenario it is "epilepsy"}

\item{B}{Numeric value to represent the number of times boostrapping should be done}

\item{seed}{optional argument to set a seed before bootstrapping so that results can be saved}
}
\value{
A list of standard errors for each variable in the dataset
}
\description{
Conduct bootstrapping to find standard error using point estimates of variables
}
\examples{
btsp(data = epilepsy,example = "epilepsy",B = 20,seed = 1)
btsp(data = epilepsy,example = "epilepsy",B = 20)
}
