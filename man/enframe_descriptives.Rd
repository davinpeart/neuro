% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayesian_functions.R
\name{enframe_descriptives}
\alias{enframe_descriptives}
\title{Calculate desriptives of a vector of observations or matrix of posterior predicive distributions}
\usage{
enframe_descriptives(y, stat1 = mean, stat2 = var)
}
\arguments{
\item{y}{Vector or matrix of values to be described.}

\item{stat1}{Function to calculate descriptive 1.}

\item{stat2}{Function to calculate descriptive 2.}
}
\value{
Dataframe containing stat1 and stat2 of y.
}
\description{
Calculate desriptives of a vector of observations or matrix of posterior predicive distributions
}
\examples{
data <- data.frame(y = rnorm(n = 100))
enframe_descriptives(y = data$y)
}
