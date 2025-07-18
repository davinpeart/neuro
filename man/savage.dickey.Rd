% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayesian_functions.R
\name{savage.dickey}
\alias{savage.dickey}
\title{Bayesian hypothesis testing using Savage-Dickey density ratio as Bayes Factor.
Prior samples and saving of parameters when fitting required. See example code}
\usage{
savage.dickey(
  brmsfit,
  ...,
  point.null = 0,
  plot = T,
  gamma = 1000,
  colour_scheme = "blue"
)
}
\arguments{
\item{brmsfit}{An object of class brmsfit.}

\item{...}{An expression indicating the posterior to be tested against the null.}

\item{point.null}{The value of the point null hypothesis.}

\item{plot}{Logical indicating whether to return a plot of the hypothesis. If F, returns a Bayes Factor in facour of the null.}

\item{gamma}{Resolution of sampling from the prior and posterior distributions for plotting.}

\item{colour_scheme}{String indicating a single value to be supplied to nature_pallete(). Optionally
a string containing two values separated by "-". See ?nature_pallete for acceptable inputs.}
}
\value{
See argument plot.
}
\description{
Bayesian hypothesis testing using Savage-Dickey density ratio as Bayes Factor.
Prior samples and saving of parameters when fitting required. See example code
}
\examples{
data <- data.frame(a = rep(c("a", "b"), each = 30), b = c(rnorm(n = 30), rnorm(n = 30, mean = 3)))
fit <- brms::brm(b ~ 1 + a, data = data, sample_prior = "yes", save_pars = save_pars(all = TRUE))
savage.dickey(brmsfit = fit, b_a)
}
