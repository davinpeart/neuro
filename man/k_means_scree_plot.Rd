% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlation_functions.R
\name{k_means_scree_plot}
\alias{k_means_scree_plot}
\title{Scree plot of SS of k-means up to a specified value of k}
\usage{
k_means_scree_plot(matrix, max.k)
}
\arguments{
\item{matrix}{Numeric matrix of correlation coefficients.}

\item{max.k}{Maximum value of k up to which SS will be displayed.}
}
\value{
A ggplot scree plot of SS or k-means outputs up to max.k.
}
\description{
Scree plot of SS of k-means up to a specified value of k
}
\examples{
matrix <- cor(x = data.frame(a = rnorm(n = 20), b = rnorm(n = 20), c = rnorm(n = 20), d = rnorm(n = 20)))
k_means_scree_plot(matrix = matrix, max.k = 7)
}
