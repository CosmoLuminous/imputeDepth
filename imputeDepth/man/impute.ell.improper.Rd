% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imputation.R
\name{impute.ell.improper}
\alias{impute.ell.improper}
\title{Improper imputation for elliptical family}
\usage{
impute.ell.improper(X, num = 5, iter.burnin = 10, iter.skip = 10,
  verbosity = 1)
}
\arguments{
\item{iter.burnin}{Number of burn-in iterations.}

\item{iter.skip}{Sampling period.}

\item{verbosity}{Verbosity parameter: \code{0} for no messages at
all, \code{1} (default) for minimal messages, \code{2} for a diagnostic mode
detailing each imputation iteration, \code{3} for the graphical mode.}

\item{m}{Number of improper-imputed data sets.}
}
\value{
A list with each element being an improper-imputed data set.
}
\description{
Depth-based improper imputation for an elliptically symmetric distribution.
}
