% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pruning_factors.R
\name{DropFactorsByFit}
\alias{DropFactorsByFit}
\title{Drop factors until a specified MaxK is reached. Options are either:
"objective", which drops the factor with the smallest impact on the overall objective (i.e. hurts the least), or
"fit", which drops the factor that hurts the fit the least.}
\usage{
DropFactorsByFit(X, W, W_c, U, V, maxK, option, calc.param = "obj", scalar = 1)
}
\arguments{
\item{X}{full data matrix}

\item{W}{standard errors}

\item{W_c}{Whitening covariance matrix}

\item{U}{Predicted U}

\item{V}{Predicted V}

\item{maxK}{K to parse down to}

\item{option}{std options}

\item{calc.param}{which parameter to use as metric, either "fit" or "obj" (which includes sparsity terms)}

\item{scalar}{for calculating objective...}
}
\value{
list containing, U, V and K
}
\description{
Drop factors until a specified MaxK is reached. Options are either:
"objective", which drops the factor with the smallest impact on the overall objective (i.e. hurts the least), or
"fit", which drops the factor that hurts the fit the least.
}
