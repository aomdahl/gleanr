% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_V.R
\name{completeUExpand}
\alias{completeUExpand}
\title{Expand out U weighted by W and W_c}
\usage{
completeUExpand(joined.weights, N, M, K, U)
}
\arguments{
\item{joined.weights}{a list of weights to scale out U by}

\item{N}{the number of SNPs}

\item{M}{the number of studies}

\item{K}{the number of factors}

\item{U}{current estimate of U}
}
\value{
a sparse matrix form of U for regression
}
\description{
Expand out U weighted by W and W_c
}
