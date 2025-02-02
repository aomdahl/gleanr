% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bic_autofit_functions.R
\name{trackMinParam}
\alias{trackMinParam}
\title{Updated tracking of the optimal learned parameters from gleaner}
\usage{
trackMinParam(
  min.dat,
  U,
  V,
  alpha,
  lambda,
  bic_a,
  bic_l,
  bic_a_dat,
  bic_l_dat,
  iter_tag,
  init = FALSE
)
}
\arguments{
\item{min.dat}{tracking data so far}

\item{U}{learned at settings minimizing BIC_alpha at iteration \code{iter_tag}. Should NOT have 0 columns dropped, but be in K dimension the same as V}

\item{V}{learned at settings minimizing BIC_lambda at iteration \code{iter_tag} Should NOT have 0 columns dropped, but be in K dimension the same as U}

\item{alpha}{setting minimizing BIC_alpha at iteration \code{iter_tag}}

\item{lambda}{setting minimizing BIC_lambda at iteration \code{iter_tag}}

\item{bic_a}{minimal score at iteration \code{iter_tag}}

\item{bic_l}{minimal score at iteration \code{iter_tag}}

\item{bic_a_dat}{all the data for that score, including the scaling term and the other pieces}

\item{bic_l_dat}{all the data for that score, scaling term and other pieces}

\item{iter_tag}{a string value tag to specify which iteration we are on (e.g. V1, U2, etc.)}

\item{init}{indicate if this is being called for the first time (TRUE) or not (FALSE) to initialize the object.}
}
\value{
a list containing the optimal learned parameters (U,V,alpha,lambda) and the associated BIC
}
\description{
Updated tracking of the optimal learned parameters from gleaner
}
