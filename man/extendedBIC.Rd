% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bic_autofit_functions.R
\name{extendedBIC}
\alias{extendedBIC}
\title{Calculate the extended term which accounts for model size, based on work by Chen et al (2008)}
\usage{
extendedBIC(bic.dat, fit, scale = FALSE)
}
\arguments{
\item{bic.dat}{object to update with extension}

\item{fit}{glmnet object}
}
\value{
extended
}
\description{
Calculate the extended term which accounts for model size, based on work by Chen et al (2008)
}
