% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pve.R
\name{linearPVE}
\alias{linearPVE}
\title{Very basic estimation of PVE, based on R2}
\usage{
linearPVE(V, U, X, W, W_c, option, K = NULL, D = NULL, jointly = TRUE)
}
\arguments{
\item{X}{adjusted for everything already}

\item{K}{number of factors}

\item{D}{for methods that have it}

\item{jointly}{if you want the PVE additively or, each factor separately (factors aren't uncorrrelated, so this may exceed actual total PVE)}
}
\value{
vector of PVE per factor.
}
\description{
Very basic estimation of PVE, based on R2
}
\examples{
#PercentVarEx(as.matrix(X)*as.matrix(W) \%*\% W_c,return.dat$V, return.dat$U)
}
