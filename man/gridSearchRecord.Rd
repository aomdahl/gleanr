% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/K_selection_grid_search.R
\name{gridSearchRecord}
\alias{gridSearchRecord}
\title{Helper function to track the results of run from one to the next. This re-updates all the information each time, so a
bit less efficient than it could be, but needed to ensure the BICs are on the right global scale.}
\usage{
gridSearchRecord(gs_object, params, record_obj)
}
\arguments{
\item{gs_object}{return from grid search object}

\item{params}{which settings of K we tried at}

\item{record_obj}{the object this is being stored in}
}
\value{
list object containing a table of all BIC score information and the data from the test
}
\description{
Helper function to track the results of run from one to the next. This re-updates all the information each time, so a
bit less efficient than it could be, but needed to ensure the BICs are on the right global scale.
}
