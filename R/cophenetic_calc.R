
#' Calculate cophenetic correlation coefficient
#'
#' @param M a list of matrices
#'
#' @return coph a numerical correlation score in (0,1)
#' @export
#'

copheneticAssessment <- function(M){
  #compute the consensus matrix C, allowing samples to be assigned to one factor
  #Question: how does this handle cases where stuff drops out
  if(length(M) < 2)
  {
    message("Unable to calculated cophenetic correlation on just 1 sample")
    return(NA)
  }
  C = list();
  D = nrow(M[[1]])
  #I think to do this accurately you would need to scale across factors first, huh.
  for(m in 1:length(M)){
    #Scale factor to unit norm....
    scaledM <- apply(M[[m]], 2, function(x) x/norm(x, type = "2"))#
    #or should it be max scaled?
    #doesn't matter, amounts to the same thing.
    assignment = apply(scaledM, 1, function(x) which.max(abs(x))) #get the maximum column per row. Needed to adjust this for absolute value.
    ag_rep_rows = matrix(rep(assignment, D), ncol=D) #Repeat the assignment for each SNP
    C[[m]] = (ag_rep_rows == t(ag_rep_rows)) * 1; #This tells us if 2 phenotypes share the same top assignment pairwise(i.e. cell 1,2 ==1 means phenotypes 1 and 2 are most strongly loaded on factor 1) )
  }
  #take the average to get the fraction of consistent assignment
  Cbar = plyr::aaply(plyr::laply(.data=C, .fun=as.matrix), c(2, 3), mean) #Slick way of doing this- go along the entries, and
  #compute the correlation between original distances and the cophenetic distances from a hierachical clustering based on average linkage
  if(sum(Cbar != 1) == 0){
    return(1)
  }
  Ds = Cbar - diag(rep(1,D));
  if(sum(Ds) == 0){
    return(1)
  }
  d1 = dist(Ds)
  hc = hclust(d1, "ave")
  d2 = cophenetic(hc)
  coph = cor(d1, d2)

  return(coph)
}

#' Calculate the mean Frobenius norm of list of correlation matrices
#'
#' @param M list of factor matrices
#'
#' @return average of frobenius norm of correlation matrices (numeric)
#' @export
#'
correlationAssessment <- function(M)
{
  return(scaledCorrelationAssessment(M, type = "unscaled"))
}

#' Various methods to scale and aggregate the magnitude of correlation internal to several matrices
#'
#' @param M list of matrices to asses (i.e. factor matrix)
#' @param type Specify type of scaling to do: options are squared (divides by ncol^2),
#' average_r2 (takes median of average of covariance per matrix ), or unscaled (no scaling by # cols).
#' Default is just to scale by # cols in each matrix
#'
#' @return scaled/adjusted correlation structure score
#' @export
#'
scaledCorrelationAssessment <- function(M, type = "squared") {
  factor_corrs = lapply(M, cor)
  corr_frobs = unlist(lapply(factor_corrs, function(x) norm(x, 'F')))
  if(type == "squared")
  {
    mat.dims <- unlist(lapply(M, function(x) ncol(x)^2))
  }else if(type == "average_r2")
  {
    corr_avgs = unlist(lapply(factor_corrs, function(x) mean(x[upper.tri(x)]^2))) #Get the average correlation (off diagonal) per matrix
    return(median(corr_avgs))
  }else if(type == "unscaled")
  {
    mat.dims = 1
  }
  else
  {
    mat.dims <- unlist(lapply(M, function(x) ncol(x)))
  }

  return(mean(corr_frobs/mat.dims))
}
