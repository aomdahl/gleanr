################################################################################################################################
## Ashton Omdahl, November 2024
## Core functions to estimate V matrix in GLEANR's generalized alternating least squares framework
## Input:
##		- X:
##		- U:
##		- W:
##    - W_c:
##		- option: a list with parameters including lambda1, the l1 penalty parameter for the factor matrix (F)
##
## Return:
##		- A factor matrix V that minimize_F ||W_c((X - UV')* W)'||_F^2 + lambda1*|V|_1, F is non-negative
##
## Example of usage:
##
## Todo fill this out

################################################################################################################################

#' Update the factor matrix V during the alternating least squares
#'
#' @param X NxM matrix of SNP effect sizes
#' @param W NxM matrix of SNP weightings (1/SE)
#' @param W_c MxM decorrelating transformation matrix
#' @param U NxK matrix of curent estimated SNP loadings
#' @param option list of options for gleanr
#' @param formerV - option to submit previous iteration's V. Not used
#' @param reg.elements - list with pre-transformed outcome variable and weighting data. Helps with memory and speed
#'
#' @return list containing updated V and several additional key statistics, including BIC information
#' @export
#'
#' @examples
FitVGlobal <- function(X, W, W_c, U, option, formerV = NULL, reg.elements=NULL)
{
  #Prep data for regression
  K <- ncol(U)
  M = ncol(X)
  N=nrow(X)
  sparsity.est <- NA; penalty = NA; ll = NA; l = NA; FactorM = c(); penalty = NA
  if(is.null(reg.elements))
  {
   long.x <- c(t(W_c) %*% t(X*W)) #Stacked by SNP
   if(option$std_y) { long.x <- mleStdVector(long.x)}
   #Weights to multiple the reordered elements of U by
   joined.weights <- lapply(1:nrow(X), function(i) t(W_c) %*% (diag(W[i,])))
   }else   {
     long.x <- reg.elements$long.x
     joined.weights <-  reg.elements$joined.weights
   }
  s = 1
  #Expand to the full regression U.
  long.u <- completeUExpand(joined.weights, nrow(U),M,K,U)
  nopen.cols <- sapply(1:ncol(long.u), function(x) x %% K)

  if( option$regression_method == "glmnet")
  {
    lasso.active.set <- rep(1, ncol(long.u))
    #Specify the unregularized columns, if any.
    if(option$fixed_ubiq) {lasso.active.set[nopen.cols == 1] <- 0 }

    #In the model selection phase
       if(is.na(option$lambda1)) #we are still parameter searching
       {
         fit <- glmnet::glmnet(x = long.u, y = long.x, family = "gaussian", alpha = 1,
                              intercept = FALSE, standardize = option$std_coef,
                              penalty.factor = lasso.active.set,  trace.it=1) #lambda = option[['alpha1']],
         lambda.list <- fit$lambda
         penalty <- fit$penalty
        #Calculate the BIC and store the associated data
         bic.dat <- calculateBIC(fit, long.u, long.x, option$bic.var)
         bic.list <- bic.dat$bic.list
         min.index <- which.min(bic.list)
         bic.dat <- extractMinBicDat(bic.dat, min.index) #To track our progress

         v.ret = matrix(coef(fit, s = lambda.list[min.index])[-1], nrow = ncol(X), ncol = ncol(U),byrow = TRUE)
          #Make sure things happen in the right order
         stopifnot(!is.unsorted(-fit$lambda)) #make sure its ordered

         return(list("V" =  pm(v.ret),"lambda.sel"=fit$lambda[min.index],"bic"= bic.list[min.index], "sparsity_space"=max(fit$lambda),
                     "total.log.lik" = NA, "penalty" = penalty, "s"=s, "bic.dat"=bic.dat))

       }else {
         fit <- glmnet::glmnet(x = long.u, y = long.x, family = "gaussian", alpha = 1,
                              intercept = FALSE, standardize = option$std_coef, penalty.factor = lasso.active.set,
                              lambda=option$lambda1, trace.it=1)
         penalty <- fit$penalty
         v.curr = matrix(coef(fit, s = option$lambda1)[-1], nrow = ncol(X), ncol = ncol(U),byrow = TRUE)
         return(list("V" = pm(v.curr), "sparsity_space"=max(fit$lambda), "total.log.lik" = NA, "penalty" = NA, "s"=s,"SSE"=deviance(fit)))
       }

  }	else if(option$regression_method == "None"){
    FactorM = c()
  }  else   {
        message("Not implemented. Crash and burn.")
  }
  if(option$actively_calibrating_sparsity) { sparsity.est <- rowiseMaxSparsity(long.x, pen, fixed_first = FALSE)}

  return(list("V" = pm(FactorM), "sparsity_space"=sparsity.est, "total.log.lik" = ll, "penalty" = penalty, "s"=s, "bic.dat"=bic.dat))

}
#' Wrapper for the V function, not serving a purpose at the moment.
#'
#' @param X the matrix to be decomposed (N x T, N is the number of data points, T is the number of features)
#' @param W the weight matrix, same size as in X, containing standard errors of the GWAS estimates
#' @param U learned loading matrix  (N x K, N is the number of data points, K is the number of factors)
#' @param option list of all option settings (accesses V, LD, gls,
#' actively_calibrating_sparsity, posF, regression_method, epsilon,traitSpecificVar,intercept_ubiq)
#' @param formerV Optional argument to pass in previous iteration of V as prior to next one.
#'
#' @return List containing the new V ("V") and the sparsity upper limits for each regression steps ("sparsity_space")
#' @export
#'
FitVWrapper <- function(X, W,W_c, U, option, formerV = NULL,...)
{
  FitVGlobal(X, W, W_c, as.matrix(U), option, formerV = formerV,...)
}



#' Expand out U weighted by W and W_c
#'
#' @param joined.weights a list of weights to scale out U by
#' @param N the number of SNPs
#' @param M the number of studies
#' @param K the number of factors
#' @param U current estimate of U
#'
#' @return a sparse matrix form of U for regression
#' @export
#'
#' @examples
completeUExpand <- function(joined.weights, N,M,K,U)
{
  test <- do.call("rbind", lapply(1:N, function(n) matrix(c(sapply(1:M, function(m) rep(m+(n-1)*M,K)),
                                                            1:(M * K),
                                                            rep(U[n,],M)),
                                                          ncol = 3)))
  ret.dat.alt <- Matrix::sparseMatrix(i = test[,1], j = test[,2], x = test[,3], dims = c(N * M, M*K))
  Matrix::bdiag(joined.weights) %*% ret.dat.alt
}






