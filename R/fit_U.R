  ################################################################################################################################
  ## Ashton Omdahl, November 2024
  ## Update the factor matrix U during the alternative optimization
  ## Input:
  ##              - X: the matrix to be decomposed (N x M, N is the number of SNPs, M is the number of traits/studies)
  ##              - V: learned factor matrix  (M x K, M is the number of traits, K is the number of factors)
  ##              - W: the weight matrix, same size as in X
  ##              - option: a list with parameters including lambda1, the l1 penalty parameter for the factor matrix (F)
  ##
  ## Return:
  ##              - A loading matrix with mixed signs that minimize_F ||W_c((X - LF') .* W)'||_F^2 + lambda1*|L|_1
  ##
  ## Example of usage:
  ##
  ## TODO- upudate
  ##
  ################################################################################################################################


#' Performs fit step for of U all SNPs simultaneously
#'
#' Internal regression step for one SNP when learning U
#' Key options here include:
#'
#' @param X unscaled coefficients for single SNP across M studies
#' @param W Standard errors associated with row of x SNPs
#' @param W_c the decorrelation covariance matrix used in the GLS form of this software.
#' @param V Fixed factor matrix from previous iteration
#' @param option List of setting options (here, function accesses: traitSpecificVar, gls, carry_coeffs,regression_method,actively_calibrating_sparsity)
#' @param formerU Matrix of U from previous iteration to use in this one
#' @param r.v List of residual variances- for use only when accounting for study-specific variance
#' @param reg.elements Pre-weighted regression elements
#'
#' @return a list containing: "l", the current learned row of U (k x 1) of weights for a given SNP, and the maximum possible sparsity for that run
#' @export
#'
#' @examples

#how can I test and compare this?
  FitUGlobal <- function(X, W,W_c, V, option, formerU, r.v = NULL,reg.elements=NULL)
  {

    max.sparsity <- NA; penalty = NA; ll = NA; l = NA; sse=NA
    #If we haven't specified the regression data we need, build it.
    if(is.null(reg.elements))
    {
      long.x <- c(t(W_c) %*% t(X*W)) #stacks by SNP1, SNP2... update so not needing to be passed each time.
      joined.weights <- lapply(1:nrow(X), function(i) Matrix::Matrix(t(W_c) %*% (diag(W[i,]))))
      if(option$std_y) { long.x <- mleStdVector(long.x)}
    }else
    {
      long.x <- reg.elements$long.x
      joined.weights <- reg.elements$joined.weights
    }
    if(option$std_y) { long.x <- mleStdVector(long.x)}

    #Extend to create the large sparse V matrix. Tested multiple approaches, this works best.
    long.v <- Matrix::bdiag(lapply(joined.weights, function(x) x %*% V))

    s = 1
    if(!is.null(formerU))
    {
      #This is an optional approach, but does little to help the objective in our tests
      formerU.scaled <- formerU + 1e-15
      old.long.u <- c(t(formerU.scaled))
    }
    if(option$regression_method == "OLS")
    {
      fit <- RcppArmadillo::fastLmPure(as.matrix(long.v), long.x)
      l = as.vector(coef(fit)) #check this, but I think its correct
      #ll = logLik(fit)
      ll = penLL(length(long.x), long.x - as.matrix(long.v) %*% fit$coefficients) # hack for now
      l = matrix(l, nrow = nrow(X), ncol = ncol(V),byrow = TRUE)
    }else if( option$regression_method == "glmnet")
    {
      if(is.na(option$alpha1)) #we are searching for parameters, phase 1 of the work.
      {
        fit <- glmnet::glmnet(x = long.v, y = long.x, family = "gaussian", alpha = 1,
                             intercept = FALSE, standardize = option$std_coef,nlambda = 100, trace.it=1) #lambda = option[['alpha1']],
        #print(paste0("Just fit the regression: ", pryr::mem_used()))
        penalty <- fit$penalty
        alpha.list <- fit$lambda

        #Get BIC scores
        bic.dat <- calculateBIC(fit, long.v, long.x, option$bic.var)
        bic.list <- bic.dat$bic.list
        min.index <- which.min(bic.list)
        bic.dat <- extractMinBicDat(bic.dat, min.index)

        u.ret = matrix(coef(fit, s = alpha.list[min.index])[-1], nrow = nrow(X), ncol = ncol(V),byrow = TRUE)

        stopifnot(!is.unsorted(-fit$lambda)) #make sure its ordered

        #ProposeNewSparsityParams(bic.list, fit$lambda, (fit$lambda), 2, n.points = 20) #need to modify this to allow for conditions.....
        #Check: is in the realm of those picked by CV?
        return(list("U" = pm(u.ret),"alpha.sel"=fit$lambda[min.index],
                    "bic"= bic.list[min.index], "sparsity_space"=max(fit$lambda), "total.log.lik" = NA,
                    "penalty" = penalty, "s"=s, "SSE"=sse, "bic.dat"=bic.dat))


      }
      else{ #just fitting a single instance:
        message("Fitting model")
        fit <- glmnet::glmnet(x = long.v, y = long.x, family = "gaussian", alpha = 1,
                             intercept = FALSE, standardize = option$std_coef,lambda = option$alpha1,trace.it=1)
        penalty <- fit$penalty
        u.ret = matrix(coef(fit, s = option$alpha1)[-1], nrow = nrow(X), ncol = ncol(V),byrow = TRUE)
        return(list("U" = u.ret, "sparsity_space"=u.ret, "total.log.lik" = NA, "penalty" = penalty, "s"=s, "SSE"=deviance(fit)))
      }
    } else if(option$regression_method == "None"){
     l = c()
    }  else
    {
      message("Nothing implemented.")
    }

      #convert back to a U:

      if(option$actively_calibrating_sparsity) { max.sparsity <- rowiseMaxSparsity(Matrix::Matrix(long.x), (long.v))}
      return(list("U" = pm(l), "sparsity_space"=max.sparsity, "total.log.lik" = ll, "penalty" = penalty, "s"=s, "SSE"=sse))
  }

#' Wrapper for the U-fitting step
#'
#' @param X the matrix to be decomposed (N x T, N is the number of data points, T is the number of features)
#' @param W the weight matrix, same size as in X, containing standard errors of the GWAS estimates
#' @param W_c the whitening matrix, associated with covariance structure, a M x M
#' @param V learned factor matrix  (M x K, M is the number of features, K is the number of factors)
#' @param option list with parameters including lambda1, the l1 penalty parameter for the U matrix
#' (relevant arguments include: ncores,)
#' @param ... arguments passed on to FitU
#'
#' @return A list containing the return matrix U, the distribution of maximum sparsity parameters, and any columns dropped because they contained NAs
#' @export
#'
  FitUWrapper <- function(X,W,W_c,V, option, prevU = NULL,...)
  {
      U = FitUGlobal(X, W,W_c, V, option,formerU = prevU, ...);
    return(U)
  }

