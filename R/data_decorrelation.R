######
## Functions for adjusting for covariance structure by building a decorrelation matrix (to be used as in GLS)
######

#' Force a covariance matrix to be block diagonal
#' This matches our assumption that studies share covariance with a block diagonal structure (i.e. limited numbers are grouped)
#'
#' @param covar The actual covariance matrix; a matrix object
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
blockifyCovarianceMatrix <- function(blocks, covar)
{
  rownames(covar) <- colnames(covar)
  blockify(covar, blocks)
}

#From Guanghao fast asset
create_blocks <- function(cormat, cor_thr=0.2){
  # Hierarchical clustering
  corrdist <- stats::as.dist(1-abs(cormat))
  hc <- stats::hclust(corrdist)
  htree <- stats::cutree(hc, h=1-cor_thr) # Criterion: corr>0.2
  block <- as.integer(names(table(htree))[table(htree)>=2])
  block <- lapply(block, function(x) names(htree)[htree==x])

  return(block)
}

#Test each block

getPDStatus <- function(blocks, covar)
{
  unlist(lapply(blocks, function(n) matrixcalc::is.positive.definite(covar[n,n])))
}
containsBadBlocks <- function(blocks,covar)
{
  i=1
  pd.status <- getPDStatus(blocks,covar)
  if(length(pd.status) < 1) {return(FALSE)}
  for(bad.blocks in which(!pd.status))
  {
    message("ERROR: Cohort overlap data indicates the following phenotypes with highly similar overlap effects:")
   for(pheno in blocks[[bad.blocks]]){message("     -",pheno)}
    message("Covariance block is not PD and cannot be used in Cholesky decomposition")
    message("Consider removing one (or multiple) of these redundant phenotypes")
  }
  if(!all(pd.status)) {return(TRUE)}
  FALSE
}

blockify <- function(cormat, blocks)
{
  ret.mat <- diag(diag(cormat))
  colnames(ret.mat) <- colnames(cormat); rownames(ret.mat) <- rownames(cormat);
  if(length(blocks) != 0)
  {

    for(b in blocks)
    {
      if(length(b) >= 2)
      {
        for(i in b)
        {
          for(j in b)
          {
            ret.mat[i,j] <- cormat[i,j]
          }
        }
      }
    }
  }
  return(ret.mat)

}


#buildWhiteningMatrix(decorrelate)
# W_c <- buildWhiteningMatrix(C, ncol(X), blockfiy = TRUE)
buildWhiteningMatrix<- function(covar, dim, blockify = 0.2,...)
{
  if(is.null(covar) | all(covar == diag(dim))){
    #return(diag((dim))) #Just return an identity matrix, we aren't going to do anything....
    return(list("W_c" = diag((dim)),"C_block"=covar))
  }
  blocks <- create_blocks(covar,cor_thr=blockify,...)
  if(blockify > 0)
  {
    covar <- blockifyCovarianceMatrix(blocks, covar)
  }
  if(!isSymmetric(covar))
  {
    covar <- as.matrix(Matrix::forceSymmetric(covar))
  }

  if(containsBadBlocks(blocks, covar))
  {
    message("Program will now terminate")
    quit()
  }
  #chol.no.diag <- (chol(covar, pivot=TRUE))
  #make sure covar is symmetric and pd.
  #if(!lqmm::is.positive.definite(covar))
  if(!matrixcalc::is.positive.definite(covar))
  {
    message("C is not PD, terminating now...")
    quit()
    #covar <- lqmm::make.positive.definite(covar)
    covar <- as.matrix(Matrix::nearPD(covar, corr = TRUE, ensureSymmetry = TRUE)$mat)
    #TRY Matrix::nearPD
    #https://www.rdocumentation.org/packages/Matrix/versions/1.5-4/topics/nearPD
  }
  #chol.pd <- (chol(covar, pivot=TRUE))
  #chol.pd.nopivot <- (chol(covar))
#sometimes this transformation might modify the diagonal- force it to be 1
  #diag(covar) <- 1
  #chol.diag <- (chol(covar, pivot=TRUE))
  #I think this was wrong...
  #solve(t(chol(covar)))
  return(list("W_c" = solve((chol(covar))),"C_block"=covar))
}
#Non PD covariance matrix.....
#possible fixes:
#just force it to be PD and move on
#use the pivot feature to re-order everything, then it should work.


#Account for whitening, if sppecified:
#Covar- the covariance you want to account for
#If this is null, no change made
#whitener- don't calculate the whitening transofrm new, assume its been done already.
#adjustMatrixAsNeeded(diag(w) %*% X, decorrelate, whitener = precalc.U)
adjustMatrixAsNeeded <- function(X, covar, whitener = NULL)
{
  X <- as.matrix(X)
  if(!is.null(covar))
  {
    if(is.null(whitener))
    {
      decor.dat <- buildWhiteningMatrix(covar, ncol(X))
      U_t_inv <- decor.dat$W_c
    }else
    {
      U_t_inv <- whitener
    }
    #check the dimensions. Only need the double transpose if the dimensions aren't aligned.
    #In the case that they are, just leave alone
    if(nrow(X) == nrow(U_t_inv))
    {
      Xin <- U_t_inv %*% (X)
    }else #this is the joined matrix case, some transposition here to make it work.
    {
      #message("NEeding to transposes")
      #print(head(U_t_inv))
      #print(head(X))
      Xin <- t(U_t_inv %*% t(X))
      #Sanity check
      C2 <- var(X)
      U_t_inv <- solve(t(chol(C2)))
      Xin2 <- t(U_t_inv %*% t(X))
      if(any(round(var(Xin2), digits = 3) != diag(nrow(covar))))
      {
        message("Possible dangerous territory, error with whitening....")
        non.compliant <- which(round(var(Xin2), digits = 3) != 1 & round(var(Xin2), digits = 3) != 0)
        print(var(Xin2)[non.compliant])
      }
    }
    stopifnot(all(dim(Xin) == dim(X)))

    return(Xin)
  } else
  {
    if(!is.null(whitener))
    {
      message("No covariance structure passed, but whitener is specified.")
      message("NO CORRECTION BEING MADE.")
    }
    return(X)
  }
}



#' Simplifed WL shrinkage using a
#'
#' @param sample_cov_mat the input covariance matrix used, from LDSC
#' @param gamma the setting at which to run this
#'
#' @return shrunk estimate
#' @export
#'
#' @examples
linearShrinkLWSimple <- function(sample_cov_mat, gamma)
{
  stopifnot(gamma >= 0 & gamma <= 1)
  message("gamma set to ", gamma)
  # get the number of variables and observations
  p_n <- ncol(sample_cov_mat)

  # define the identity
  idn_pn <- diag(p_n)

  # estimate the mean scalar (average across diagonal, should just be 1?)
  m_n <- matrixStats::sum2(sample_cov_mat * idn_pn) / p_n

  estimate <- gamma * m_n * idn_pn +
    (1-gamma) * sample_cov_mat

  return(estimate)
}

#'Implementation of the covarinage shrinkage proposed by Schaefer and Strimmer in 2005,
#' "A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics"
#' (2005 "Statistical Applications in Genetics and Molecular Biology")
#' @param covar- covariance matrix, M x M
#' @param covar_se- SE for covariance estimates, M x M
#' @param sd.scaling- Scaling matrix based on standard deviation of Z-scores
#'
#' @return
#' @export
#'
#' @examples
strimmerCovShrinkage <- function(args, covar, covar_se, sd.scaling=1)
{
  if(all(covar == diag(nrow(covar))))
  {
    message("No covariance structure detected.")
    return(covar)
  }
  diag(covar_se) <- 0
  covar_se <- covar_se * sd.scaling
  covar_se[covar == 0] <- 0

  offdiag <- covar; diag(offdiag) <- 0
  gamma = sum(covar_se^2) / sum(offdiag^2) #squaring covar_se since we want variance, not SE
  print(gamma)
  args$WLgamma <- gamma
  message("Selected gamma is:", gamma)
  message("Norm prior to adjustment: ", norm(covar, "F"))
  linearShrinkLWSimple(covar, gamma)
}

#test
#cor.comp <- matrix(rnorm(100), ncol=2)
#Z <- matrix(rnorm(1000, sd= 2), ncol = 10)
#' Quickly estimate the correlation matrix that comes from non-genetic effects
#' Baesd on the method described here by Sebanti:https://deepblue.lib.umich.edu/bitstream/handle/2027.42/143992/sebanti_1.pdf?sequence=1
#' @param Z matrix of Z scores with dimensions n x m
#' @param z_thresh threshold of Z scores to omit, default is 2
#'
#' @return an
#' @export
#'
#' @examples
truncatedCovarMatrix <- function(Z, z_thresh=2)
{
  cor.mat <- diag(ncol(Z))
  for(i in 1:ncol(Z))
  {
    for(j in i:ncol(Z))
    {
      cor.comp <- Z[,c(i,j)]
      cor.comp[abs(cor.comp) < z_thresh] <- NA
      cor.comp <- data.frame(cor.comp) %>% drop_na()
      cor.mat[i,j] <- cor(cor.comp)[1,2]
      cor.mat[j,i] <- cor(cor.comp)[1,2]
    }
  }
  cor.mat
}
