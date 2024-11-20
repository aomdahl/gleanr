##############################################################################################################
## Compute objective:
##############################################################################################################




getQuickObj <- function(sse,n, U,V,option)
{
  1/(2*n) * (sse) + getUPenalty(U, option) + getVPenalty(V, option)
}
#' Helper function to get the residuals across all terms
#'
#' @param X Matrix of raw effect sizes (N x M)
#' @param W Matrix of weights corresponding to effect sizes
#' @param U Current estimate of loadings (U, N x K)
#' @param V Currest estimate of factors (V, M x K)
#' @param W_c Inverse cholesky decomposition decorrelation matrix (M x M), optional
#' @param fixed_first Do we omit the first factor?
#' @param scalar is the fit scaled in any way?
#'
#' @return the matrix of residuals, N x M
#' @export
calcGlobalResiduals <- function(X,W,U,V, W_cov = NULL, fixed_first = FALSE, std.explanatory = FALSE, scalar = 1)
{

  #if(fixed_first)
  #{
  #  #pull out the effects of the 1st factor here, since we aren't interested in the fitting benefit of that
  #  #Residuals should be exactly the same at the end of the day.
  #  X <- X -  (U[,1] %*% t(V[,1]))
  #  V <- V[,-1]
  #  U <- U[,-1]
    #message("Need to think about how Covariance structure applies here...")
  #}

  if(!is.null(W_cov))
  {

    if(ncol(W_cov) != ncol(X))
    {
      message("Dimensions mismatch. Bug detected")
      return(NA)
    }
    outcome.var <-t(W_cov) %*% t(W * X)
    if(std.explanatory)
    {
      outcome.var <- mleScaleMatrix(outcome.var)
    }
    xhat <- (t(W_cov) %*% t(W * (U %*% t(V))))/scalar
    return((outcome.var -xhat) %>% t())
    #return(t(W_cov) %*% t(W * X - W * (U %*% t(V))) %>% t()) #Final transpose at end to switch the dimensions correctly. Verified this matches expected objective on 4/10
  }
  message("This case should not occur")
  explanatory.var <- (W * X)
  if(std.explanatory)
  {
    explanatory.var <- mleScaleMatrix(explanatory.var)
  }
  return((explanatory.var - (W * (U %*% t(V)))/scalar))
}

mleScaleMatrix <- function(explanatory.var)
{
  message("Updated for simplification")
  v <- scale(c(explanatory.var))
  r <- matrix(v, nrow = nrow(explanatory.var))
  #i have verified this orders the matrix correctly, rest assured.
  return(r)
  #mu = mean(explanatory.var)
  #var <- (1/(ncol(explanatory.var) * nrow(explanatory.var))) * sum((explanatory.var-mu)^2)
  #explanatory.var/sqrt(var)
}

mleScaleVector <- function(x)
{
  n <- length(x)
  mle.var <- unlist(var(x) * (n-1)/n)
  x /sqrt(mle.var)[1]
}

mleStdVector <- function(x)
{
	#simplified implementation as of Aug 31, 2024
      as.vector(scale(x))
	#x_bar<- mean(x)
 	# n <- length(x)
  #mle.x <- sum((x - x_bar)^2)/n
  #(x-x_bar) /sqrt(mle.x)[1]
}

mleStdMatrix <- function(X)
{
  n <- nrow(X)
  #pulled from https://statisticaloddsandends.wordpress.com/2018/11/15/a-deep-dive-into-glmnet-standardize/
    X_centered <- apply(X, 2, function(x) x -mean(x))
    apply(X_centered, 2, function(x) x / sqrt(sum(x^2) / n))
}

penalizedLogLik <- function(X,W,W_c, U,V,...)
{
  #simple way
  #NOTE: this will be different from the simple sum of provided terms because those calculated at a given iteration use the previous step's U to get the liklihood when estimating V.
  #This calculation gets the liklihood as is.
  penalizedLogLikV(X,W, U,V,...) + penalizedLogLikU(X,W,W_c,U,V,...)
}

stdLogNormalFit <- function(residuals, resdidual.variance = 1)
{
  -(resdidual.variance/2) * sum(residuals^2)-0.5*log(2*pi*resdidual.variance)
}

fitTermGLMNet <- function(residuals, M,N)
{
  (1/(2*M*N)) * sum(residuals^2)
}

penalizedLogLikV <- function(X,W,U,V, use.resid = NULL,...)
{
  #simple way
  total.log.lik <- 0
  message("For now, not including the covariance term in the residuals WRT V- they have no bearing on the regression step.")
  residuals = calcGlobalResiduals(X,W,U,V,...)
  if(!is.null(use.resid))
  {
    residuals = use.resid
  }
  for(j in 1:ncol(X))
  {
    n <- nrow(X)
    ss <- sum(residuals[,j]^2)
    loglik <- (-n/2) * (log(2*pi/n) + 1 + log(ss + .Machine$double.xmin))
    total.log.lik <- total.log.lik + loglik
  }
  total.log.lik
}

#' Calculate the log-liklihood of a data fit as given by penalized
#'
#' @param n the number of samples involved
#' @param resids the residuals involved, in a list
#'
#' @return the log-liklihood
#' @export
penLL <- function(n, resids)
{
  ss <- sum(resids^2)
  (-n/2) * (log(2*pi/n) + 1 + log(ss + .Machine$double.xmin))
}


#' Calculate the log-liklihood of the data fit as actually optimized by penalized
#' This is the log-normal distribution, where we assume variance =1
#' @param n the number of samples involved
#' @param resids the residuals involved, in a list
#'
#' @return the log-liklihood
#' @export
penLLEmpirical <- function(n, resids)
{
  (-1/2)*sum(resids^2) - 0.5*log(2*pi)
}


penLLSimp <- function(N, M, residuals)
{
  -M * N / 2 * log(sum(residuals^2))
}
penYuan <-function(residuals)
{
  -sum(sum(residuals^2))
}

penalizedLogLikU <- function(X,W,W_c, U,V, use.resid = NULL,...)
{
  #simple way
  total.log.lik <- 0
  fixed_first = FALSE
  residuals = calcGlobalResiduals(X,W,U,V, W_cov=W_c,...)
  if(!is.null(use.resid))
  {
    residuals = use.resid
  }
  for(i in 1:nrow(X))
  {
    n <- ncol(X)
    ss <- sum(residuals[i,]^2)
    loglik <- (-n/2) * (log(2*pi/n) + 1 + log(ss + .Machine$double.xmin))
    loglik2 = penLL(n, residuals[i,])
    stopifnot(loglik == loglik2)
    total.log.lik <- total.log.lik + loglik
  }
  total.log.lik
}

getUPenalty <- function(U, option)
{
  if(option$std_coef) {
    U <- mleStdMatrix(U)
  }
  #if(option$scale)
  #{
  #  message("not scaling U- this is the output and is directly associated with the given lambda")
  #  #U <-FrobScale(U)$m.scaled
  #}
  option[['alpha1']]* sum(abs(U)) + 0 * sqrt(sum(U^2));

}

getVPenalty <- function(V, option)
{
  if(option$std_coef) {
    V <- mleStdMatrix(V)
  }
  #if(option$scale)
  #{
  #  message("not scaling V- this is the direct output from fit and is associated with the given lambda")
 #   #V <-FrobScale(V)$m.scaled
  #}
  FactorV_penalty = option[['lambda1']] * sum(abs(V)) + 0 * sqrt(sum(V ^ 2));
  if(option$fixed_ubiq)
  {
    FactorV_penalty = 0
    if(ncol(Matrix::Matrix(V)) > 1) #only non-zero if valid scores
    {
      FactorV_penalty = option[['lambda1']] * sum(abs(V[,-1])) + 0 * sqrt(sum(V[,-1] ^ 2));

    }
  }
  return(FactorV_penalty)
}
#modified to global LL!
compute_obj <- function(X, W, W_c, L, FactorM, option, decomp = FALSE, loglik = TRUE, globalLL=TRUE, scalar = 1){

	if(is.null(FactorM) | is.null(L))
	{
		message("Unable to calculate objective for this run, because no valid matrix generated")
		return(NA)
	}
  if(is.null(option$alpha1) | is.null(option$lambda1))
  {
    message("Unable to calculate objective, no valid sparsity parameters given.")
    return(NA)
  }
  if(is.na(option$alpha1) | is.na(option$lambda1))
  {
    message("Unable to calculate objective, no valid sparsity parameters given.")
    return(NA)
  }
  M = ncol(X)
  N = nrow(X)
  #message("Scaling out L, in DEVELOPMENT")
  #L <- apply(L,2, function(x) x/norm(x, "2"))
  residuals = calcGlobalResiduals(X,W,L, FactorM, W_cov =W_c, fixed_first = FALSE, scalar =scalar, std.explanatory=option$std_y);
	#Residual_penalty = sum(sum(residual^2));
#	Residual_penalty = sum(residual^2) / (2 * nrow(L));
#After some quick evaluation, we either scale everything or don't scale at all. For simplicity we scale nothing.

	# the l2 penalty is added in the fitting step to avoid singularity, and is accounted for here
	L_penalty = getUPenalty(L, option)
	FactorM_penalty = getVPenalty(FactorM, option)
  if(!is.null(loglik))
  {
    #Residual_penalty = -loglik #passed
    #message("Replacing with my calculated log lik")
    if(globalLL)
    {
      if(option[["regression_method"]] == "glmnet" )
      {
        #negative b/c swapped below
        mine <- -fitTermGLMNet(residuals, ncol(residuals),nrow(residuals))
      }else
      {
        mine = stdLogNormalFit(residuals)
      }
    }else
    {
      mine = penalizedLogLik(X,W,W_c,L,FactorM, fixed_first = FALSE, scalar=scalar)
    }

    #message("Mine:  ", mine)
    #message("penalized_summed: ",loglik )
    Residual_penalty = -mine
    #todo= implement this manually. Can do, just wanted to "play it safe"
  } else
  {
    #derived version
    Residual_penalty = (M+N)/(2*M*N)*sum(residuals^2)
  }

  obj = Residual_penalty + FactorM_penalty + L_penalty;
  #AVerage version
  avg_v  = (1/M) * ((1/2*N) * sum(residuals)^2 + FactorM_penalty)
  avg_u = (1/N)  * ((1/2*M) * sum(residuals)^2 + L_penalty)
  #print(paste0('obj decomposition: ', (Residual_penalty), '; ', (L_penalty), '; ', (FactorM_penalty)));
  #print(paste0("Average terms: V: ", avg_v, "  U: ", avg_u))
  #reportObjectivePerFactor(X,W,L,FactorM, option)
  if(decomp)
  {
    return(list("Fit.p" = Residual_penalty, "V.p"=FactorM_penalty, "U.p"= L_penalty, "Tot" = obj))
  }
	return(obj)
}
reportObjectivePerFactor <- function(X,W,W_c,L,FactorM, option)
{
  init = 1
  if(option$fixed_ubiq)
  {
    init = 2
  }
  #Reorder factors by PVE
  #compute_obj <- function(X, W, W_c, L, FactorM, option, decomp = FALSE, loglik = TRUE, globalLL=TRUE, scalar = 1){
  #July 25, 2024- updated to a more interpretable way.
  pve.rank=order(PercentVarEx(as.matrix(X)*as.matrix(W) %*% W_c, V,U, K=ncol(V)), decreasing = TRUE)
  curr.indices <- c()
  for(i in init:ncol(FactorM))
  {
    index = which(pve.rank == i)
    curr.indices <- c(curr.indices, index)
    print(curr.indices)
    #does adding this in increase our objective
    #What do we do with the objective and ubiq 1?
    print(compute_obj(X, W,W_c, L[,curr.indices], FactorM[,curr.indices], option))
  }

}
#option <- list()
#option$alpha1 <- f.qd$autofit_alpha[1]
#option$lambda1 <- f.qd$autofit_lambda[1]
#option$fixed_ubiq <- FALSE


#  Residual_penalty = (M+N)/(2*M*N)*sum(residuals^2)
#residuals = W * (X - L %*% t(FactorM));

calcMSE <- function(X,W,Fa, L,...)
{
  message("Not accounting for covariance here, just looking at mse with data.")
  sum(calcGlobalResiduals(X,W,L,Fa)^2,...)
}

#Beta test function- can we pick out the best factors by how they contribute (or don't) to the objective?
#This will balance in sparsity better than PVE does.
DropFactorsByObjective <- function(X,W,W_c,U,V, maxK, option, minK = 0)
{

	if(is.null(maxK) | length(maxK) == 0)
	{
		maxK = ncol(X) - 1
	  message("Warning case- beware.. max k is off.")
	}
  if(maxK == 0)
  {
    message("Pruning in case where no specified limit on K.")
    maxK = ncol(V)
  }

  #1) End if only 1 column left doing fixed ubiq version
  if(ncol(as.matrix(V)) == 1 && option$fixed_ubiq)
  {
	  return(list("U"=U, "V" = V, "K" = ncol(V)))
  }
  #2) End if there are no more columns left to remove
  if(ncol(as.matrix(V)) < 1)
  {
	  message("All factors removed; none contribute to objective")
	  return(list("U"=U, "V" = V, "K" = NA))
  }

  remaining.set <- 1:ncol(V)
  #Repeat the process until conditions:
  #we reach our minimum k size threshold
  if(length(remaining.set) <=  minK)
  {
    return(list("U"=U, "V" = V, "K" = ncol(V)))
  }
  # function(X, W, L, FactorM, option)
  init.obj <- compute_obj(X,W,W_c, U, V, option)
  drop.set <- remaining.set
  if(option$fixed_ubiq)
  {
    drop.set <- 2:ncol(U)
  }
  min.obj <- init.obj
  min.index <- -1
  #Try dropping each one
  #if removing one reduces the objective, drop that one that results in biggest reduction.
  #A greedy approach to remove factors- drop the one that minimizes the objective overall
  for(i in drop.set) #we don't drop ubiq one
  {
    new.obj <- compute_obj(X,W,W_c, U[,-i], V[,-i], option, loglik = TRUE)
     if(new.obj < min.obj)
    {
      min.index <- i
      min.obj <- new.obj
    }
  }

  #Removing any objective increases (or stays the same) rather than decreases.
  if((min.obj == init.obj | min.index == -1) & (ncol(U) > maxK))
  {
    message("No factor reduces objective; none dropped.")
    return(list("U"=U, "V" = V, "K" = ncol(V)))
  }
  #3) End if dropping any more factors increases the objective and we have reached the max number of K
  else if((min.obj == init.obj | min.index == -1) & (ncol(U) <= maxK))
  {
    return(list("U"=U, "V" = V, "K" = ncol(V)))
  }else
  {
    return(DropFactorsByObjective(X,W,W_c,U[,-min.index],V[,-min.index], maxK, option, minK=minK))
  }
}

#' Calculate the change in objective with each individual matrix update
#' DEPRECATED- this doesn't make sense anymore, since we aren't updating one column or row at a time
#' @param X data
#' @param W uncertainty
#' @param W_c decorrelation matrix
#' @param old.mat the old version of the matrix you just updated
#' @param new.mat the new version of the matrix you just updated
#' @param companion.mat the complement of the matrix you just updated
#' @param fixed.term which term the complement is (u or v)
#' @param option list of options
#' @return at list of updated objectives with each change
#' @export
GetStepWiseObjective <- function(X,W,W_c,old.mat,new.mat, companion.mat,fixed.term, option)
{
  return(1)
  if(FALSE){
    iteration.objectives <- c()
    temp.mat <- old.mat
    if(ncol(old.mat) > ncol(new.mat))
    {
      for(i in 1:(ncol(old.mat) - ncol(new.mat)))
      {
        new.mat <- cbind(new.mat, rep(0, nrow(new.mat)))
      }
    }
    for(r in 1:nrow(new.mat))
    {
      temp.mat[r,] <- new.mat[r,] #update at each step.
      if(fixed.term == "U") #we were learning V
      {
        iteration.objectives <- c(iteration.objectives, compute_obj(X,W,W_c,L = companion.mat, FactorM = temp.mat, option, globalLL = TRUE))
      }else
      {
        iteration.objectives <- c(iteration.objectives, compute_obj(X,W,W_c,L = temp.mat, FactorM = companion.mat,option, globalLL = TRUE))
      }
    }
    iteration.objectives
  }

  }




#' Wrapper on objective update, for in between each iteration- adds a check for alignment that shouldn't ever get triggered.
#'
#' @param X NxM matrix of SNP effect sizes
#' @param W NxM matrix of SNP weights (1/SE)
#' @param W_c MxM matrix for decorrelating transformation
#' @param U NxK matrix of SNP loadings
#' @param V MxK matrix of trait loadings
#' @param option list of gleanr options
#' @param alpha sparsity parameter corresponding to U
#' @param lambda sparsity parameter corresponding to V
#' @param ... other arguments for compute_obj function
#'
#' @return
#' @export
#'
#' @examples
computeObjIntermediate <- function(X,W,W_c,U,V,option,alpha,lambda,...)
{
  if(ncol(U) != ncol(V))
  {
    message("Calculating objective when the matrices haven't been aligned yet. This should not happen")
    quit()
  }
  option$lambda1 <- lambda
  option$alpha1 <- alpha
  compute_obj(X, W,W_c, U, V, option, ...)
}
