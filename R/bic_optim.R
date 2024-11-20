#BIC optimization function

GetVBICOptim <- function(par, X,W,W_c,initU, option, weighted = TRUE, bic.method = 4,...)
{
  lambda <- par
  option$lambda1 <- lambda
  #fit V according to the given lambda
  #curr.v <- fit_V(X, W, initU, option, formerV = NULL)
  curr.v <- FitVWrapper(X, W,W_c, initU, option)
  V <- curr.v$V
  if(option$scale)
  {
    V <- as.matrix(curr.v$V/ curr.v$s)
  }
  fixed_first = option$fixed_ubiq
  df.dat=MatrixDFU(curr.v$V,fixed_first=fixed_first)
  bic = switch(
    bic.method,
    "1" = CalcMatrixBIC(X,W,initU,V,df=df.dat,fixed_first = fixed_first, weighted = weighted,learning = "V",...),
    '2' =  CalcMatrixBIC.loglikversion(X,W,initU,V, which.learning = "V", df = df.dat, fixed_first=fixed_first),
    '3' = NA,
    "4" = CalcMatrixBIC.loglikGLOBALversion(X,W,initU,V, W_cov = W_c, which.learning = "V", df = df.dat,fixed_first=fixed_first)
  )
  return(bic)
}


#' BIC using global fit calculation instead of sum of each regression's liklihood
#' Note that this doesn't require any fancy counting on the residuals, since all are just lumped into one regression
#' @param X original data matrix (SNPs x studies)
#' @param W uncertainty associated with X (SNPs x studies)
#' @param U SNPs x K matrix
#' @param V Studies x K matrix
#' @param W_cov Covariance matrix to correct for. For now just implemented in U.
#' @param which.learning if we are learning "V" or "U"
#' @param df degrees of freedom associated with matrix currently learning
#' @param lm.fit.residuals Specify these to scale the log-liklihood by the optimal OLS fit. Not recommended.
#'
#' @return
#' @export
CalcMatrixBIC.loglikGLOBALversion <- function(X,W,U,V, which.learning, W_cov = NULL, df = NULL, lm.fit.residuals = NULL, decomp = FALSE,...)
{
  `%>%` <- magrittr::`%>%`
  n=nrow(X)
  d = ncol(X)
  k = ncol(U)
  #Residuals are the same regardless of
  #extend_the_weighting
  #remove conditions below
  resids = calcGlobalResiduals(X,W,U,V,W_cov=W_cov,...)

  #if(which.learning == "U")
  #{
  #  resids = calcGlobalResiduals(X,W,U,V,W_cov=W_cov,...)
  #}else
  #{
  #why is this condition here? might just make thigns complex.
  #  resids = calcGlobalResiduals(X,W,U,V, W_cov = NULL,...)
  #}

  #model.ll = penLL(n*d, resids)
  model.ll = penLLSimp(n,d,resids)
  #The following was tested on 3/17. It doesn't work, my derivation is wrong. I am missing something in the log.
  #message("Compare this to...")
  #print(penLLEmpirical(n*d, resids))
  #message("Trying to look at things differently")
  #model.ll=penLLEmpirical(n*d, resids)
  if(is.null(df))
  {
    df <- MatrixDFU(V)
  }
  #message("Not scaling by n*d here anymore...")
  #ret <- -model.ll + (log(n*d)/(n*d))*df
  if(decomp)
  {
    message("log fit: ", -2*model.ll)
    message("df term: ", log(n*d)*df)
  }
  if(which.learning == "U" & df == 0)
  {
    message("Red alert: we have zeroed out everything")
    message("This means our model thinks there is no signal whatsoever. Also possible.")
  }
  ret <- -2*model.ll + log(n*d)*df
  if(!is.null(lm.fit.residuals))
  {
    lm.ll = penLL(n*d, lm.fit.residuals)
    ret <- lm.ll/(model.ll) + (log(n*d)/(n*d))*df
  }
  return(ret)
}


#' BIC using global fit calculation and the optimizing log-liklihood (a normal form)
#' Note that this doesn't require any fancy counting on the residuals, since all are just lumped into one regression
#' @param X original data matrix (SNPs x studies)
#' @param W uncertainty associated with X (SNPs x studies)
#' @param U SNPs x K matrix
#' @param V Studies x K matrix
#' @param W_cov Covariance matrix to correct for. For now just implemented in U.
#' @param which.learning if we are learning "V" or "U"
#' @param df degrees of freedom associated with matrix currently learning
#' @param lm.fit.residuals Specify these to scale the log-liklihood by the optimal OLS fit. Not recommended.
#'
#' @return
#' @export
CalcMatrixBIC.NormalLoglikGLOBALversion <- function(X,W,U,V, which.learning, W_cov = NULL, df = NULL, calc.residual.variance = FALSE, decomp = FALSE,...)
{
  `%>%` <- magrittr::`%>%`
  n=nrow(X)
  d = ncol(X)
  k = ncol(U)
  #Residuals are the same regardless of
  if(which.learning == "U")
  {
    resids = calcGlobalResiduals(X,W,U,V,W_cov=W_cov,...)
  }else
  {
    resids = calcGlobalResiduals(X,W,U,V, W_cov = NULL,...)
  }
  residual.variance =1
  if(calc.residual.variance)
  {
    residual.variance = CalcVariance(n*d, resids)
  }

  model.ll = stdLogNormalFit(resids, residual.variance=residual.variance)
  if(is.null(df))
  {
    df <- MatrixDFU(V)
  }

  if(decomp)
  {
    message("log fit: ", -2*model.ll)
    message("df term: ", log(n*d)*df)
  }
  if(which.learning == "U" & df == 0)
  {
    message("Red alert: we have zeroed out everything")
    message("This means our model thinks there is no signal whatsoever. Also possible.")
  }
  ret <- -2*model.ll + log(n*d)*df
  if(!is.null(lm.fit.residuals))
  {
    lm.ll = penLL(n*d, lm.fit.residuals)
    ret <- lm.ll/(model.ll) + (log(n*d)/(n*d))*df
  }
  return(ret)
}



GetUBICOptim <- function(par, X,W,W_c,initV, option, weighted = TRUE, bic.method = 4,...)
{
  option$alpha1 <- par
  #fit U according to the given alpha
  #curr.u <- fit_U(X, W,W_c, initV, option)
  curr.u <- FitUWrapper(X, W,W_c, initV, option)
  U <- curr.u$U
  if(option$scale)
  {
    U  <- curr.u$U /  curr.u$s
  }
  df.dat=MatrixDFU(U,fixed_first=FALSE) #don't account for with U
  bic = switch(
    bic.method,
    "1" = CalcMatrixBIC(t(X),t(W),initV,U,df=df.dat,weighted =weighted,learning = "U",...), #if you do this, more nuanced things needed.
    '2' =  CalcMatrixBIC.loglikversion(X,W,U,initV, which.learning = "U",W_cov = W_c, df = df.dat),
    '3' = NA,
    "4" = CalcMatrixBIC.loglikGLOBALversion(X,W,U,initV, which.learning = "U", W_cov = W_c, df = df.dat)
  )
  return(bic)
}


#Functions called here, not used elsewhere. To be retired:
CalcMatrixBIC <- function(X,W,U,V, W_cov = NULL, ev="std", weighted = FALSE, df = NULL, fixed_first = FALSE, learning = "V")
{
  `%>%` <- magrittr::`%>%`
  #Calculated based on my understanding of what is given by equations 12 and 7 in Lee et al 2010
  #We assume U (First term) is known/fixed, V is our predictor we evaluate here
  #11/17 UPDATE: df is calculated for U and V!
  if(!weighted)
  {
    message("unweighted, alert.")
  }
  n=nrow(X)
  d = ncol(X)
  #if(fixed_first) FF$
  #This shouldn't make a difference, resiudals will be unchanged.
  #TODO: remove this.
  if(fixed_first)  #if(FALSE) FF$
  {
    #message("Removing first factor because its fixed")
    #Regress out of X the first factor effects.
    #Remove the first factor from the downstream steps
    X <- X -  (U[,1] %*% t(V[,1]))
    V <- V[,-1]
    U <- U[,-1]
  }

  #its possible the matrices are different sizes, depending on which stage of the fitting they get passed in. Correct for this:
  #aligned <- AlignFactorMatrices(X,W,U,V); U <- aligned$U; V<- aligned$V
  k = ncol(U)
  if(!weighted)
  {
    message("not weighting here..")
    W <- matrix(1, nrow = nrow(X), ncol = ncol(X))
  }
  if(ev =="std" )
  {
    errvar <- AltVar(X,U)
  }else if(ev == "ols" & weighted)
  {
    errvar <- WeightedAltVar(X,W,U, var.method = bic.var, W_cov = W_cov)
  }
  else # anuumber gets passed in.
  {
    message("in correct place...")
    #print(ev)
    errvar = ev
  }

  if(is.null(df))
  {
    message("something went wrong,..")
    df <- MatrixDFU(V)
  }

  if(is.null(W_cov))
  {
    W_cov <- diag(nrow(X))
  }

  #I was doing this wrong previously.
  #If we reotated things, we need to rotate them back....
  if(learning == "U")
  {
    u.old <- U;
    U <- V ; V <- u.old
    X <- t(X); W <- t(W);
  }
  if(learning == "V")
  {
    message("Currently no covariance on U")
    W_cov = NULL
  }
  resid.matrix <- calcGlobalResiduals(X,W,U,V, W_cov = W_cov, fixed_first = fixed_first)
  stopifnot(any(!is.null(resid.matrix)))
  ret <- norm(resid.matrix, type = "F")^2/(n*d * errvar) + (log(n*d)/(n*d))*df
  #message("Total ", ret)
  #message(". ")

  return(ret)
}


#Determine the variance using OLS. Basically, this is the "best" residual variance we can get.
#Determine the variance using OLS. Basically, this is the "best" residual variance we can get.

AltVar <- function(X,U)
{
  resids <- c()
  for(col in 1:ncol(X))
  {
    fit <- lm(X[,col]~U)
    resids <- c(resids, resid(fit))
  }
  n = length(resids)
  (var(resids)* (n - 1)) / n
}



WeightedAltVar <- function(X,W,U, var.method = "mle", fit.resids = FALSE, W_cov = NULL)
{
  if(is.null(W_cov))
  {
    W_cov <- diag(nrow(X))
    message("filling in with blank W")
    print(dim(W_cov))
  }
  n <- nrow(X) * ncol(X)
  resids <- c()
  fit.mat.dat <- NULL
  for(col in 1:ncol(X))
  {
    w <- unlist(W[,col])
    wx <- W_cov %*% (w*X[,col])
    wu <- W_cov %*% (w * U)
    fit <- lm(wx~wu + 0)
    resids <- c(resids, resid(fit))
    if(fit.resids)
    {
      fit.mat.dat <- rbind(fit.mat.dat, unlist(residuals(fit)))
    }
  }
  stopifnot(length(resids) == n)
  r = 0
  if(var.method == "unbiased")
  {
    message("unbiased")
    p <- ncol(U) * ncol(X)
    r = (1/(n-p))*sum(resids * resids)
  }else if(var.method == "map")
  {
    message("map")
    alpha <- 0.001
    beta <- 0.001
    den <- (alpha + n/2 - 1)
    num <- beta + 0.5*sum(resids * resids)
    r <- num/den
  }else if(var.method == "mle") #MLE
  {
    #message("mle")
    p <- 0
    r = sum(resids * resids)/n
  } else #var method = std
  {
    message("std method.")
    r <- var(resids)
  }

  #print("my calculated r is:")
  #print(r)
  #print("mle var")
  #print((var(resids) * (n-1) )/ n)
  #print(r)
  #This is equivilant to all ther terms being put into one matrix
  if((sum(resids * resids) == 0))
  {
    message("No residual variance with current fit.")
    print(r)
    message("Be wary here- not clear what to do. Going to just approximate to small number")
    message("This is strong evidence of overfitting- we recommend dropping current factors.")
    message("Current U dims:")
    print(dim(U))
    #this would be easy enough to do, if we had everything in "PVE" order.
    r = 1e-10
  }
  if(fit.resids)
  {
    return(list("resid.var" = r, "resids" = as.matrix(fit.mat.dat)))
  }
  return(r)
}



WeightColumnHelper <- function(w,x,U)
{
  return(list("wx"=w*x, "wu" = w*U))
}

