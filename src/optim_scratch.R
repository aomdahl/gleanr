#BIC optimization function

GetVBICOptim <- function(par, X,W,initU, option, weighted = TRUE, bic.method = 4,...)
{
  lambda <- par
  option$lambda1 <- lambda
  #fit V according to the given lambda
  curr.v <- fit_V(X, W, initU, option, formerV = NULL)
  V <- curr.v$V
  fixed_first = option$fixed_ubiq
  df.dat=MatrixDFU(curr.v$V,fixed_first=fixed_first)
    bic = switch(
      bic.method,
      "1" = CalcMatrixBIC(X,W,initU,V,df=df.dat,fixed_first = fixed_first,...),
      '2' =  CalcMatrixBIC.loglikversion(X,W,initU,V, which.learning = "V", df = df.dat, fixed_first=fixed_first),
      '3' = NA,
      "4" = CalcMatrixBIC.loglikGLOBALversion(X,W,initU,V, which.learning = "V", df = df.dat,fixed_first=fixed_first)
    )
  return(bic)
}

GetUBICOptim <- function(par, X,W,initV, option, weighted = TRUE, bic.method = 4,...)
{
  option$alpha1 <- par
  #fit U according to the given alpha
  curr.u <- fit_U(X, W, initV, option)
  U <- curr.u$U
  fixed_first = option$fixed_ubiq
  df.dat=MatrixDFU(U,fixed_first=FALSE) #don't account for with U
  bic = switch(
    bic.method,
    "1" = CalcMatrixBIC(t(X),t(W),initV,U,df=df.dat,...), #if you do this, more nuanced things needed.
    '2' =  CalcMatrixBIC.loglikversion(X,W,U,initV, which.learning = "U", df = df.dat),
    '3' = NA,
    "4" = CalcMatrixBIC.loglikGLOBALversion(X,W,U,initV, which.learning = "U", df = df.dat)
  )
  return(bic)
}
