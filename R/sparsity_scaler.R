################################################################################################################################
## March 2022
## Ashton Omdahl
## Series of scripts for approximating the sparsity parameters, as well as
## estimating factor sparsity. This allows a user to specify a paramter space between 0 and 1, rather than exploring a wide range
##
################################################################################################################################

# Helpful stuff for debugging in Rstudio.

find_mode <- function(x) {
  u <- unique(x)
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}

FrobScale <- function(m)
{
  s <- Matrix::norm(m, type = "F")
  return(list("m.scaled" =m / s, "s"=s))
}

##Autofiting scripts
### MAP estimation approach
#4/4 observation- the problem here is that one might start to go up, pusing all the sparsity into the other one.
MAPfitLambda <- function(matin, dim, option)
{
  if(dim == 1){
    return(singleDimMAPFit(matin, option$max_lambda, option$lambda1))
  }else if(dim == 0){
    return(singleDimDropMAPFit(matin, option$max_lambda, option$lambda1))
  }
  else{
    print("Not yet implemented")
  }
}

MAPfitAlpha <- function(matin, dim, option)
{
  if(dim == 1){
    return(singleDimMAPFit(matin, option$max_alpha, option$alpha1))
 }else if(dim == 0){
    return(singleDimDropMAPFit(matin, option$max_alpha, option$alpha1))
  }
  else{
    print("Not yet implemented")
  }
}
  #modified to control past going down into
  #this is going to the minimum at each update
  #not sure if this is the best way to do it; an alternative would be to step along the gradient, as in
  #lambda - 0.01*(d/dlambda)- assume true lambda is unknown, so step towards it?
  #This update only makes sense to prevent major jumps in this, prevent us from rapidly getting too sparse or not...
  #I might like this other way better- it actually allows us to balance
  #matin: either F or L matrix
  # max: the reported MAx sparsity that would 0 everything out. Note that this should change from iteration to iteration
  #ensure that this gets updated, and isn't the same on each
#Ensuire k is updated

    singleDimMAPFit <- function(matin, max, curr, step_rate = 0.5)
    {
      #MLe approach directly
      K = ncol(matin)
      #If current is bigger than the new one, we step down. If map is bigger, we step up.
      map <- (K * nrow(matin)) / (sum(abs(matin)))
      return(updateSparsityEstimate(curr,map,step_rate, max))
    }

    singleDimDropMAPFit <- function(matin, max, curr, step_rate = 0.5)
    {
      #MAP approach, except we only count non-zero elements.
      #intuition here is that if we keep counting the same number of elements as the matrix shrinks, sparsity will keep going up.
      #in practice we don't see this though.
      #If current is bigger than the new one, we step down. If map is bigger, we step up.

      M_k_mod <- sum(matin != 0)
      message("Full size ", nrow(matin) * ncol(matin))
      message("New size ", M_k_mod)
      map <- M_k_mod / (sum(abs(matin)))

      return(updateSparsityEstimate(curr,map,step_rate, max))
    }

    updateSparsityEstimate <- function(curr, map, step_rate, max)
    {
      f <- curr - sign(curr - map) * step_rate * (map)

      #Gradient approach:
      #d_dx <- (K * nrow(matin))/curr - (sum(abs(matin)))
      #f <- curr - step_rate * d_dx
      if(is.na(f) | is.nan(f))
      {
        message("calibrated sparsity not areal number now (?) ")
        print(paste("Curr:", curr))
        print(paste("map": map))
        return(curr)
      }
      else if(f > max)
      {
        message("Warning calibrated sparsity exceeds estimated maximum...")
        print(paste("Max:", max))
        print(paste("Estimate:", f))
        return(f)
      } else{
        return(f)
      }
    }



cor2 <- function(x) {
1/(NROW(x)-1) * crossprod(scale(x, TRUE, TRUE))
}

##Using the top sparsity approach...

#11/10 version now

#Get the mode of the data based on the density function
DensityMode <- function(ln)
{
  #got this off some R help website or something
  density.vals <- density(ln)
  density.vals$x[which.max(density.vals$y)]
}

#Wrapper for quickly selecting coarse sparsity params from scratch, with burn in.
GetCoarseSparsityParams <- function(X,W,W_c,option, burn.in.reps = 1,...)
{
  message("In here...")
  param.space <- DefineSparsitySpaceInit(X, W,W_c,NULL, option, burn.in = burn.in.reps)
  return(SelectCoarseSparsityParams(param.space,burn.in.reps,...))
}

SelectCoarseSparsityParams <- function(param.space, num.reps,...)
{
  if(num.reps < 1){ num.reps <- 1 }
    #1/2 min, min, and mode
  lambda.grid <- SelectCoarseSparsityParamsGlobal(param.space$max_sparsity_params[[num.reps]]$lambda,...)

  #min, mode, and 3rd quantile.
  alpha.grid <- SelectCoarseSparsityParamsGlobal(param.space$max_sparsity_params[[num.reps]]$alpha,...)
  return(list("alphas"=alpha.grid, "lambdas"=lambda.grid))
  #list("alphas"=alpha.grid, "lambdas"=lambda.grid)
}
##NOTE: for these functions, we err on the side of density for the moment, because sparsity in practice is compounding.
#Problem arises here if we are running full iterations or just the adaptive one...
SelectCoarseSparsityParamsU <- function(param.space,n.points=3)
{
  las <- summary(param.space)
  #(c(las[1],DensityMode(param.space), las[2]))
  c(min(param.space) * 0.5, min(param.space), DensityMode(param.space))
}

SelectCoarseSparsityParamsGlobal <- function(param.space,n.points = 3, logs = TRUE)
{
  if(is.null(param.space))
  {
    return(NULL)
  }
  if(all(param.space == 0))
  {
    message("All parameters set to 0. Terminate soon.")
    return(rep(1e-10, n.points))
  }
  #The median- half of the phenotypes are zeroed out
  #the Mode:
  #mean- what it takes to zero out on average.
  #dm <- DensityMode(param.space)
  #mode is not high enough is some cases, want something more extreme. Just do the max.
  dm <- max(param.space)
  #maybe try someting a bit simpler.....
  if(!logs)
  {
    if(n.points == 5)
    {
      c(min(param.space) * 0.01, min(param.space) * 0.1, min(param.space), mean(min(param.space),dm), dm)
    } else if(n.points == 7)
    {
      div = abs(dm - min(param.space))/4
      c(min(param.space) * 0.01, min(param.space) * 0.1, min(param.space), min(param.space)+ div,min(param.space)+ 2*div, dm-div, dm)
      #c(min(param.space) * 0.15, min(param.space) * 0.30, min(param.space) * 0.45, min(param.space), min(param.space)+ div, dm-div, dm)
    }else if(n.points > 7)
    {
      n.divs <- n.points - 4
      div = abs(dm - min(param.space))/n.divs
      #seq(min(param.space), dm, by = div)
      #c(min(param.space) * 0.001, min(param.space) * 0.01, min(param.space)* 0.1, seq(min(param.space), dm, by = div))
      #10^(seq(log10(min(param.space) * 0.01), dm, length.out = n.points))
      (seq((min(param.space)), (dm), length.out = n.points))
    }
    else
    {
      c(min(param.space) * 0.1, min(param.space), dm)
    }
  }else
  {
    min.factor <- 0.1
    if(length(param.space) == 1)
    {
      min.factor <- 0.001
    }
    #Issue- if they are small, we are increasing, not decreasing
    #Update to log 10
    #s <- seq(log10(1e-5), log10(dm), length.out = (n.points + 1))
    s <- 10^(seq(log10(min(param.space)* min.factor), log10(dm), length.out = n.points))
    s
  }

}


#This is the top-level function for approximating sparsity. Also will give a starting point for K if its needed.
#@param X: input betas of effect sizes
#@param W: weights for the betas
#@param option: options associated with the funciton call
#@return max sparsity setting for L and for F
#@Todo: allow for flexibility on the sparsity settings.
approximateSparsity <- function(X, W, option){
  message('this method is deprecated dont use it anymore please')
Z <- as.matrix(X * W)
#print(is.na(Z))
#If k is unspecified, do that here too

print(option$K)
if(option$K == 0 | option$K == "kaiser" | option$K == "CnG" | option$K == "avg")
{
  library(nFactors)
  decomp <- svd(Z)
  if(option$K == 0 | option$K == "avg")
  {
    #logr::log_print("Approximating number of factors based on SVD PVE >  average PVE")
    message("Approximating number of factors based on SVD PVE >  average PVE")

    pve <- decomp$d^2 / sum(decomp$d^2)
    option$K <- sum(pve >= (mean(pve)))

  } else if(option$K == "CnG")
  {
    message("Using the top performing CnG estimate.")
    library(nFactors)
    option$K <- nCng(decomp$d)$nFactors
  } else if(option$K == "kaiser")
  {
    ops <- nScree(decomp$d)
    option$K <- ops$Components[4]

  }else{
    #logr::log_print("Using pre-specified number of components.")
    message("Using pre-specified number of components.")

  }
 #logr::log_print(paste0("Proceeding with ", option$K, " latent factors"))
 message(paste0("Proceeding with ", option$K, " latent factors"))


}else{
  decomp <- svd(Z,nu = option$K, nv = option$K)
}
  message("Estimating sparsity parameters with SVD")
  L.mat <- decomp$u[,1:option$K] %*% diag(decomp$d[1:option$K])
  F.mat <- t(diag(decomp$d[1:option$K]) %*% t(decomp$v))
  zcor <- cor2(Z)
  first <- svd(zcor)$u[,1]
  F.mat <- cbind(sign(first), decomp$v[,2:(option$K)])
  lsparsity <- sparsityParamsL(Z, F.mat, option)
  fsparsity <- sparsityParamsF(Z, L.mat, option)
  return(list("alpha" = min(lsparsity), "lambda" = min(fsparsity), "newK" = option$K, "zcor" = zcor))
}

#Estimate the MAX sparsity paramters for the loading matrix
#@param Z: input Z scores (scaled to what will be in run)
#@param FactorM: The F matrix we regress on
#@return max sparsity settings for each row of Z
sparsityParamsL<- function(Z, FactorM, option){
  L = NULL
  tS = Sys.time()
  for(row in seq(1,nrow(Z))){
    z = Z[row, ];
    l = rowiseMaxSparsity(z, as.matrix(FactorM));
    L = rbind(L, l);
  }
  updateLog("Sparsities for L estimated.", option)
  return(L)
}

#Helper function for individual row-wise sparsity; called by sparsityParamsL
#@param z: single row ofo the Z scores
#@param FactorM: The F matrix we regress on
#@return: the recommend range for a single row
rowiseMaxSparsity <- function(z, FactorM, fixed_first = FALSE){

  if(fixed_first)
  {
    #2-24 changes
    #actually from yesterday
    #may need adust this- we aren't actually going strictly against X, we maybe need to regress out 1st col?
    #message("Correcting for fixed first factor effects!")
    #fit <- lm(z ~ FactorM[,1] + 0)
    FactorM <- FactorM[,-1]
    #z <- z - fitted(fit)

  }
  #try.1 <- Matrix::t(FactorM) %*% z
  try.2 <- Matrix::t(Matrix::crossprod(z, FactorM))
  #t1norm <- norm(as.matrix(try.1), type = "I")
  t2norm <- norm(as.matrix(try.2), type = "I")
  #stopifnot(t1norm == t1norm)
  t2norm
}

#Estimate the MAX sparsity paramters for the Factor matrix
#@param Z: input Z scores (scaled to what will be in run)
#@param L: The Floadingmatrix we regress on
#@return max sparsity settings for each col of Z
sparsityParamsF <- function(Z, L, option){
  FactorM  = NULL;
  r.v <- c()
  lambda_list <- c()
  ## fit each factor one by one -- because of the element-wise multiplication from weights!
  for(col in 1:ncol(Z)){
    xp = Z[, col];
    f <- norm(t(L) %*% xp, type = "I")
    FactorM = rbind(FactorM, f);
  }
  updateLog("Sparsities for F estimated.", option)
  return(FactorM)
}



#Range recommender function. The workhorse that suggests the max
# Code was yanked from penalized
#@param response: the response variable
#@param penalized- which variables have weights on them
#@param unpenalized- which have no weights on them


recommendRange <- function (response, penalized, unpenalized, lambda1 = 100, lambda2 = 0, data,  startbeta, startgamma,
                         steps = 1, epsilon = 1e-10, maxiter, positive = FALSE, params= option)
{
  if(params$regression_method == "penalized")
  {
    trace <- FALSE
    standardize <- FALSE
    park <- FALSE
    steps = 100
    fusedl = FALSE
    if (missing(maxiter))
      maxiter <- if (lambda1 == 0 && lambda2 == 0 && !positive)
        25
    else Inf
    prep <- penalized:::.checkinput(match.call(), parent.frame())
    if (ncol(prep$X) >= nrow(prep$X) && all(lambda1 == 0) &&
        all(lambda2 == 0) && !any(prep$positive))
      stop("High-dimensional data require a penalized model. Please supply lambda1 or lambda2.",
          call. = FALSE)
    fit <- penalized:::.modelswitch(prep$model, prep$response, prep$offset,
                                    prep$strata)$fit
    pu <- length(prep$nullgamma)
    pp <- ncol(prep$X) - pu
    n <- nrow(prep$X)
    nr <- nrow(prep$X)
    fusedl <- prep$fusedl
    if (length(lambda1) == pp && (!all(lambda1 == 0))) {
      wl1 <- c(numeric(pu), lambda1)
      lambda1 <- 1
    }
    else {
      wl1 <- 1
    }

    if (park || steps > 1 && fusedl == FALSE) {
      if (pu > 0)
        lp <- drop(prep$X[, 1:pu, drop = FALSE] %*% prep$nullgamma)
      else lp <- numeric(n)
      chck <- (wl1 > 0) & c(rep(FALSE, pu), rep(TRUE, pp))
      gradient <- drop(crossprod(prep$X[, chck, drop = FALSE],
                                fit(lp)$residuals))
      if (length(wl1) > 1) {
        rel <- gradient/(wl1[chck] * prep$baselambda1[chck])
      }
      else {
        rel <- gradient/(wl1 * prep$baselambda1[chck])
      }
      from <- max(ifelse(prep$positive[chck], rel, abs(rel)))
    }
    return(from)
  }
  if(params$regression_method == "glmnet")
  {
    #TODO test and debug this
    #simply run it once with no labmda specified and get the
    penalties <- c(0, rep(1, (ncol(penalized) - 1)))
    r <- glmnet(x = penalized, y = data[,1], alpha = 1,
                    intercept = FALSE, penalty.factor = penalties)$lambda
     return(r[1])
  }
}



#Helper function for looking ath teguts of the penalized function.
#copied from their code base
penalizedDEBUG <- function (response, penalized, unpenalized, lambda1 = 0, lambda2 = 0,
                            positive = FALSE, data, fusedl = FALSE, model = c("cox",
                                                                              "logistic", "linear", "poisson"), startbeta, startgamma,
                            steps = 1, epsilon = 1e-10, maxiter, standardize = FALSE,
                            trace = FALSE)
{
  message("here")
  if (missing(maxiter))
    maxiter <- if (lambda1 == 0 && lambda2 == 0 && !positive)
      25
  else Inf
  if (steps == "Park" || steps == "park") {
    steps <- 1
    park <- TRUE
  }
  else park <- FALSE
  prep <- penalized:::.checkinput(match.call(), parent.frame())
  if (ncol(prep$X) >= nrow(prep$X) && all(lambda1 == 0) &&
      all(lambda2 == 0) && !any(prep$positive))
    stop("High-dimensional data require a penalized model. Please supply lambda1 or lambda2.",
         call. = FALSE)
  fit <- penalized:::.modelswitch(prep$model, prep$response, prep$offset,
                                  prep$strata)$fit
  pu <- length(prep$nullgamma)
  pp <- ncol(prep$X) - pu
  n <- nrow(prep$X)
  nr <- nrow(prep$X)
  fusedl <- prep$fusedl
  if (length(lambda1) == pp && (!all(lambda1 == 0))) {
    wl1 <- c(numeric(pu), lambda1)
    lambda1 <- 1
  }
  else {
    wl1 <- 1
  }
  if (length(lambda2) == pp)
    lambda2 <- c(numeric(pu), lambda2)
  print(park)
  print(steps)
  if (park || steps > 1 && fusedl == FALSE) {
    if (pu > 0)
      lp <- drop(prep$X[, 1:pu, drop = FALSE] %*% prep$nullgamma)
    else lp <- numeric(n)
    chck <- (wl1 > 0) & c(rep(FALSE, pu), rep(TRUE, pp))
    gradient <- drop(crossprod(prep$X[, chck, drop = FALSE],
                               fit(lp)$residuals))
    print(gradient)
    if (length(wl1) > 1) {
      rel <- gradient/(wl1[chck] * prep$baselambda1[chck])
    }
    else {
      rel <- gradient/(wl1 * prep$baselambda1[chck])
    }
    from <- max(ifelse(prep$positive[chck], rel, abs(rel)))
    message("Max is:")
    message(from)
    #readline()
    if (from < lambda1) {
      warning("Chosen lambda1 greater than maximal lambda1: \"steps\" argument ignored")
      steps <- 1
      park <- FALSE
      from <- lambda1
    }
  }
  else {
    message("didn't go...")
    from <- lambda1
  }
  print(lambda1)
  print(from)
  message("here?")
  lambda1s <- seq(from, lambda1, length.out = steps)
  beta <- prep$beta
  louts <- if (park)
    4 * pp
  else length(lambda1s)
  outs <- vector("list", louts)
  rellambda1 <- lambda1s[1]
  ready <- FALSE
  i <- 0
  while (!ready) {
    ready <- (rellambda1 == lambda1)
    i <- i + 1
    if (!fusedl) {
      if (rellambda1 != 0 || any(prep$positive)) {
        if (all(lambda2 == 0)) {
          out <- penalized:::.steplasso(beta = beta, lambda = rellambda1 *
                                          wl1 * prep$baselambda1, lambda2 = 0, positive = prep$positive,
                                        X = prep$X, fit = fit, trace = trace, epsilon = epsilon,
                                        maxiter = maxiter)
        }
        else {
          out <-  penalized:::.lasso(beta = beta, lambda = rellambda1 *
                                       wl1 * prep$baselambda1, lambda2 = lambda2 *
                                       prep$baselambda2, positive = prep$positive,
                                     X = prep$X, fit = fit, trace = trace, epsilon = epsilon,
                                     maxiter = maxiter)
        }
      }
      else {
        if (pp > n) {
          P <-  penalized:::.makeP(prep$X, lambda2 * prep$baselambda2)
          gams <-  penalized:::.solve(crossprod(t(P)), P %*% beta)
          PX <- P %*% t(prep$X)
          Pl <- P * matrix(sqrt(lambda2 * prep$baselambda2),
                           nrow(P), ncol(P), byrow = TRUE)
          PlP <- crossprod(t(Pl))
          out <-  penalized:::.ridge(beta = gams, Lambda = PlP, X = t(PX),
                                     fit = fit, trace = trace, epsilon = epsilon,
                                     maxiter = maxiter)
          out$beta <- drop(crossprod(P, out$beta))
        }
        else {
          out <-  penalized:::.ridge(beta = beta, Lambda = lambda2 *
                                       prep$baselambda2, X = prep$X, fit = fit,
                                     trace = trace, epsilon = epsilon, maxiter = maxiter)
        }
      }
    }
    if (fusedl) {
      out <-  penalized:::.flasso(beta = beta, lambda1 = rellambda1 *
                                    wl1 * prep$baselambda1, lambda2 = lambda2 * prep$baselambda2,
                                  chr = prep$chr, positive = prep$positive, X = prep$X,
                                  fit = fit, trace = trace, epsilon = epsilon,
                                  maxiter = maxiter)
    }
    if (trace)
      cat("\n")
    beta <- out$beta
    if (!ready) {
      if (!fusedl) {
        if (park) {
          newpark <-  penalized:::.park(beta = beta, lambda = rellambda1 *
                                          wl1 * prep$baselambda1, lambda2 = 0, positive = prep$positive,
                                        X = prep$X, fit = out$fit)
          rellambda1 <- rellambda1 * (1 - newpark$hh)
          if (rellambda1 < lambda1 || rellambda1 == Inf) {
            rellambda1 <- lambda1
            beta <- out$beta
          }
          else {
            beta <- newpark$beta
          }
          lambda1s <- c(lambda1s, rellambda1)
        }
        else {
          rellambda1 <- lambda1s[i + 1]
          beta <- out$beta
        }
      }
      else {
        rellambda1 <- lambda1s[i + 1]
        beta <- out$beta
      }
    }
    outs[[i]] <- out
  }
  if (length(lambda2) > 1)
    lambda2 <- lambda2[pu + 1:pp]
  outs <- sapply(1:i, function(nr) {
    thislambda1 <- lambda1s[[nr]] * ifelse(length(wl1) >
                                             1, wl1[pu + 1:pp], wl1)
    penalized:::.makepenfit(outs[[nr]], pu, fusedl = fusedl, prep$model,
                            thislambda1, lambda2, prep$orthogonalizer, prep$weights,
                            prep$formula, rownames(prep$X))
  })
  if (length(outs) == 1)
    outs <- outs[[1]]
  outs
}





#Copied over from Yuan:
nonEmpty <- function(v,u, iter =0 ){
	#must be empty in both u and v
  if(is.null(u) | is.null(v))
  {
    return(NA)
  }
   if(length(u) == 0 | length(v) == 0)
  {
    message("the 0 case....")
    return(0)
  }
  non_empty_v = which(apply(v, 2, function(x) sum(x!=0)) > 0)
  non_empty_u = which(apply(u, 2, function(x) sum(x!=0)) > 0)
  length(intersect(non_empty_v,non_empty_u))
}

nonEmptyAvg <- function(v,u){
  l = length(v)
  if(l != length(u))
  {
    message("Problem")
  }
  ret <- sapply(1:l, function(i) nonEmpty(v[[i]], u[[i]]))
  if(max(ret) == min(ret))
  {
    return(max(ret))
  } else{
    message("variation in repeated runs... giving an average")
    mean(ret)
  }

}



matrixSparsity <- function(m, initK, thresh = 0, wrt.init = FALSE)
{
  #I don't think this is the right way to calc it..
  #r <- sum(abs(m) <= thresh)/initK/nrow(m)
  r <- sum(abs(m) <= thresh)/(nrow(m) * ncol(m)) #This gives the local sparsity of what is there.
  if(wrt.init)
  {
    missing.cells <- (initK - ncol(m)) * nrow(m)
    r <- (sum(abs(m) <= thresh) + missing.cells) / (initK* nrow(m))
  }
  r
}


matrixSparsityAvg <- function(m,initK,...)
{
  #m is a list of matrices
  b = sapply(m, function(x) matrixSparsity(x,initK,...))
  mean(b)
}


##Helping with scaling:
getColScales <- function(matin)
{
  apply(matin,2,function(x) norm(x, "2"))
}
unitScaleColumns <- function(matin, colnorms = NA)
{
  #see https://stats.stackexchange.com/questions/8605/column-wise-matrix-normalization-in-r for speed options.
  #wordspace::normalize.cols(matin)
  if(any(is.na(colnorms)))
  {
    colnorms <- getColScales(matin)
  }

  matin %*% diag(1/colnorms)
}
