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



#' Calculate the proportion of empty cells in a matrix
#'
#' @param m matrix to estimate sparsity of
#' @param initK inital K, if the matrix is a sub-matrix (e.g. via dropped columns)
#' @param thresh threshold for counting sparsity
#' @param wrt.init if you want sparsity with respect to full inintal size of matrix (e.g. N x K)
#'
#' @return Proportion of cells with absolute values below threshold
#' @export
#'
#' @examples
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
