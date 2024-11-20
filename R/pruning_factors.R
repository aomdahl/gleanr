#####
## Functions for pruning of factors.
######



pm <- function(matrix_in)
{
  if(is.null(matrix_in))
  {
    return(matrix_in)
  }
  m_out <-as.matrix(matrix_in)
    rownames(m_out) <- NULL
  colnames(m_out) <- NULL
  m_out
}

#' Cut out factors based on objective or some threshold (maxK)
#'
#' @param X full effect size data
#' @param W effect size uncertainty
#' @param W_c covariance structure of effect sizes
#' @param reg.run output from update_UV; full ALS output
#' @param minK Minimum allowed # of factors. NOT USED
#' @param maxK Maximum allowecx # of factors
#' @param option
#'
#' @return updated ALS output list (same as update_UV output)
#' @export
#'
PruneNumberOfFactors <- function(X,W,W_c,reg.run, minK, maxK, option)
{

  if(maxK == 0)
  {
    #this means you need to detect the number of factors)
    maxK = ncol(reg.run$V)
  }

  if(minK > maxK)
  {
    message("Invalid case encountered- minK > maxK")
    #If this is encountered, it means we've already cut out more factors than we wanted to.
    return(reg.run)
    #quit()
  }

  ret.dat <- reg.run
  r <- DropFactorsByObjective(X,W,W_c,ret.dat$U,ret.dat$V, maxK, option, minK = minK, scalar = ret.dat$final.scaling) #want it to be framed at 5, for now.
  ret.dat$V <- r$V
  ret.dat$U <- r$U
  ret.dat$K <- r$K
  if(ncol(reg.run$V) > maxK)
  {
    #print("in this case...")
    r <- DropFactorsByFit(X,W,W_c,ret.dat$U,ret.dat$V,maxK, option, scalar = ret.dat$final.scaling) #want it to be famed at 5?
    ret.dat$V <- r$V
    ret.dat$U <- r$U
    ret.dat$K <- r$K
  }
  ret.dat
}


#This will cut out factors by some specified parameter (either "fit" term or "objective") until the minimum threshold is reached.

#' Drop factors until a specified MaxK is reached. Options are either:
#'  "objective", which drops the factor with the smallest impact on the overall objective (i.e. hurts the least), or
#'  "fit", which drops the factor that hurts the fit the least.
#'
#' @param X full data matrix
#' @param W standard errors
#' @param W_c Whitening covariance matrix
#' @param U Predicted U
#' @param V Predicted V
#' @param maxK K to parse down to
#' @param option std options
#' @param calc.param which parameter to use as metric, either "fit" or "obj" (which includes sparsity terms)
#' @param scalar for calculating objective...
#'
#' @return list containing, U, V and K
#' @export
DropFactorsByFit <- function(X,W,W_c,U,V, maxK, option, calc.param="obj", scalar = 1)
{

  if(is.null(maxK) | length(maxK) == 0 | maxK == 0)
  {
    maxK = ncol(X) - 1
    message("Warning case- beware.. max k is off.")
  }

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
  if(length(remaining.set) <=  maxK)
  {
    return(list("U"=U, "V" = V, "K" = ncol(V)))
  }
  init.obj.fit <- switch(calc.param,
         "fit" = compute_obj(X,W,W_c, U, V, option, decomp = TRUE, scalar =  scalar)$Fit.p,
         "obj" = compute_obj(X,W,W_c, U, V, option, scalar =  scalar))
  drop.set <- remaining.set
  min.obj <- init.obj.fit
  min.index <- -1
  min.increase = Inf


  if(ncol(V) > maxK)
  {
    if(option$fixed_ubiq)
    {
      #Need to account for a special case here- its possible to drop the ubiquitous factor
      sub.opt <- option
      sub.opt$fixed_ubiq <- FALSE
      UFI=1 #ubiquitous factor index
      new.obj.fit <- switch(calc.param,
                            "fit" = compute_obj(X,W,W_c, U[,-UFI], V[,-UFI], sub.opt, decomp = TRUE, scalar =  scalar)$Fit.p,
                            "obj" =  compute_obj(X,W,W_c, U[,-UFI], V[,-UFI], sub.opt, scalar =  scalar))
      diff = new.obj.fit - init.obj.fit
      if(diff < min.increase)
      {
        min.index <- UFI
        min.increase <- diff
      }
      drop.set <- 2:ncol(U)
    }

    for(i in drop.set) #we don't drop ubiq one
    {
      new.obj.fit <- switch(calc.param,
                            "fit" = compute_obj(X,W,W_c, U[,-i], V[,-i], option, decomp = TRUE, scalar =  scalar)$Fit.p,
                            "obj" =  compute_obj(X,W,W_c, U[,-i], V[,-i], option, scalar =  scalar))

      diff = new.obj.fit - init.obj.fit
      if(diff < min.increase)
      {
        min.index <- i
        min.increase <- diff
      }
    }
    if(!is.infinite(min.increase))
    {
      message("Calling the next iteration")
      message("Dropping factor number: ", min.index)
      #print(dim(U))
      #print(dim(V))
      if(min.index > ncol(U))
      {
        message("You see the problem...")
        stopifnot(ncol(U) < min.index)
      }
      return(DropFactorsByFit(X,W,W_c,U[,-min.index],V[,-min.index], maxK, option,scalar = scalar))
    }else
    {
      return(list("U"=U, "V" = V, "K" = ncol(V)))
    }

  }else
  {
    return(list("U"=U, "V" = V, "K" = ncol(V)))
  }
}
#Beta test function- can we pick out the best factors by how they contribute (or don't) to the objective?
#This will balance in sparsity better than PVE does.
#' A recursive function to drop factors that by so doing improve (that is, reduce) the overall objective
#'
#' @param X full data matrix
#' @param W standard errors
#' @param W_c Whitening covariance matrix
#' @param U Predicted U
#' @param V Predicted V
#' @param maxK maximum K to include
#' @param option GLEANER options object
#' @param minK fewest number of factors to allow; default is 0
#'
#' @return updated U,V,K
#' @export

DropFactorsByObjective <- function(X,W,W_c,U,V, maxK, option, minK = 0, scalar = 1)
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
  init.obj <- compute_obj(X,W,W_c, as.matrix(U), as.matrix(V), option, scalar =  scalar)
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
    new.obj <- compute_obj(X,W,W_c, U[,-i], V[,-i], option, loglik = TRUE, scalar =  scalar)
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
    message("Dropping factor: ", min.index)
    message("Min.obj:", min.obj)
    message("original obj:", init.obj)

    return(DropFactorsByObjective(X,W,W_c,U[,-min.index],V[,-min.index], maxK, option, minK=minK))
  }
}


PruneFactorsByObjective <- function(X,W,U,V, minK, option)
{
  keep.dat <- DropFactorsByObjective(X,W,W_c,U,V, minK, option)
  return(list("U"))
}


#' Order a return object so factors in order by PVE
#'
#' @param ret the object return from ALS, with a V, U and PVE object ($)
#' @param X,W to calculate the new PVE (its possible V,U have changed..)
#'
#' @return return.dat, ret object updated by order of PVE
#' @export
#'
OrderEverythingByPVE <- function(ret,X,W,W_c,options,terminal.attempt = FALSE,...)
{
  return.dat <- ret
  if(all(ret$V == 0) | all(ret$U == 0))
  {
	  message("returned matrices empty,no reorder will be performed.")
	  return(ret)
  }
  if(ncol(ret$V) == 1 & ncol(ret$U) == 1)
  {
    message("Only one column remains, no reorder to do.")
    return(ret)
  }
  #PercentVarEx <- function(x,v,u, K=NULL)

  pve.ord = getPVEOrder(as.matrix(X), pm(ret$V), pm(ret$U),W,W_c,options,...)

  #pve=PercentVarEx(as.matrix(X)*as.matrix(W) %*% W_c, pm(ret$V), pm(ret$U),...)
  #pveo <- order(pve, decreasing=TRUE)
  return.dat$V<- pm(ret$V)[,pve.ord]
  return.dat$U <- pm(ret$U)[,pve.ord]
  return.dat$PVE <- PercentVarEx(X,return.dat$V, return.dat$U,W,W_c, options)
  #save(X,W,W_c, ret, return.dat, file= "/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/run_scripts/final_sim_june_2024/debug.RData")
  #if((!all(sort(return.dat$PVE, decreasing=TRUE) == return.dat$PVE)) & !terminal.attempt) #terminal attempt ensures it onlyl tries to correct once.
  #{
  #  message("Not a descending setting. Odd. Try again.")
  #  print(return.dat$PVE)
  #  print(sort(return.dat$PVE, decreasing=TRUE))

  #  return(OrderEverythingByPVE(return.dat,X,W,W_c,terminal.attempt = TRUE,...))
  #}
  #need go recalculate since its order dependent- adding them in one by one- and they aren't necessarily totally independent.
  return.dat
}

