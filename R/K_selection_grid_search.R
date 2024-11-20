
#' Find the K that minimize the BIC via grid search
#'
#' @param opath to write out the K grid search data object
#' @param option to specify gleaner running settings
#' @param X NxM matrix of SNP effect sizes
#' @param W NxM matrix of SNP weights
#' @param W_c MxM matrix decorrelating transform
#' @param all_ids SNP ids corresponding to X
#' @param names Trait names corresponding to X
#' @param step.limit - what is the maximum number of grid look steps you are willing to take?
#' @param init.range - how broad should the initial search range be? Default is 4
#' @param parallel - options are "none" or "multicore" or "snow"
#' @param ncpu - "how many CPUs to use if going parallel
#' @param ...
#'
#' @return the bic_dat object with the best BIC
#' @export
#'
#' @examples
gridSearchK <- function(opath, option,X,W,W_c,all_ids,names,step.limit = 8,init.range=4,...)
{
  #parallel settings
  parallel="no"
  if(option$ncores > 1){parallel = "multicore"; message("Running in parallel on ", option$ncores, " cores.")}
  #Specify a valid search range
  M = ncol(X)
  k.range <- 1:(M-1)

  #Special case- k.range <= 4:
  #Choose init params, 4:
  init.params <- list("K" = floor(quantile(k.range,seq(0, 1, 1/init.range)))[-1])

  #Edge case- make sure no bad initializations
  if(any(init.params$K <= 0))
  {
    d <- which(init.params$K ==0)
    init.params$K<- init.params$K[-d]
  }
  #Initial grid search
  #optimizeK <- function(K, opath, option, X_in, W_in, W_c, all_ids, names, reg.vect,...)
  init.results <- gridSearch(as.list(init.params$K), option$ncores, opath, option, X, W, W_c, all_ids, names,...)

  #init.results <- paramtest::grid_search(optimizeK, params = init.params, n.iter = 1 , boot = FALSE, bootParams = NULL, parallel = parallel, ncpus = option$ncores,
  #                                       opath = opath,
  #                                       option = option,
  #                                       X_in = X,
  #                                       W_in = W,
  #                                       W_c = W_c,
  #                                       all_ids = all_ids,
  #                                       names = names)

  grid.search.record <- gridSearchRecord(init.results, init.params,  NULL)
  sdv <- gatherSearchData(init.results,k.range,grid.search.record)

  if(sdv$next_params$terminate)
  {
    message("All BICs the same. Grid search is complete.")
    if(all(is.infinite(sdv$query_matrix$BIC)))
    {
      message("All entries are 0'd out. Program will end.")
    }
    return(sdv$min_result)
  }
  save(sdv,grid.search.record, file = paste0(opath, "_K_search.RData"))
  #If we aren't looking at immediatley adjacent K and we haven't exceeded our step limit
  while(length(sdv$query_matrix$K) < step.limit & sdv$k_diff >= 1)
  {
    #Keep looking
    next.grid <- gridSearch(as.list(sdv$next_params$K), option$ncores, opath, option, X, W, W_c, all_ids, names,...)
    #next.grid <- paramtest::grid_search(optimizeK, params = sdv$next_params, n.iter = 1 , boot = FALSE, bootParams = NULL,
    #                                    parallel = parallel, ncpus =option$ncores,
    #                                    opath = opath,
    #                                    option = option,
    #                                    X_in = X,
    #                                    W_in = W,
    #                                    W_c = W_c,
    #                                    all_ids = all_ids,
    #                                    names = names)
    #New paradigm- we need this to update off the full list
    grid.search.record <- gridSearchRecord(next.grid, sdv$next_params, grid.search.record) #use this last entry in the search
    sdv <- gatherSearchData(next.grid,k.range,grid.search.record, curr_grid=sdv$query_matrix,curr_best = sdv$min_result ) #need previous dat here: UPDATE
    if(sdv$next_params$terminate)
    {
      message("All BICs the same. Grid search is complete.")
      if(all(is.infinite(sdv$query_matrix$BIC)))
      {
        message("All entries are 0'd out. Program will end.")
        return(sdv$min_result)
      }
    }

    save(sdv,grid.search.record, file = paste0(opath, "_K_search.RData"))
  #query.matrix <- storeGridResults(next.grid, curr_grid=query.matrix)
  #min.result <- getBestRun(next.grid,curr_best = min.result)
  #For the next iteration:
  #curr_best_K = query.matrix$K[which.min(query.matrix$BIC)]
  #next.params <- chooseNextParams(curr_best_K,query.matrix,k.range)
  #check we aren't looking at super close things
  #k.diff = checkK(c(curr_best_K,next.params$K))
  }
  #save(sdv, file="/scratch16/abattle4/ashton/DEBUG.RData")
  save(sdv,grid.search.record, file = paste0(opath, "_K_search.RData"))
  sdv$min_result
}

#' Manages results of grid search and nominates next parameters to look at
#'
#' @param gs_object - output of grid_search function
#' @param k_range - full range of possible K to consider searching
#' @param ...
#' May be curr_grid=sdv$query_matrix, the current matrix of searches, or
#' curr_best = sdv$query_matrix$min_result, the best bic_dat object so far
#' @return a list containing
#' @export
#'
#' @examples
#' #gatherSearchData(next.grid,k.range, curr_grid=sdv$query_matrix,curr_best = sdv$min_result )
gatherSearchData <- function(gs_object,k_range,grid_record,curr_grid=NULL,curr_best=NULL, version="global")
{
  if(version != "global")
  {
    #Store results so we don't repeat tests
    query.matrix <- storeGridResults(gs_object,curr_grid = curr_grid) #update the list with BICs and Ks
    #Which one is best
    min.result <- getBestRun(gs_object,curr_best= curr_best) #pick out the best result

    #fix starting here
    curr.best.K = query.matrix$K[which.min(query.matrix$BIC)] #and its corresponding K

  }else
  {
    #alternate version
    query.matrix <- list("K"=grid_record$curr_runs$Kinit, "BIC"=grid_record$curr_runs$bic_global)
    #Alternate version
    min.score <- min(grid_record$curr_runs$bic_global)
    min_i = which(grid_record$curr_runs$bic_global == min.score)
    if(length(min_i) > 1)
    {
      min.k.i <- which.min(grid_record$curr_runs$Kinit[min_i])
      min_i= min_i[min.k.i]
    }
    min.result <- grid_record$test_dat[[min_i]]

    curr.best.K <- grid_record$curr_runs$Kinit[min_i]
  }

  #What parameters should we consider next? Binary search approach
  next.params <- chooseNextParams(curr.best.K,query.matrix,k_range) #pick out next params
  #How close are the next Ks we are considering to existing Ks? Make sure we don't repeat anything.
  k.diff = checkK(next.params$K)
  return(list("query_matrix"=query.matrix, "min_result"=min.result,
              "curr_best_K"=curr.best.K,
              "next_params" = next.params, "k_diff"=k.diff))
}

#' Helper function to track the results of run from one to the next. This re-updates all the information each time, so a
#' bit less efficient than it could be, but needed to ensure the BICs are on the right global scale.
#' @param gs_object return from grid search object
#' @param params which settings of K we tried at
#' @param record_obj the object this is being stored in
#'
#' @return list object containing a table of all BIC score information and the data from the test
#' @export
#'
gridSearchRecord <- function(gs_object,params, record_obj)
{
  if(is.null(record_obj))
  {
    record_obj <- list("test_dat"=list(), "curr_runs"=NULL)
  }
  record_obj$test_dat <- c(record_obj$test_dat, gs_object$results) #concatenate all the new information together

  bics <- sapply(record_obj$test_dat, function(x) (x$min.dat$min_sum))
  #bics_global procedure
  #Get the scalar that is max across all for the variance comparison; generally should be M-1.
  bics_global <- getGlobalBICSum(record_obj$test_dat) #need to include the previous information

  alphas <- sapply(record_obj$test_dat, function(x) (x$min.dat$alpha))
  lambdas <- sapply(record_obj$test_dat, function(x) (x$min.dat$lambda))
  K_end <- sapply(record_obj$test_dat, function(x) ncol(x$optimal.v))
  K_end[is.null(K_end)] <- 0

 # record_obj$curr_runs <- rbind(record_obj$curr_runs,
  record_obj$curr_runs <-        data.frame("Kinit"=c(record_obj$curr_runs$Kinit, params$K),
                                           "K_end"=K_end,
                                           "bic_local"=bics,
                                           "bic_global"=bics_global,
                                           "alpha"=alphas,
                                           "lambda"=lambdas) #)

  record_obj
}

#' Figure out which points to test next
#' Currently limited to BINARY SEARCH- only searches 2 points next, at the midpoint between points that have been previously tested (integers only)
#' @param best_K- the K that is currently best
#' @param curr_grid-a list with all the BICs and all the Ks
#' @param all_K- all possible Ks to consider
#'
#' @return  A list with the next K to test, and a logical saying if the search should end
#' @export
#'
#' @examples
chooseNextParams <- function(best_K, curr_grid,all_K)
{
  end = FALSE
  if(all(curr_grid$BIC == curr_grid$BIC[which(curr_grid$K == best_K)]))
    {
      #All the BICs are the same, stop.
      end = TRUE
  }
  testbounds <- getNewTestBounds(best_K, curr_grid,all_K)
  above <- testbounds$upper; below <- testbounds$lower

    #A number in the range above and below:
    above.test <- floor(quantile(best_K:above,probs = c(0.5)))
    below.test <- floor(quantile(below:best_K,probs = c(0.5)))
    test.list <- c(above.test, below.test)

  #Ensure there are no duplicates in new list:
  if(any(duplicated(test.list)))
  {
    message("duplicates in list- not good.")
    test.list <- test.list[!duplicated(test.list)] #Drop duplicates
  }
  #Ensure we aren't repeating previous tests
  if(any(test.list %in% curr_grid$K))
  {
    message("Caught a repeat test, avoiding")
    repeats <- which(test.list %in% curr_grid$K)
    test.list <- test.list[-repeats]
  }

  list("K"=test.list, "terminate"=end)

}

#' Helper function for chooseNextParams to get the bounds of the space to test next.
#'
#' @param best_K K setting which minimzes the BIC
#' @param curr_grid list of previously evaluated K points
#' @param all_K range of valid K for consideration
#'
#' @return upper (above) and lower (below)
#' @export
#'
#' @examples
getNewTestBounds <- function(best_K, curr_grid,all_K)
{
  diffs <- (curr_grid$K - best_K)
  if(best_K == max(curr_grid$K))
  {
    above <- max(all_K)
    #Special case- we are on the edge of what we can try
    if(best_K == max(all_K))  {    above <- max(all_K) - 1    }
    below <-  curr_grid$K[which(diffs == max(diffs[diffs < 0]))]
  } else if(best_K == min(curr_grid$K))
  {
    above <- curr_grid$K[which(diffs == min(diffs[diffs > 0]))]
    below <- min(all_K)
    #Special case- we are on the edge of what we can try
    if(best_K == min(all_K))  {    below <- min(all_K) + 1    }
  }else
  {
    #closest above (omitting the 0 case)
    above <- curr_grid$K[which(diffs == min(diffs[diffs > 0]))]
    below <-  curr_grid$K[which(diffs == max(diffs[diffs < 0]))]
  }
  list("upper"=above, "lower"=below)
}

#' Update search grid with newest search data.
#'
#' @param gs_object list of bic_autofit objects
#' @param curr_grid current data containing a list of Kinits and global BIC scores
#'
#' @return an updated grid with the latest looks.
#' @export
#'
#' @examples
storeGridResults <- function(gs_object, curr_grid=NULL)
{
  tested_k <- gs_object$tests[,2]
  bic_scores <- sapply(gs_object$results, function(x) x$min.dat$min_sum)
  bic_scores <- getGlobalBICSum(gs_object$results)
  ret <- list("K"=tested_k, "BIC"=bic_scores)
  if(!is.null(curr_grid))
  {
    ret <- curr_grid
    ret$K <- c(ret$K,tested_k)
    ret$BIC <- c(ret$BIC,bic_scores)
  }
  ret
}
#
getBestRun <- function(gs_object, curr_best = NULL)
{
  #$$todo
  # Instead of just looking at the best score, we need to rescale all to be on the same variance scale.
  #Do this by getting the fit.scalar for the object with the largest K across the set under consideration
  #then get new BIC scores for each with fit.term/global.fit.scalar + df.term + addends
  #Note that for sklearn, the addends will be problematic since they are subject to the sse of that calculation.
  #They also contain ebic
  #consider just dropping this
  #From all of these rescaled BIC terms, pick the one that minimizes
  #make sure to test in zou case:
  #we pick a terrible zou because it is scaled to look good even though its terrible.
  #list("bic.list" = BIC,
  #     "fit.term" = deviance(fit),
  #     "df.term"=  log(n)*k,
  #     "fit.scaler"=1,
  #     "addends" =0)

  #min_i <- which.min(sapply(gs_object$results, function(x) x$min.dat$min_sum))
  #Its possible that there are matching minimums.
  #In this case, choose the one that has the smaller Kinit score
  find_min_of <- getGlobalBICSum(gs_object$results)
  min_i = which(find_min_of == min(find_min_of, na.rm = TRUE))
  if(length(min_i) > 1)
  {
    min.k.i <- which.min(gs_object$tests$K[min_i])
    min_i= min_i[min.k.i]
  }
  #min_i <- which.min(getGlobalBICSum(gs_object$results))
  new.best <- gs_object$results[[min_i]]
  if(!is.null(curr_best))
  {
    if(curr_best$min.dat$min_sum < new.best$min.dat$min_sum) #if the current best is better than the new best
    {
      new.best <- curr_best #don't change the new best
      message("no change!")
    }
  }
  new.best
}

  checkK <- function(ks)
  {
    #Have to check, if we are testing the same K twice, need to stop
    if(any(duplicated(ks)))
    {
      #return(0)
      ks <- ks[!duplicated(ks)] #Drop duplicates
    }
    if(length(ks) == 0)
    {
      return(0)
    }
    if(length(ks) == 1) #We only have 1 to test
    {
      return(1)
    }
    #However, in the grid any item - itself is 0.
    tests <- expand.grid(ks,ks)
    nonzero.diffs <- abs(tests[,1] - tests[,2])
    z <- which(nonzero.diffs == 0)

    min(nonzero.diffs[-z])

  }
# Load the rgenoud package
# Define a wrapper function for genoud to minimize bic.sum by modifying K
optimizeK <- function(K, opath, option, X_in, W_in, W_c, all_ids, names,...) {
  option$K <- K
  bic.dat <- getBICWorkhorse(opath, option, X_in, W_in, W_c, all_ids, names, ...)
  message("K: ",K, ", BIC: ", bic.dat$min.dat$min_sum)
  return(bic.dat)
}


gridSearch <- function(iter_params, ncpus, opath, option, X_in, W_in, W_c, all_ids, names,...)
{
  ret.list <- list("tests"=data.frame("iter"=1:length(iter_params), "K"=unlist(iter_params)))
  if(ncpus == 1 | length(iter_params) == 1)
  {
    message("Not running in parallel this time...")
    ret.list$results <- lapply(iter_params, function(kinit) optimizeK(kinit, opath, option, X_in, W_in, W_c, all_ids, names,...))
  }else
  {
    message("Now running in parallel across ", ncpus, " cores.")
    ret.list$results <- parallel::mclapply(iter_params, function(kinit)
    {
      optimizeK(kinit, opath, option, X_in, W_in, W_c, all_ids, names,...)
    }, mc.cores = ncpus)
  }
  ret.list
}

getGlobalBICs <- function(fit_vect, source)
{
  #
  max.coef.index <- which.max(sapply(fit_vect, function(x) (x$min.dat[[source]]$n.coef)))
  scalar_global <- fit_vect[[max.coef.index]]$min.dat[[source]]$fit.scaler #
  #Make a check- if this is quite a bit smaller or bigger than other scalars, issue a warning. This might result in some irregular parameter selection:
  all.other.scalars <- sapply(fit_vect, function(x) (x$min.dat[[source]]$fit.scaler))[-max.coef.index]


  #If this value REALLY off from all the others, that is a sign of a potentially bad choice. Could happen if we get into a high-overfitting case:. Throw a warning
  #2 tests- its an order of magnitude off
  checkOverfitVar(scalar_global, all.other.scalars,fit_vect,max.coef.index)
  #Get the scalar that is max across all for the variance comparison; generally should be M-1.
 sapply(fit_vect, function(x) x$min.dat[[source]]$fit.term/scalar_global + x$min.dat[[source]]$df.term + x$min.dat[[source]]$addends)
}


getGlobalBICSum <- function(fit_vect)
{
  round(getGlobalBICs(fit_vect, "bic_a_dat") + getGlobalBICs(fit_vect, "bic_l_dat"), digits=8) #need to round because R counts 1e-12 errors as differences.
}

#' Check if the current variance estimate is really off from all the others and give a warning.
#'
#' @param global - globally minimial variance estimate
#' @param remaining - all other variance estimates
#' @param fit_vect - fit data for all entries
#' @param min_index - inde corresponding to the globally minimal variance estiamte
#'
#' @return nothing if selected variance estimate looks reasonable; otherwise it throws a warning
#' @export
#'
#' @examples
checkOverfitVar <- function(globalmin, remaining,fit_vect,min_index)
{
  ks.all <- sapply(fit_vect, function(x) ncol(x$optimal.v))
  #We are only really concerned about it if it results in many more K being selected than any other option
  if(all(ks.all[min_index] > 2*max(ks.all[-min_index])))  #Got stuck in some local valley with 2 more K than any other option
     {
       if(all(abs(round(log10(globalmin)-log10(remaining),digits=1)) >= 1))
       {
         warning("Selected K has 2x more factors than any other option and a variance estimate at least one order of magnitude different from all others. Proceed with caution, and consider an alternative Kinit strategy.")
       }
     }
     #Alternative test- it doesn't match the distribution of expected variances.
     if(length(unique(remaining)) > 1 & FALSE) #This condition throws a bug if other options are too close- Find something more reliable
     {
       fit = MASS::fitdistr(1/remaining, "gamma")
       prob.not <- min(invgamma::pinvgamma(globalmin,shape = fit$estimate[1],rate = fit$estimate[2], lower.tail = FALSE),
                       invgamma::pinvgamma(globalmin,shape = fit$estimate[1],rate = fit$estimate[2]))
       if(prob.not < 0.005)
       {
         warning("Selected K for variance scaling has variance estimate beyond what we expect by chance (alpha < 0.01).\nProceed with caution, consider an alternative Kinit strategy.")
       }
     }
}

