#Parameter selection helper functions
singlePointGridExpansion <- function(sparsity.params, optimal.sparsity.param, n.points)
{
  sorted.sparsity.params <- sort(sparsity.params, index.return = TRUE)
  ordered.list <- sorted.sparsity.params$x
  optimal.index = which.min(abs(optimal.sparsity.param-sparsity.params))
  sorted.index <- which(sorted.sparsity.params$ix == optimal.index)
  #Had a bug- not exact values being compared. Distance is safer than "==" or dealing with all.equal tolerances
  diffs= ordered.list - optimal.sparsity.param
  optimal.sparsity.param <- ordered.list[which(diffs == min(diffs))] #The optimal sparsity parameter is the one closest to the full precision optimal sparsity parameter.
  if(min(ordered.list) == optimal.sparsity.param)
  {
    above <- ordered.list[sorted.index + 1]
    below <- optimal.sparsity.param - (above - optimal.sparsity.param)
    if(below < 0)
    {
      below = 0
      lower = optimal.sparsity.param * (1/(n.points - 1))
      #seq(0,optimal.sparsity.param, length.out=n.points-1 )
    }
  } else if(max(ordered.list) == optimal.sparsity.param) #its the largest parameter tested
  {
    below <- ordered.list[sorted.index - 1]
    above <- optimal.sparsity.param + (optimal.sparsity.param-below)
  }else {
    #Its bounded- our estimates should be between the one immediately above and below
    above <- ordered.list[sorted.index + 1]
    below <- ordered.list[sorted.index - 1]
  }
  if(is.null(above) | is.null(below))
  {
    print("Error:proposed new paramters are not possible")
    return(NA)
    #quit()
  }
  if(is.na(above) | is.na(below))
  {
    print("Error:proposed new paramters are not possible")
    #quit()
  }
  if(above == below)
  {
    message("Converged on single solution")
    return(below)
  }
  new.list <- seq((below),(above),length.out=n.points+2)
  new.list <- new.list[-c(1, length(new.list))]
  if(any(new.list <= 0))
  {
    print("removing negatives...")
    rep.list <- new.list[which(new.list > 0)]
    if(any(new.list == 0))
    {
      rep.list <- c(min(rep.list)/2,rep.list)
    }
    new.list <- rep.list
  }
  new.list
}

#we don't want to be redundant, so omit the ones that are repeated...
#alphas: the optimal parameters to consisder
#all.a: the full list of parameters from which they are drawn
#n.points: how many new points to pick from.
ProposeSparsityParamsFromGrid <- function(alphas,all.a, n.points)
{
  #N-points will include the original ones...
  if(length(alphas) == 1)
  {
    new.list <- singlePointGridExpansion(all.a[order(all.a)], alphas, n.points)
  }else
  {
    range.a <- range(alphas)
    new.list <- (seq(range.a[1],range.a[2],length = n.points))
  }
  if(any(is.na(new.list)))
  {
    message("detected possible ERROR in proposed parameters")
    print(new.list)
    message("Dropping NAs")
    new.list <- new.list[!is.na(new.list)]
  }
  if(any(new.list < 0))
  {
    message('Detected negatives in proposed parameters; dropping these')
    new.list <- new.list[new.list > 0]
  }
  if(length(new.list) == 0 | all(is.null(new.list)) | all(is.na(new.list)))
  {
    message("Possible error, no terms detected")
    message("Retruning original input.")
    return(alphas)
  }
return(new.list)

}


ProposeNewSparsityParamsFromScore <- function(bic.list,sparsity.params, n.points = 7, no.score = FALSE)
{
  #we want to keep it griddy. So a bit further downstream compare this to the BIC, I think this actually the one that we want, but okay.
  if(no.score)
  {
    #then bic.list is the optimal one; generate fake scores
    fake.scores <- rep(100,length(sparsity.params))
    fake.scores[which(sparsity.params == bic.list)] <- -1
    bic.list <- fake.scores
  }
  optimal.sparsity.param <- sparsity.params[which.min(bic.list)]
  message("optimal sparsity param is ", optimal.sparsity.param)
  ordered.list <- sparsity.params[order(sparsity.params)] #order the list
  #New paradigm: always look above and below,
  #If its the smallest paramter tested
  singlePointGridExpansion(sparsity.params[order(sparsity.params)], optimal.sparsity.param, n.points)
}

#not implemented, too early
ScaleInput <- function(mod.mat, test.method = "preWeight")
{
  s = 1
  test.method = "preWeight"
  if(test.method == "preWeight")
  {
    message("not implemented")
    s <- apply(mod.mat, 2, function(x) norm(x, "2"))
    scaled.m <- sweep(mod.mat,2,s,FUN="/")
    stopifnot(norm(scaled.m[,2], "2") == 1)
    weighted.copies <- lapply(1:nrow(X), function(i) W_c %*% diag(W[i,]) %*% scaled.v)
    long.v <- Matrix::bdiag(weighted.copies) #weight this ish too you silly man.
  }
  if(test.method == "postWeight")
  {
    s = apply(long.v, 2, function(x) norm(x, "2"))
    long.v <- sweep(long.v,2,s,FUN="/")
    #long.v <- long.v / s
    stopifnot(norm(long.v[,3], "2") == 1)
  }
  return(list("s" = s, "updated.mat" = long.mat))
}


ProposeNewSparsityParams <- function(bic.list,sparsity.params, curr.dist, curr.iter, n.points = 7, no.score = FALSE, one.SD.rule = FALSE, no.drop = FALSE)
{
  curr.mode = DensityMode(curr.dist)
  global.min <- min(curr.dist)
  if(length(bic.list) == 1)
  {
    message("No list to choose from- have zeroed all out..?")
    #Go from cuyrrent value to the mode, give a spread
    return(sort(10^seq(log10(sparsity.params),log10(curr.mode),length.out=n.points)))
  }
  if(no.score)
  {
    #then bic.list is the optimal one; generate fake scores
    fake.scores <- rep(100,length(sparsity.params))
    fake.scores[which(sparsity.params == bic.list)] <- -1
    bic.list <- fake.scores
  }
  optimal.index <- selectOptimalScoreIndex(bic.list, sparsity.params, one.SD.rule)
  #cases with redundancy are complex.
  optimal.sparsity.param <- sparsity.params[optimal.index]
  sorted.sparsity.params <- sort(sparsity.params, index.return = TRUE)
  ordered.list <- sorted.sparsity.params$x
  sorted.index <- which(sorted.sparsity.params$ix == optimal.index)
  #what is the index int eh sorted list of my optimal sparsty parameter?
  #New paradigm: always look above and below,
  #If its the smallest paramter tested
  if(min(ordered.list) == optimal.sparsity.param)
  {
    message('best case is the minimum..')
    if(no.drop & curr.iter > 2)
    {
      message("Not allowing a decrease in score if we are after the first iteration")
      message("Method needs to deal with the minimum as the option")
      above <- ordered.list[sorted.index + 1]
      below <- optimal.sparsity.param
    } else
    {
      above <- ordered.list[sorted.index + 1]
      #below <- optimal.sparsity.param - (above - optimal.sparsity.param)
      #simplify this: we are stil searching, so look orders of magnitude
      #below <- 1e-10
      below <- global.min #maybe a better way to do this
      if(below > above)
      {
        message("Global min param of distribution is greater than current one.")
        message("Setting new minimum to 0.1 of current parameter")
        #print(above)
        #print(below)
        below <- optimal.sparsity.param * 0.1
      }
    }

  } else if(max(ordered.list) == optimal.sparsity.param) #its the largest parameter tested
  {
    message('best case is the maximum')
    below <- ordered.list[sorted.index - 1]
    above <- curr.mode
    #new.list <- 10^seq(log10(below),log10(above),length.out=n.points)
  }else {
    #Its bounded- our estimates should be between the one immediately above and below
    above <- ordered.list[sorted.index + 1]
    below <- ordered.list[sorted.index - 1]
    #new.list <- seq(below,above,length.out = n.points)
  }
  if(length(above) > 1 | length(below) < 1)
  {
    message("WHAT is going on...")
    print(above)
    print(below)
    print(bic.list)
    print(sparsity.params)
    readline()
    quit()
  }
  if(is.na(above) | is.na(below))
  {
    print("proposed new paramters are not possible")
    quit()
  }
  if(above == below)
  {
    message("Converged on single solution")
    return(below)
  }
  new.list <- 10^seq(log10(below),log10(above),length.out=n.points)

  #Ensure none of them are less than 0


  if(any(new.list <= 0))
  {
    rep.list <- new.list[which(new.list > 0)]
    if(any(new.list == 0))
    {
      rep.list <- c(min(rep.list)/2,rep.list)
    }
    new.list <- rep.list
  }

  unique(sort(c(optimal.sparsity.param, new.list)))
}


oneSDRule <- function(bics, params)
{
  if(Inf %in% params)
  {
    message("issue here....")
  }
  if(length(unique(bics)) == 1)
  {
    #all the parameters are the same
    message("BIC score for all parameters are same. Likely that all factors have been 0'd out")
    message("Seek recourse.")
    return(which.min(params))
  }
  sd <- sd(bics)
  opt <- min(bics)
  in.range <-bics[(bics > (opt - sd)) & (bics < (opt+sd))]
  if(FALSE)
  {
    print("Optimal")
    print(opt)
    print("SD")
    print(sd)
    print(in.range)
    #top.indices <- which(bics %in% in.range)
    print("all BICs")
    print(bics)
    print("tops")
    print(top.indices)
    print(params)
    #optimal.l <- max(params[top.indices])
    print("Selecting params:")
    print(optimal.l)
  }
  top.indices <- which(bics %in% in.range)
  optimal.l <- max(params[top.indices])
  return(which(params == optimal.l))
}


#This gets the index for the optimal
selectOptimalScoreIndex <- function(bic, params, oneSD, ndigits = 6)
{
  bic <- round(bic, digits = ndigits)
  if(oneSD)#just testing our the one sd rule
  {
    optimal.index <- oneSDRule(bic,params)
  } else
  {
    best.score <- min(bic)
    optimal.index <- which.min(bic)
    if(length(which(bic == best.score)) > 1)
    {
      optimal.index <- SelectBICFromIdenticalScores(bic, params)
      if(length(optimal.index) > 1)
      {
        message("Warning: proposing redundant values.")
        optimal.index <- optimal.index[1]
      }
    }
  }
  optimal.index
}

#Deals with cases if redundant scores.
#If these are at the upper extreme of the parameter list (likely occurs when all terms have been zeroed out), pick the SMALLEST parameter
#if these are at the lower extreme of the parameter list (likely occurs when the terms are fully dense), pick the largest parameter
SelectBICFromIdenticalScores <- function(bic, params)
{
  best.score <- min(bic)
  optimal.index <- which.min(bic)
  #message("Warning- BIC scores are unchanging for certain settings. This is likely evidence of no sparsity.")

  #Choose next lowest bic score index
  i <- sort(bic, index.return = TRUE)

  sorted.bic <- bic[i$ix]
  params.sorted.by.bic <- params[i$ix]
  matching.score.indices <- which(sorted.bic == best.score)

  #If the scores are on the  bigger end of scale
  #if all scores yield the same, its not obvious if we have 0d out or total density. Pick the one closest ot he averagge
  if(length(unique(bic)) == 1)
  {
    message("All scores yield the same BIC. Unclear what to do...")
    mid = abs(params - mean(params))
    optimal.param <- params[which(min(mid) == mid)][1]
  }
  else if(all(min(params.sorted.by.bic[matching.score.indices]) > params.sorted.by.bic[-matching.score.indices]))
  {
    message("Suspect that scores are zeroing out the results, picking the smallest parameter with low BIC")
    optimal.param <- min(params.sorted.by.bic[matching.score.indices])
  } else if(all(max(params.sorted.by.bic[matching.score.indices]) < params.sorted.by.bic[-matching.score.indices]))
  {
    message("Suspect that scores are inducing no sparsity, picking the largest parameter with low BIC")
    optimal.param <- max(params.sorted.by.bic[matching.score.indices])
  }
  else
  {
    #weird case, in the middle. #This means that all of them are equally goood?
    #I this case, I want to pick the most spare one actually
    #print("Beware, unusual case...")
    #optimal.param <- params.sorted.by.bic[matching.score.indices][ceiling(length(matching.score.indices)/2)] #get the middle one
    optimal.param <- max(params.sorted.by.bic[matching.score.indices])
    #print(sorted.bic)
    #print(optimal.param)
    #print(params.sorted.by.bic)
    #print("")
    message("Unusual case with the center, pick the middle. Likely swung too far above or below. Hope is lost :(")

  }

  which(params == optimal.param) #return the optimal index

}
#From new distribution and current list, how to pick the new ones?
#Bic. list: list of BIc scores for all choices
#optimal.sparsity.param- the top parameter chosen
#new.dist: distribution of all the ne lambda parameter space
#@param curr.mode= the mode of the sparsity space based on the current V and U settings
#@return a list of new sparsity points to try out.



