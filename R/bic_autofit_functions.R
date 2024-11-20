################################################################################################################################
## Ashton Omdahl, November 2024
## Machinery for performing model fitting in GLEANR using the BIC, as well as tracking relevant data.
##
################################################################################################################################


########################################### Functions for tracking convergence/model fitting progression ###########################################

#' Update sparsity parameter list with next selection
#' Notify user if there was a large jump in the sparsity parameters.
#' @param prev.list list of sparsity parameters so far
#' @param new new parameter proposed
#' @param errorReport Setting if we want to give extra information
#'
#' @return updated list of sparsity parameters
#' @export
#'
UpdateAndCheckSparsityParam <-  function(prev.list, new, errorReport = FALSE)
{
 if(!is.null(prev.list))
    {

	 prev <- prev.list[length(prev.list)]
	  new <- new[1]
	if( errorReport & !is.null(prev.list) > 1 & (abs(new-prev)/prev > 1.5))
	  {
	    message("WARNING: large jump in sparsity parameter encountered....")
	    message("Current: ", prev)
	    message("New: ", new)
	  }
	}
  return(c(prev.list, new))
}

#' Updated tracking of the optimal learned parameters from gleaner
#'
#' @param U learned at settings minimizing BIC_alpha at iteration `iter_tag`. Should NOT have 0 columns dropped, but be in K dimension the same as V
#' @param V learned at settings minimizing BIC_lambda at iteration `iter_tag` Should NOT have 0 columns dropped, but be in K dimension the same as U
#' @param alpha setting minimizing BIC_alpha at iteration `iter_tag`
#' @param lambda setting minimizing BIC_lambda at iteration `iter_tag`
#' @param bic_a minimal score at iteration `iter_tag`
#' @param bic_l minimal score at iteration `iter_tag`
#' @param bic_a_dat: all the data for that score, including the scaling term and the other pieces
#' @param bic_l_dat: all the data for that score, scaling term and other pieces
#' @param iter_tag
#'
#' @return a list containing the optimal learned parameters (U,V,alpha,lambda) and the associated BIC
#' @export
#'
#' @examples
trackMinParam <- function(min.dat,U,V,alpha,lambda,bic_a,bic_l,bic_a_dat, bic_l_dat, iter_tag, init=FALSE)
{
  if(init)
  {
    message("Initializing parameter tracking")
    return(list("optimal.u"=NA, "optimal.v"=NA,
                "alpha"=NA, "lambda"=NA, "bic_a"=Inf,
                "bic_l"=Inf, min_sum=Inf, iter="",
                "bic_a_dat"=NA, "bic_l_dat"=NA))
  }else {
    stopifnot(ncol(U) == ncol(V)) #this should be true.
    if(length(min.dat)==1)
    {
      message("Error- should not occur")
      quit()
    } else {
      #Special case: all of U or all of V is empty.
      if(CheckUEmpty(U) & is.infinite(bic_l)) #U is empty and we have no score for V
      {
        #we just zeroed out of all of U
        bic_l = bic_a_dat$fit.term/1 + 0 +  bic_a_dat$addends
        #In all cases, scaling term is 1 b/c p= 0, and y is scaled to have variance 1
        if(is.na(bic_l_dat) & !is.na(bic_a_dat)[[1]]) #subcase- we never got to update v
        {
          bic_l_dat <- list("bic.list" = bic_a,
                            "fit.term" = bic_a_dat$fit.term,
                            "df.term"=  0,
                            "fit.scaler"= 1, #what would the fit be with no predictors? either 1 or less than 1, depending on p and what version; depends on the bic version ug.
                            "addends" = bic_a_dat$addends, "n.coef"= 0) #not sure how to handle this
        }
      }
      if(CheckVEmpty(V) & is.infinite(bic_a)) #V is empty and we have no score for a
      {
        bic_a = bic_l_dat$fit.term/1 + 0 +  bic_l_dat$addends
        #In all cases, scaling term is 1 b/c p= 0, and y is scaled to have variance 1

        if(is.na(bic_a_dat) & !is.na(bic_l_dat)) #subcase- we never got to update v
        {
          bic_a_dat <- list("bic.list" = bic_l,
                            "fit.term" = bic_l_dat$fit.term,
                            "df.term"=  0,
                            "fit.scaler"= 1,
                            "addends" =bic_l_dat$addends, "n.coef"= 0) #not sure how to handle this
        }

      }
      if(CheckUEmpty(U) & CheckVEmpty(V))
      {
        #if both of these guys are 0'd out
        message("Unexpected edge case")
      }
      new_bic = bic_a+bic_l
      if(new_bic < min.dat$min_sum)
      {
        message("Found better setting, updating now...")
        return(list("optimal.u"=U, "optimal.v"=V,
                    "alpha"=alpha, "lambda"=lambda, "bic_a"=bic_a,
                    "bic_l"=bic_l, "min_sum"=new_bic, "iter"=iter_tag,
                    "bic_a_dat"=bic_a_dat, "bic_l_dat"=bic_l_dat))

      }else
      {
        return(min.dat)
      }
    }
  }

}


#' Wrapper for returning program data
#'
#' @param min.dat tracker for minimal/optimal data
#' @param rec.dat record data from iteration to iteratino
#' @param mat.fit u.fit data object
#' @param conv.options record of convergence data
#' @param general.options list of program options.
#'
#' @return list with optimal U and V, optimal lambda and alpha, and some other helpful stuff.
#' @export
#'
#' @examples
returnCompletionDat <- function(min.dat, rec.dat, mat.fit, conv.options, general.options)
{
  if(is.null(mat.fit) | length(mat.fit) == 1) #may need to check these conditions
  {
    mat.fit <- list()
    mat.fit$resid_var <- NA
  }
  #K
  if(all(is.na(min.dat$optimal.v)) | CheckVEmpty(min.dat$optimal.v) | CheckUEmpty(min.dat$optimal.u))
  {
    kout=NA
  }
  else
  {
    kout=ifelse(!is.null(ncol(min.dat$optimal.v)), ncol(min.dat$optimal.v),NA)
  }
  list("optimal.v" = min.dat$optimal.v,"resid.var" = mat.fit$resid_var,
       "rec.dat" = rec.dat, "lambda"=min.dat$lambda, "alpha"=min.dat$alpha, "options" = general.options,
       "K"= kout, "alpha.path" = rec.dat$alpha.s, "lambda.path" = rec.dat$lambda.s, "optimal.u" = min.dat$optimal.u, "convergence.options" = conv.options,
       "min.dat"=min.dat)
}


#' Helper function to update record data showing sparsity and bic changes
#'
#' @param rec.dat current record data, a list of lists
#' @param mat.fit the fit of the most recent step, either u.fit or v.fit
#' @param matrix_on which matrix you just updated, either "U" or "V"
#' @param iter which iteration you're on
#' @param X NxM matrix of SNP effect sizes
#' @param W NxM matrix of SNP weights (1/SE)
#' @param W_c MxM matrix for decorrelating transformation
#' @param optimal_complement the opposite matrix of matrix_on (actual matrix object)
#' @param option
#'
#' @return updated rec.dat
#' @export
#'
#' @examples
updateRecDat <- function(rec.dat, mat.fit, matrix_on, iter,X,W,W_c, optimal_complement,option)
{

  if(matrix_on == "U")
  {
    curr.lambda <- rec.dat$lambda.s[length(rec.dat$lambda.s)]
    rec.dat$alphas <- UpdateAndCheckSparsityParam(rec.dat$alphas, mat.fit$alpha.sel, errorReport = TRUE)
    rec.dat$bic.a <- c(rec.dat$bic.a,mat.fit$bic); rec.dat$alpha.s <- c(rec.dat$alpha.s,mat.fit$alpha.sel)
    rec.dat$U_sparsities = c(rec.dat$U_sparsities, matrixSparsity(mat.fit$U, ncol(X)));
    rec.dat$sparsity.obj[[paste0("U_", iter)]] <- computeObjIntermediate(X, W,W_c, mat.fit$U, optimal_complement, option,mat.fit$alpha.sel,curr.lambda, decomp = TRUE, loglik = TRUE) #should happen before drop empty colums
    rec.dat$Ks <- c(rec.dat$Ks, ncol(DropEmptyColumns(mat.fit$U,option)$updatedMatrix))
    rec.dat$bic_sum <- c(rec.dat$bic_sum,mat.fit$bic + rec.dat$bic.l[length(rec.dat$bic.l)] )
    }
  if(matrix_on == "V")
  {
    curr.alpha <- rec.dat$alpha.s[length(rec.dat$alpha.s)]
    rec.dat$lambda.s <- c(rec.dat$lambda.s,mat.fit$lambda.sel) #selected lambda to path
    rec.dat$V_sparsities = c(rec.dat$V_sparsities, matrixSparsity(mat.fit$V, ncol(X))); #Sparsity of V
    rec.dat$lambdas <- UpdateAndCheckSparsityParam(rec.dat$lambdas, mat.fit$lambda.sel, errorReport = TRUE); rec.dat$bic.l <- c(rec.dat$bic.l,mat.fit$bic) #lambdas and bic scores
    #computeObjIntermediate <- function(X,W,W_c,U,V,option,alpha,lambda,...)computeObjIntermediate <- function(X,W,W_c,U,V,option,alpha,lambda,...)
    rec.dat$sparsity.obj[[paste0("V_", iter)]] <- computeObjIntermediate(X, W,W_c, optimal_complement, mat.fit$V,option,curr.alpha,mat.fit$lambda.sel, decomp = TRUE, loglik = TRUE) #Objective calculation
    rec.dat$Ks <- c(rec.dat$Ks, ncol(DropEmptyColumns(mat.fit$V,option)$updatedMatrix)) #Track K HEREEEe
    rec.dat$bic_sum <- c(rec.dat$bic_sum,mat.fit$bic + rec.dat$bic.a[length(rec.dat$bic.a)] )
  }

  return(rec.dat)
}

#Convergence criteria for the BIC ssearch
#Converges when K is unchanging from one run to the next, and the percentage size change in the alpha/lambda paramters is less than 5%
#Might consider making this more generous- if it stays on the same log scale, then that is probably good enough....
checkConvergenceBICSearch <- function(index, record.data, conv.perc.thresh = 0.05, conv.mode = "SEPARATE")
{
  if(index > 10)
  {
    message("Late stage convergence, terminate soon....")
  }
  conv.K <- record.data$Ks[[index]] == record.data$Ks[[index-1]]
  if(conv.mode == "COMB.SUM")
  {
    message("Need to test and debug still, but okay")
    EPSILON = 1e-4
    sum.prev = record.data$alpha.s[[index-1]] + record.data$lambda.s[[index-1]]
    sum.curr = record.data$alpha.s[[index]] + record.data$lambda.s[[index]]

    return(  conv.K & (sum.prev - sum.curr) < EPSILON)

  }else if(conv.mode == "BIC.change")
  {
    #Goes until the BIC stops dropped, or the change crosses 0.
    #I think maybe you need a better condition.. if you stop too soon you miss a nice minimum somewhere else.
    sum.bic.score <- record.data$bic.a + record.data$bic.l
    #at which iteration is our search no longer decreasing? TRUE = still decreasing
    delta.bic.score <- sapply(2:length(sum.bic.score), function(i) sum.bic.score[i-1] > sum.bic.score[i])
    li <- length(delta.bic.score)
    #if(any((!delta.bic.score)))
    if(!delta.bic.score[li]) #The last element is NOT decreasing relative to the previous
    {
      #CAse 1- we've only been increasing so far, give it some more time
      if(all(!delta.bic.score)) #all have been increasing
      {
        message("BIC has only increased after minimum iterations. Continuing search")
        return(FALSE)
      }else{      return(TRUE) } #Case 2 we end here, not searching anymore.
    }else
    {
      return(FALSE)
    }


  }  else if(conv.mode == "mat.change")
  {
    if(norm(record.data$Vs[[index]], "f") -  norm(record.data$Vs[[index-1]], "f") < 1e-4){return(TRUE)}
    return(FALSE)
  }
  else
  {
    queries <- c(conv.K,
                 abs(record.data$alpha.s[[index]] - record.data$alpha.s[[index-1]])/record.data$alpha.s[[index-1]] < conv.perc.thresh,
                 abs(record.data$lambda.s[[index]] - record.data$lambda.s[[index-1]])/record.data$lambda.s[[index-1]] < conv.perc.thresh)
    return(all(queries))
  }




}


#' Track alternative source of convergence. Primarily experimental, to identify alternatives to the objective
#'
#' @param X NxM matrix of SNP effect sizes
#' @param W NxM matrix of SNP weights (1/SE)
#' @param W_c MxM matrix for decorrelating transformation
#' @param U NxK matrix of SNP loadings
#' @param V MxK matrix of trait loadings
#' @param index whcih iteration we are currently on
#' @param record.data
#' @param u.alphamaxthe maximum possible sparsity parameter for a given input X, V (||y'X||_INF)
#' @param v.lambdamax  the maximum possible sparsity parameter for a given input X, U (||y'X||_INF)
#' @param option_old current option settings
#'
#' @return track.modes list containing possible convergence criteria, including V and U norms, the objective function, BIC scores, their sum, and the U and V matrices
#' @export
#'
#' @examples
trackAlternativeConvergences <- function(X,W,W_c, U,V, index, record.data, u.alphamax, v.lambdamax, option_old)
{
  track.modes <- list()
  #simple sum mode
  EPSILON = 1e-4
  track.modes[["simple.sum"]] = record.data$alpha.s[[index]] + record.data$lambda.s[[index]]
  track.modes[["U.norm"]] = norm(U, "F")
  track.modes[["V.norm"]] = norm(V, "F")
  track.modes[["U"]] = U #For memory?
  track.modes[["V"]] = V #For memory?
  #scaled sum mode
  track.modes[["lambda.max"]] <-  v.lambdamax
  track.modes[["alpha.max"]] <- u.alphamax
  track.modes[["alpha.scaled"]] <- record.data$alpha.s[[index]]/u.alphamax
  track.modes[["lambda.scaled"]] <- record.data$lambda.s[[index]]/v.lambdamax
  track.modes[["scaled.sum"]] = record.data$alpha.s[[index]]/u.alphamax + record.data$lambda.s[[index]]/v.lambdamax
  track.modes[["bic.alpha"]] = record.data$bic.a[[index]]
  track.modes[["bic.lambda"]] = record.data$bic.l[[index]]

  option <- option_old
  option$fixed_ubiq <- TRUE
  option$alpha1 <- record.data$alpha.s[[index]]
  option$lambda1 <-  record.data$lambda.s[[index]]
  track.modes[["decomp.obj"]] <- compute_obj(X, W, W_c, U, V, option, decomp = TRUE)
  track.modes
}


################# Miscellaneous helper functions ###################################################

#' Estimator of the lasso degrees of freedom as the # of non-zero coefficients
#'
#' @param mat_in to get DF of
#' @param fixed_first if we count the first factor (applies to V only)
#'
#' @return
#' @export
#'
#' @examples
MatrixDFU <- function(mat_in,fixed_first = FALSE)
{
  #this means we are learning U
  #Use the degrees of freedom: number of non-zero coefficients
  #sum(round(mat_in, digits = 4) != 0)
  #if(fixed_first) FF$
  if(fixed_first)  #if(FALSE) FF$
  {
	  sum(mat_in[,-1] != 0)
  }
  else
  {
	  sum(mat_in != 0)
  }
  #This seems like it would favor sparsity?
}

#Empty means all the terms are 0.
DropEmptyColumns <- function(matin,option)
{
  #DTTFT
  message("In the right version of this?")
  matin <- as.matrix(matin)
  c <- colSums(matin != 0)
  ret.option <- option
  if(any(c == 0))
  {
    #print("dropping")
    drops <- which(c == 0)
    if(1 %in% drops & option$fixed_ubiq){
      #message("Ubiqitous factor was removed. ")
      #message("gleaner now lifting the regularization constraint on 1st factor. Keep an eye on the objective.")
      ret.option$fixed_ubiq <- FALSE
    }
    return(list("updatedMatrix" = as.matrix(matin[,-drops]), "updatedOptions"=ret.option))
  }
  return(list("updatedMatrix" = matin, "updatedOptions"=ret.option))
}

#Same as the above, but for magrittr piping (not actually used.)
DropEmptyColumnsPipe <- function(lin,options)
{
  ret.lin <- lin
  ret.lin[[1]] <- DropEmptyColumns(lin[[1]],options)
  return(ret.lin)

}


#' Utilities for loading GLEANR by function, as opposed to the command line tools for reading everything in
#'
#' @param X NxM matrix of SNP effect sizes
#' @param W- weights from Standard errors
#' @param C - cohort overlap estimates
#' @param snp.ids - list of SNP ids, corresponds to order in X
#' @param trait.names - traitnames to use, corresponds to order in X
#' @param K - specify a K if you would like.
#' @param init.mat- which matrix to initialize to
#' @param covar_shrinkage- Which covariance matrix shrinkage type to apply, default is variance based (Strimmer)
#' @param enforce_blocks - Force the matrix to be block matrix
#' @param covar_se
#' @param ...
#'
#' @return a list containing the objects needed to run gleaner, including the decorrelating transformation matrix andthe options matrix.
#' @export
#'
#' @examples
initializeGLEANR <- function(X,W,C,snp.ids, trait.names, K=0, init.mat = "V", covar_shrinkage=-1,enforce_blocks=TRUE,covar_se=NULL, ...)
{
  #A few quick sanity checks:
  if(all(apply(X, 2, var) == 0) | all(apply(X, 1, var) == 0))
  {
    message("Matrix of effect sizes contains no variance. Likely not a valid matrix- program will terminate")
    quit()
  }
  if(all(apply(W, 2, var) == 0))
  {
    warning("Columns of standard errors contain no variance. May not be valid")
  }

  message("This is an archaic initialization; recommend doing away with this...")
  args <- defaultSettings(K=K,init.mat = init.mat,...)
  args$pve_init <- FALSE
  option <- readInSettings(args)
  option$swap <- FALSE
  option$alpha1 <- 1e-10
  option$lambda1 <- 1e-10
  option$block_covar <- 0.2
  output <- args$output
  #log.path <- paste0(args$output, "gwasMF_log.", Sys.Date(), ".txt")
  #lf <- logr::log_open(log.path, show_notes = FALSE)
  #options("logr.compact" = TRUE)
  #options("logr.notes" = FALSE)
  #Read in the hyperparameters to explore
  hp <- readInParamterSpace(args)
  all_ids <-snp.ids; names <- trait.names
  initk <- option$K
  if(initk == 0)
  {
    message('Initializing X to the max -1')
    option$K <- ncol(X)-1
  }
  #Run the bic thing...
  option$V <- FALSE
  option$fixed_ubiqs <- TRUE

  if(enforce_blocks)
  {
    blocks <- create_blocks(C,cor_thr=0.2)
    covar <- blockifyCovarianceMatrix(blocks, C)
  }else
  {
    message("Not enforcing a block structure")
    covar <- C
  }

  #if specified, srhink

  if(covar_shrinkage == "strimmer" | covar_shrinkage == "STRIMMER")
  {
    covar <- strimmerCovShrinkage(args, covar, covar_se, sd.scaling=1)
  }
  else if(covar_shrinkage > -1)
  {
    message("Performing WL shrinkage, as desired")

    covar <- linearShrinkLWSimple(covar, args$WLgamma)

  }
  else{
    message("Not applying shrinkage on covar")
  }
  W_c_dat <- buildWhiteningMatrix(covar, ncol(X),blockify = -1)
  W_c = W_c_dat$W_c
  list("args" = args, "options" = option,
       "hp" = hp, "all_ids" = snp.ids, "names" = trait.names, "initk" = initk, "W_c" = W_c)
}


########################################### Functions calculating the BIC ####################################################
#' Wrapper for all BIC methods under consideration.
#' TODO: In some instances, sklearn_ebic converts input matrices into "dense" matrix object. Want to avoid this.
#' @param fit glmnet object
#' @param covars long.v or long.u
#' @param y long.x- the explanatory data (weighted)
#' @param bic.mode which version to use; options are Zou, dev, sklearn, or the extended variants of them.
#' @return BIC score for each call (list);
#' Object is a list with bic.list (list of scores), fit.term (list of fit terms), df.term (degrees of freedom term), fit.scaler (denominator for fit term)
#' and addends (non-consant vector of additional terms, usually eBIC. Should not include constant single terms.), and n.coef- how many coefficients in the maximally testeed model.
#'
#' @export
#'
calculateBIC <- function(fit, covars,y, bic.mode)
{
    switch(
    bic.mode,
    "sklearn_eBIC"= extendedBIC(sklearnBIC(fit,covars,y), fit),
    "sklearn"= sklearnBIC(fit,covars,y),
    "Zou_eBIC"= extendedBIC(ZouBIC(fit, covars, y), fit, scale = TRUE),
    "Zou"= ZouBIC(fit, covars, y),
    "dev_eBIC"= extendedBIC(BICglm(fit), fit),
    "dev"= BICglm(fit),
    "std"= stdBIC(fit,covars,y),
    sklearnBIC(fit,covars,y)
          )
}

#' Calculate the extended term which accounts for model size, based on work by Chen et al (2008)
#'
#' @param fit glmnet object
#' @param bic.dat object to update with extension
#'
#' @return extended
#' @export
extendedBIC <- function(bic.dat, fit, scale = FALSE)
{
    #p is the number of parameters under consideration
    #n is the
    p = dim(fit$beta)[1] #all the possible parameters to learn
    n = nobs(fit)
    gamma = 1-1/(2 * log(p) / log(n))
    r <- 2 * gamma * lchoose(p, fit$df)
    if(scale)
    {
      r <- (r/nobs(fit))
    }
  ret.dat <-  bic.dat
  ret.dat$bic.list <- bic.dat$bic.list + r
  ret.dat$addends <-  bic.dat$addends + r
  ret.dat

}

#' Calculate the BIC using a vanilla, textbook standard formulation from a GLMNET fit object
#'
#' @param fit the GLMNET fit object
#' @param long.v the covariates matrix, sparse
#' @param long.x the outcome variable
#'
#' @return a vector of BIC scores across multiple different lambda settings
#' @export
#'
stdBIC <- function(fit, long.v, long.x)
{
  all.preds <- glmnet::predict(fit, newx = long.v)
  #I've already extended the dims here, so n = MN
  p = 1
  n = glmnet::nobs(fit)
  BIC= sapply(1:ncol(all.preds), function(i) -2*penLLSimp(n, p, long.x - all.preds[,i]) + log(n) * fit$df[i])

  list("bic.list" = BIC,
       "fit.term" = sapply(1:ncol(all.preds), function(i) -2*penLLSimp(n, p, long.x - all.preds[,i])),
       "df.term"=  log(n) * fit$df,
       "fit.scaler"=1,
       "addends" =0, "n.coef"=nrow(stats::coef(fit)[-1,]))

}

#' BIc calculation for use with glmnet
#'
#' @param fit glmnet object
#'
#' @return BIC score for each call
#' @export
BICglm <- function(fit, bic.mode = ""){
  #based on code from https://stackoverflow.com/questions/40920051/r-getting-aic-bic-likelihood-from-glmnet
  #Note that this isn't exactly the BIC, but the BIC +a constant (2L(saturated model)), but this doesn't impact model selection
  #I revisited this model and explanation on Aug 29, 2024. Should
  tLL <- -stats::deviance(fit) # 2*(loglike_sat - loglike_fit)
  k <- fit$df
  n <- nobs(fit)
  BIC <- log(n)*k - tLL
  BIC

  list("bic.list" = BIC,
       "fit.term" = stats::deviance(fit),
       "df.term"=  log(n)*k,
       "fit.scaler"=1,
       "addends" =0,
       "n.coef"=nrow(stats::coef(fit)[-1,]))
}

#' Calculate BIC based on what is in the sci-kit learn package (see https://scikit-learn.org/stable/modules/linear_model.html#lasso-lars-ic)
#' this differs from the traditional derivation, in that the use the unbiased estimate for s rather than the MLE.
#' Source code is here https://github.com/scikit-learn/scikit-learn/blob/872124551/sklearn/linear_model/_least_angle.py#L1991
#'
#' @param fit glmnet fit object of current regression run across many sparsity settings
#' @param x expanded sparse block matrix corresponding to U or V used in regression by FitVGlobal or FitUGlobal
#' @param y long vector of all stacked, weighted, and transformed SNP effect sizes
#'
#' @return
#' @export
#'
#' @examples
sklearnBIC <- function(fit,x,y)
{
  p = ncol(x)
  n = length(y)
  #OLS fit
  #This is expensive to do afeter taking so long to fit a goshdarn lasso. Need to find an alternative option?
  lm.fit <- glmnet::bigGlm(y=y , x=x, family = "gaussian",intercept = FALSE,trace.it = 1)

  #our fit
  coef = stats::coef(fit)
  #Should be using t <- predict(fit, newx=x)
  yhat=x%*%coef[-1,] #Why is there an intercept at all? all my runs request none.
  residuals = (as.vector(y)- yhat)
  sse = apply(residuals, 2, function(x) sum(x^2))

  if(n==p)
  {
    warning("Unable to calculated BIC when initialized # factors the same as M. Output results will not make sense.")
  }
  self.noise.variance <- stats::deviance(lm.fit)/(n - p)
  scikitlearn.bic <- n*log(2*pi*self.noise.variance) + sse/self.noise.variance + log(n)*fit$df
  list("bic.list" = scikitlearn.bic,
       "fit.term" =sse,
       "df.term"=  log(n)*fit$df,
       "fit.scaler"=self.noise.variance,
  "addends" =0, "n.coef"=p)
  #ADDENDS used to be n*log(2*pi*self.noise.variance), but this is constanat per score so doesn't matter. Can verify empirically.


}

ZouBIC <- function(fit, x, y, var.meth = "mle")
{
  n <- glmnet::nobs(fit)
  p=ncol(x)
  #lm.fit <- predict(fit, newx = x, s= 0)
  lm.fit <- glmnet::bigGlm(x, y, family = "gaussian", intercept = FALSE)
  all.preds <- glmnet::predict(fit, newx = x) #pretty slow step, surprisingly.
  bic.new <- c()
  fit.new <-  c()
  df.term <- c()
  #best.fit <- norm(y-lm.fit, "2")^2/n
  best.fit <- stats::deviance(lm.fit)/n
  if(var.meth == "unbiased")
  {
    best.fit <- stats::deviance(lm.fit)/(n - length(x))
  }
  for(i in 1:length(fit$df))
  {
    curr.fit <- norm(y-all.preds[,i], type ="2")^2/n
    left.term <- curr.fit / best.fit
    fit.new <- c(fit.new, curr.fit)
    right.term <- (log(n) * fit$df[i])/n
    df.term <- c(df.term, right.term)
    bic.new = c(bic.new, left.term + right.term)
  }
  if(is.na(curr.fit))
  {
    message("here")
  }
  list("bic.list" = bic.new,
       "fit.term" =fit.new,
       "df.term"=  (log(n) * fit$df)/n,
       "fit.scaler"=best.fit,
  "addends" =0, "n.coef"=p)

}

extractMinBicDat <- function(bic.dat, min.index)
{
  list("bic.list" = bic.dat$bic.list[min.index],
       "fit.term" = bic.dat$fit.term[min.index],
       "df.term"=  bic.dat$df.term[min.index],
       "fit.scaler"=bic.dat$fit.scaler,
       "addends" = if(length(bic.dat$addends) > 1)
       {bic.dat$addends[min.index]} else {
         bic.dat$addends},
       "n.coef"=bic.dat$n.coef
  )
}


#' Model fitting routine, using the updateUV function
#'
#' @param option list of runtime options
#' @param X NxM matrix of SNP effect sizes
#' @param W
#' @param W_c
#' @param optimal.init V at which to initialize the run, typically the output from a model selection step
#' @param maxK specify how many K max to include. Factors will be removed if they exceed this.
#' @param opath optional path to write out to for saving RData
#' @param no.prune if no pruning of factors is to occur
#' @param reg.elements list with pre-weighted regression variables
#' @param ...
#'
#' @return list object with final U and V and convergence information
#' @export
#'
#' @examples
gwasML_ALS_Routine <- function(option, X, W,W_c, optimal.init, maxK=0, opath = "", no.prune = FALSE,reg.elements=NULL,...)
{

  message(paste0("Starting at k:", option$K))
  option$debug <- TRUE
  reg.run <- update_UV(X, W, W_c, option, preV = optimal.init, reg.elements =reg.elements,... )
  if(no.prune){return(reg.run)}
  if(maxK != 0)
  {
    if(ncol(reg.run$V) > maxK)
    {
      message("More columns identified than specified by the user.")
      message("gleanr now parsing down to desired number of factors")

    }
    if(ncol(reg.run$V) < maxK)
    {
      message("Resulted in fewer than the maximum number of allowed columns.")
      maxK <- ncol(reg.run$V)
    }
  } else {
    maxK <- ncol(reg.run$V) #If maxK is 0, basically ju
  }
  #Note that no pruning takes place if maxK == ncol
    reg.run <- PruneNumberOfFactors(X,W,W_c,reg.run,option$Kmin, maxK, option)
    reg.run <- OrderEverythingByPVE(reg.run,X,W, W_c,option, jointly = FALSE)

  save(reg.run, file = paste0(option$out,opath, "_gleanr_iter.Rdata" ))

return(reg.run)
}


#' Wrapper to programmatically call the BIC-based factorization
#'
#' @param X NxM matrix of SNP effect sizes
#' @param W uncertainty weights (1/SE), same dimensions as X
#' @param snp.ids the names of the SNPs in X, in order
#' @param trait.names the names of the studies in X, in order
#' @param K the target number of K to yield.
#' @param gwasmfiter how many burn-in iterations to run, only used if rep.run is TRUE
#' @param rep.run if should proceed with multiple random initializations. Otherwise, initalize with SVD
#' @param bic.var How to calculate the variance on BIC term. Deprecated as of Feb 21, 2023
#' @param use.init.k To initialize to given K or to full
#' @param init.mat Which matrix to initialize with, V or U. Default is V
#' @param is.sim Specify if this run is a simulation. Default is false
#' @param save.path save results along the way as output (optional)
#'
#' @return
#' @export
gleaner <- function(X,W, snp.ids, trait.names, C = NULL, K=0, gwasmfiter =5, rep.run = FALSE, covar_se=NULL,
                    bic.var= "sklearn", use.init.k = FALSE, init.mat = "V", is.sim = FALSE,
                    save.path = "", scale.mats = FALSE, regression.method = "glmnet", shrinkWL=-1,...)
{
  opath = save.path
  if(is.null(C))
  {
    C = diag(ncol(X))
  }
  d <- initializeGLEANR(X,W,C, snp.ids, trait.names, K=ifelse(use.init.k, K, 0),
                        init.mat=init.mat, covar_shrinkage=shrinkWL,covar_se=covar_se,...) #Either use specified, or prune down as we
  option <- d$options; args <- d$args; hp <- d$hp; all_ids <- d$all_ids; names <- d$namesl; W_c <- d$W_c

  if(is.sim)
  {
    option$Kmin <- K
  }
  option$scale <- scale.mats
  message("scaling set to: ", option$scale)
  option$svd_init <- TRUE; args$svd_init <- TRUE
  K.cap <- K
  option$bic.var <- bic.var
  #Do the bic learning...
  use.optim <- TRUE
  option$regression_method = regression.method
  reg.vect <- prepRegressionElements(X,W,W_c,option)

  if(option$regression_method == "glmnet")
  {
    bic.dat <- getBICMatricesGLMNET(opath,option,X,W,W_c, all_ids, names, min.iter=option$min.bicsearch.iter,reg.elements=reg.vect)
  }else
  {
    bic.dat <- getBICMatrices(opath,option,X,W,W_c, all_ids, names,min.iter=option$min.bicsearch.iter,reg.elements=reg.vect)
  }
  if(is.sim)
  {
    save(bic.dat, file = paste0(save.path, "bic.RData"))
  }
  if(is.null(bic.dat$K) | is.na(bic.dat$alpha) | all(bic.dat$optimal.v == 0) | is.na(bic.dat$K))
  {
    message("No signal detected under current settings. Program will end")
     return(list("V" = matrix(0, nrow=ncol(X), ncol = 1), "U" = matrix(0, nrow=nrow(X), ncol = 1), "initK" = NA, "K" =0, "obj" = c(NA), "obj_change" = c(),
                     "V_sparsities" = c(), "U_sparsities" = c(), "autofit_lambda" = c(), "autofit_alpha"=c(), "mse" = c(),
                     "V_change" = c(), "U_change" = c(), "decomp_obj" ="", "model.loglik" = c(), "Vs"=list(), "Us"=list(), "pve"=c(0)))
  }
  option <- bic.dat$options
  option$V <- TRUE
  option$alpha1 <- bic.dat$alpha #this is used to regularize U
  option$lambda1 <- bic.dat$lambda #This is used to regularize V
  option$K <- bic.dat$K

    ret <- gwasML_ALS_Routine(option, X, W, W_c, bic.dat$optimal.v, maxK=K.cap)
    if(is.sim)
    {
      save(ret, file = paste0(save.path, "global.fit.als.RData"))
    }

  ret
}

#####
#Version optimized for glmnet
#####MAJOR TODO:
## IMPLEMENT  better convergence response [X]


getBICMatricesGLMNET <- function(opath,option,X,W,W_c, all_ids, names, ...)
{

  if(option$K == "GRID")
  {
    gridSearchK(opath, option,X,W,W_c,all_ids,names,...)
  }else
  {
    getBICWorkhorse(opath,option,X,W,W_c, all_ids, names, ...)
  }
}



#########
#function(K, opath, option, X_in, W_in, W_c, all_ids, names, reg.vect,...)
#  bic.dat <- getBICWorkhorse(opath, option, X_in, W_in, W_c, all_ids, names, ...)
getBICWorkhorse <- function(opath,option,X,W,W_c, all_ids, names, min.iter = 5, max.iter = 100, burn.in.iter = 0,
                            init.mat = NULL, rg_ref=NULL,reg.elements=NULL,INCR.LIMIT=3)
{
  W_ld = NULL
  conv.options <- list()
  CONV.MODE <- option$param_conv_criteria
  #Data to track
  rec.dat <- list("alphas"=c(), "lambdas"=c(), "bic.a" = c(), "bic.l"=c(), "obj"=c(),
                  "v.sparsity" = c(), "u.sparsity"=c(), "iter"=c(), "sd.sparsity.u" = c(), "sd.sparsity.v" = c(),
                  "alpha.s" = c(), "lambda.s" = c(), "Ks" = c(), "Vs" = list(),"Us"=list(), "X_hat" = c(), "sparsity.obj" <- list(),
                  "bic_sum"=c(Inf))

  burn.in.sparsity <- DefineSparsitySpaceInit(X, W, W_c, NULL, option, burn.in = burn.in.iter, rg_ref = rg_ref,reg.elements=reg.elements) #If this finds one with NA, cut them out here, and reset K; we want to check pve here too.
  option$regression_method = "glmnet" #Just in case it wasn't set earlier

  if(option$bic.var == "NONE") #This is no longer used- tweaking settings on the variance in BIC since there isn't agreement on this in practice
    {
      return(list("optimal.v" = burn.in.sparsity$V_burn,"resid.var" = NA,
                  "rec.dat" = rec.dat, "lambda"=1e-20, "alpha"=1e-20, "options" = option,
                  "K"= burn.in.sparsity$new.k, "alpha.path" = NA, "lambda.path" = NA, "optimal.u" = NA, "convergence.options" = NA))
      #We don't need U passed in unless it is request to innitialize U, but not considering for th emoment.
    }

  #Specify key parameters to track throughout;
  optimal.u <- NA; optimal.v <- NA
  curr.alpha <- NA; curr.lambda <- NA
  curr.bic.a <- Inf; curr.bic.l <- Inf;
  u.bic.dat <- NA; v.bic.dat <- NA
  min.dat <- trackMinParam(init=TRUE)
  #If we wish to initialize with U rather than V
  if(option$u_init != "")
  {
    optimal.u <- burn.in.sparsity$U_burn
    option$K <- burn.in.sparsity$new.k
    message("initializing with U")
  }else if(!is.null(init.mat)) #Default case initializing with V from a previous iteration
  {
    optimal.v <- init.mat
    option$K <- ncol(optimal.v)
  }else
  {
    optimal.v <- burn.in.sparsity$V_burn
    option$K <- ncol(optimal.v)
  }

  #things to record
  NOT.CONVERGED <- TRUE; i = 1
  while(i < max.iter & NOT.CONVERGED){

    #Tracking all the Us and Vs, for downstream debugging
    rec.dat$Vs[[i]] <- optimal.v
    rec.dat$Us[[i]] <- optimal.u
    #Set settings so it knows to search for alpha, lambda
    option$alpha1 <- NA; option$lambda1 <- NA #important.

    #If we want to initialize with U, start here
    if(option$u_init != "")
    {
      message("initializing via U.")
      v.fits <- FitVWrapper(X,W,W_c, optimal.u,option,reg.elements=reg.elements)
      rec.dat$lambda.s <- c(rec.dat$lambda.s,v.fits$lambda.sel)
      optimal.v <- v.fits$V
      dropped.dat <- DropEmptyColumns(optimal.v,option); optimal.v <- dropped.dat$updatedMatrix; option <-dropped.dat$updatedOptions
      rec.dat$Ks <- c(rec.dat$Ks, ncol(optimal.v)) #this is based on the previous model?, we haven't pruned out empty columns yet?

      #Record new data
      rec.dat$V_sparsities = c(rec.dat$V_sparsities, matrixSparsity(optimal.v, ncol(X)));
      rec.dat$bic.l <- c(rec.dat$bic.l,unlist(v.fits$BIC[,bic.type]))
      #update the parameters for U based on the new V
      u.sparsity <- DefineSparsitySpace(X,W,W_c, optimal.v, "U", option,reg.elements=reg.elements) #Here in case we hit the max.
      alphas <- SelectCoarseSparsityParamsGlobal(u.sparsity, n.points = 15)
      min.dat <- trackMinParam(min.dat,optimal.u,optimal.v,curr.alpha,v.fits$lambda.sel,curr.bic.a,v.fits$bic,u.bic.dat, v.bic.dat, "V0")
    }

    message(paste0("Fitting U, iteration ", i))
    u.fit <- FitUWrapper(X, W, W_c, optimal.v, option,reg.elements=reg.elements)

    #Update tracking data for convergence detection
    curr.alpha <- u.fit$alpha.sel; curr.bic.a <- u.fit$bic; u.bic.dat <- u.fit$bic.dat
    rec.dat <- updateRecDat(rec.dat, u.fit, "U", i, X,W,W_c, optimal.v, option)
    message("verify this update takes place as expected.")
    if(CheckUEmpty(optimal.u)) {
      message("debug")
    }
    min.dat <- trackMinParam(min.dat,u.fit$U,optimal.v,curr.alpha,curr.lambda,curr.bic.a,curr.bic.l,u.bic.dat, v.bic.dat, paste0("U_", i))
    dropped.dat <- DropEmptyColumns(u.fit$U,option); optimal.u <- dropped.dat$updatedMatrix; option <-dropped.dat$updatedOptions
    #Do a check- is it empty?
    if(CheckUEmpty(optimal.u)) {
      message("U has no signal; ending model selection phase");
        return(returnCompletionDat(min.dat, rec.dat, u.fit,conv.options, option))
      }

    #fit V now
    message(paste0("Fitting V, iteration ", i))
    v.fits <- FitVWrapper(X,W,W_c, optimal.u,option,reg.elements=reg.elements)

    #Update the tracking data for downstream debugging
    rec.dat <- updateRecDat(rec.dat, v.fits, "V", i, X,W,W_c, optimal.u, option)
     #Update tracking data for convergence detection
    optimal.v <- v.fits$V
    curr.lambda <- v.fits$lambda.sel; curr.bic.l <- v.fits$bic; v.bic.dat <- v.fits$bic.dat
    min.dat <- trackMinParam(min.dat,optimal.u,optimal.v,curr.alpha,curr.lambda,curr.bic.a,curr.bic.l, u.bic.dat, v.bic.dat, paste0("V_", i))

    if(CheckVEmpty(optimal.v)) {
      message("V has no signal; ending model selection phase");
      return(returnCompletionDat(min.dat, rec.dat, u.fit,conv.options, option))
    }
    if(ncol(optimal.u) == ncol(optimal.v))
    {
      rec.dat$X_hat <- c(rec.dat$X_hat, norm(optimal.u %*% t(optimal.v), type = 'F'))
    }else
    {
      message("Matrices are unaligned. Visit this case to understand why.")
      print(dim(optimal.v))
      print(dim(optimal.u))
      rec.dat$X_hat <- c(rec.dat$X_hat, NA)
    }

    #Check for convergence
    if(i > min.iter){ #Require at least a minimum number of iterations

      NOT.CONVERGED <- !checkConvergenceBICSearch(i, rec.dat, conv.mode = CONV.MODE) #returns true if convergence is reached  #Default CONV.MODE is BIC.change
    }

    #Store some convergence details, for debugging
    conv.options[[i]] <- trackAlternativeConvergences(X,W,W_c, optimal.u,optimal.v, i, rec.dat, u.fit$sparsity_space, v.fits$sparsity_space, option)
    dropped.dat <- DropEmptyColumns(v.fits$V,option); optimal.v <- dropped.dat$updatedMatrix; option <-dropped.dat$updatedOptions
    #Drop empty columns- this must happen before next U iteration so we don't have an empty column, but after scoring things so matrices are still conformable
    i = i+1
  }
  #Finished parameter selection, cleanup.
  final.index <- i-1

  if(ncol(min.dat$optimal.u) != ncol(min.dat$optimal.v))
  {
    message("Warning: final top performing matrices weren't aligned. Aligning now. Note to developer: check this logical condition.")
    aligned <- AlignFactorMatrices(X,W,min.dat$optimal.u,min.dat$optimal.v);
    min.dat$optimal.u <- aligned$U; min.dat$optimal.v<- aligned$V
  }
  #This returns all the data from the last iteration
  returnCompletionDat(min.dat, rec.dat, u.fit,conv.options, option)

}




