#Helper functions for reading in and parseing code.
#TODO: develop some testing cases for each of these
#' Function to sort by lead column
#'
#' @param tab to sort
#' @param option option vector, where sort is specified
#' @param col which column to sort by
#'
#' @return sorted table
#' @export
#'
#' @examples
quickSort <- function(tab, option, col = 1)
{
  if(!option$sort)
  {
    return(tab)
  }
  #tab <- as.matrix(tab)
  tab[order(tab[,col, with = FALSE], decreasing = TRUE),] #the package version.
  #tab[order(tab[,..col], decreasing = TRUE),] #the debug version
}


#' Prepare possible parameter space for lasso sparsity parametsr
#'
#' @param args  list of arguments
#'
#' @return initial alphas and lambdas
#' @export
#'
#' @examples
readInParamterSpace <- function(args)
{
  #Read in the hyperparameters to explore
  if(args$MAP_autofit > -1 | args$auto_grid_search) #DEPRECATED
  {
    updateLog("Sparsity parameters will be proposed by software.")
    alphas_og <- c(NA)
    lambdas_og <- c(NA)
    if(args$alpha != "" | args$lambdas != "")
    {
      updateLog("Incompatible sparsity settings provided- cannot autofit and use specified. Program will terminate.")
      quit()
    }
  }else
  {
    alphas_og <- as.numeric(scan(text = args$alphas, what = character(), sep = ',', quiet = TRUE))
    lambdas_og <- as.numeric(scan(text = args$lambdas, what = character(), sep = ',', quiet = TRUE))
  }

  return(list("a" = alphas_og, "l" = lambdas_og))

}

#Holdover from previous version. May not be used
cleanUp <- function(matin, type = "beta")
{
  lister <- as.vector(unlist(matin))
  if(type == "beta")
  {
    bad <- is.na(lister)
    lister[bad] <- 0
    bad <- is.infinite(lister)
    lister[bad] <- 0
  }
  if(type == "se")
  {
    bad <- is.infinite(lister)
    lister[bad] <- 1000
    bad <- is.na(lister)
    lister[bad] <- 1
  }
  return(matrix(lister, nrow = nrow(matin)))
}


#' Clean up input GWAS data, including extreme X^2 values and entries with too many NAs
#' TODO: add inheritance for some of these arguments
#' @param X
#' @param W
#' @param id.list
#' @param na_thres
#' @param na.threshold
#'
#' @return a list of updated input items (X,W,id.list) reflecting the drops/filtering
#' @export
#'
#' @examples
matrixGWASQC <- function(X, W, id.list, na_thres = 0.5, na.threshold=0.7)
{
  #Now that we have the SE and the B, go ahead and filter out bad snps....
  #Clean up NAs, etc.
  #Sanity check for number ofo NAs. X*W gives Z scores
  #CHECK  drops by SNP first
  tnames = colnames(X)
  drop.nas.rows <- apply(X*W, 1, function(x) sum(is.na(x)))
  drops = c()
  #if(any(drop.nas.rows == ncol(X)))
  if(any(drop.nas.rows > ncol(X)*na.threshold))
  {
    drops = which(drop.nas.rows > (ncol(X)*na.threshold))
    X <- X[-drops,]
    W <- W[-drops,]
    id.list <- id.list[-drops]

    updateLog(paste0("Removed ", sum(drop.nas.rows > ncol(X)*na.threshold), " variants where ", na.threshold, " of all entries were NA..."))
  }
  X <- as.matrix(X)
  W <- as.matrix(W)
  drop.nas <- apply(X*W,2, function(x) is.na(x))
  #We aren't including this anymore
  #chi.thresh = 300
  chi.thresh = max(X*W^2, na.rm=TRUE)+10
  #drop.chi2 <- apply(X*W,2, function(x) x^2 > chi.thresh)
  #drop.chi2[drop.nas] <- FALSE
  #all.drops <- drop.chi2 | drop.nas
  all.drops <- drop.nas
  #all.drops <- drop.nas
  #do we just zero those ones out or drop them all together?
  #if more than 15% of SNPs for a particular trait are in this category, drop the trait instead
  too.many.drops <- unlist(lapply(1:ncol(all.drops), function(i) sum(drop.nas[,i])))
  drop.cols <- c()
  if(any(too.many.drops > 0.20 * nrow(X))) #greater than 20%
  {
    warning.counts <- which(too.many.drops > (0.20 * nrow(X)))
    drop.cols <- which(too.many.drops > (na_thres * nrow(X)))
    for(w in warning.counts)
    {
      na.proportion <- too.many.drops[w]/nrow(X)
      updateLog(paste0(round(too.many.drops[w]/nrow(X) * 100, digits = 2), "% of SNPs in trait ",tnames[w], " are either missing or invalid."))
    }
    updateLog(paste0("Traits with over ", na_thres*100, "% missing entries will be dropped automatically."))
  }

  #TODO: fix this report, it isn't quite right.
  updateLog("Cells with invalid entries (NA) will be given a weight of 0.")
  removed.cols <- length(drop.cols) * nrow(X)
  updateLog(paste0(sum(all.drops), " out of ", (nrow(all.drops) * ncol(all.drops)), " total cells are invalid and will be 0'd."))
  #updateLog(paste0("  DEPRECATED: Zeroing out ", sum(drop.chi2), " entries with extreme chi^2 stats > ", chi.thresh))
  updateLog(paste0("   Zeroing out ", sum(drop.nas), " entries with NA"))
  W[all.drops] <- 0
  X[drop.nas] <- 0 #technically we don't need to drop them here, if they are being given a weight of 0. But if they are NAs we need to give them a value so they don't killus.
  if(length(drops) == 0){	  drops = 0  }
  if(length(drop.cols) == 0){
    return(list("clean_X" = X, "clean_W" = W, "dropped_rows"=drops, "dropped_cols" =  0, "snp.list"=id.list))
    }
    return(list("clean_X" = X[, -drop.cols], "clean_W" = W[, -drop.cols], "dropped_rows"=drops, "dropped_cols" =  drop.cols, "snp.list"=id.list))

}



#tidy up the column-names. fread has some wierd behavior I don't like....

#Takes care of all the necessary data I/O
#Procedures include:
# Read in files
# Specify the weighting scheme
# "Drop" summary statistics with X^2 stat > 80 (set its weight to 0)
# Scale by sqrt(N), if specified
# Scale by sqrt(LDSC_int) if specified
#@param args- argument object in R
#@return a list containing the values, their corresponding weights, and the SNPS in the matching order
#readInBetas(args$gwas_effects, args)
readInBetas <- function(fpath, option)
{
  `%>%` <- magrittr::`%>%`
  data.table::fread(fpath, check.names = TRUE) %>% quickSort(.,option)
}

#Function to read in a covariance matrix.
#This code was copied from "projection_regression_helper.R", but belongs here
# readInCovariance(args$covar_matrix, names)
#Rows and columns shoulud correspond in order, and names should align
#readInCovariance(args$covar_mat, trait.list, diag_enforce = 1)
readInCovariance <- function(p, name_order, diag_enforce = 1, coerce_threshold=1)
{
  if(p == "" | is.null(p)) {return(NULL)}
  else
  {
    w.in <- as.matrix(data.table::fread(p, check.names = TRUE))
    row.names(w.in) <- colnames(w.in)
    if(!is.na(diag_enforce))
    {
      if(all(diag(as.matrix(w.in)) != diag_enforce))
      {
        stop("Diagonal elements of covariance data matrix don't match expectations. Please ensure rows and columns are in the same order such that diagonal elements correspond to the same trait.")
      }
    }

    if(any(abs(w.in) > coerce_threshold))
    {
      warning_string=paste0("Some entries in this sample sharing correlation matrix are greater than ",  coerce_threshold,
                            ". This is possible is using estimates from XT-LDSC\n",
                            "GLEANR will set these values to +/-0.98. Modify the matrix directly if you want different behavior")
	     warning(warning_string)
	     s <- sign(w.in)
	     w.in[abs(w.in) > 1] <- 0.98#Previosly set to 0.95, here to 0.98
	     w.in <- w.in * s #and then return the sign
    }
    if(!isSymmetric(w.in, tol = 1e-3))
    {
      #could just be due to numerical errors, round it
      if(!isSymmetric(round(w.in, digits = 2)))
      {
        message("WARNING: matrix may not be symmetric. Verify input covariance structure")
      }
    }
    #pick.names <- which(colnames(w.in) %in% name_order)
    #need to match the names
    as.matrix(w.in[name_order, name_order])
  }


}


readInLambdaGC <- function(fpath,X, names)
{
  `%>%` <- magrittr::`%>%`
  #note- this can come in as a matrix or as a simple list
  GC <-  data.table::fread(fpath, check.names = TRUE)
  if(nrow(GC) >= ncol(X)) #we have a single entry for each
  {
    message("Testing that the ordering is the same")
    if(nrow(GC) > ncol(X))
    {
      print(names)
      GC <- dplyr::filter(GC, phenotype %in% names)
      if(nrow(GC) != length(names))
      {
        missing <- names[which(!(names %in% GC$phenotype))]
        print("we are missing")
        print(missing)
      }
    }
    ecol = 1
    if(ncol(GC) == 2)  {   ecol = 2 } #the first one is names
    if(!all(unlist(GC[,1]) == names))
    {
      if(any(!(unlist(GC[,1]) %in% names)))
      {
        message("Mismatch in names, some don't align. Please check this")
        return(NA)
      }
      GC <- GC %>% dplyr::mutate("pheno_rank" = factor(phenotype, levels = names)) %>% arrange(pheno_rank) %>% select(-pheno_rank)
      stopifnot(GC$phenotype == names)
    }
    lambdagc <- as.matrix((do.call("rbind", lapply(1:nrow(X), function(x) unlist(GC[,ecol])))), nrow = nrow(X))
    GC <- lambdagc
  } else { #matrix version.
    GC <- GC %>% dplyr::filter(!row_number() %in% r$drops) %>% dplyr::filter(unlist(.[,1]) %in% all_ids) %>% quickSort(.,args)
    if(!all(all_ids == GC[,1])) {
      message("genomic correction values not in correct order, please advise...")
      GC <- as.matrix(GC %>% select(-1) %>% apply(., 2, as.numeric))
      GC <- GC[,-r$dropped_cols]
    }

  }
  GC
}

SpecifyWeightingScheme <- function(effects, all_ids,all_phenos, args)
{
  effects <- as.matrix(effects[,-1]) %>% orderColumnsByName(., all_phenos,force.ref = args$trait_names)
  #Look at the weighting scheme options...
  if(args$weighting_scheme == "Z" || args$weighting_scheme == "B")
  {
    message("No scaling by standard error will take place. Input to uncertainty being ignored.")
    W <- matrix(1,nrow = nrow(effects), ncol = ncol(effects))
    X <- effects

  } else if(args$weighting_scheme == "B_SE")
  {
    W_se <- data.table::fread(args$uncertainty, check.names = TRUE) %>%
      dplyr::filter(unlist(.[,1]) %in% all_ids) %>% quickSort(.,args)
    stopifnot(all(all_ids == W_se[,1]))
    W_se <- W_se[,-1] %>% orderColumnsByName(., all_phenos, force.ref = args$trait_names)
    userMessage(args$verbosity, paste0("Dims of W_se are now: ", nrow(W_se), " x ", ncol(W_se)))
    W <- 1/ W_se
    X <- effects

  } else if(args$weighting_scheme == "B_MAF")
  {
    message("Scaling by 1/var(MAF)")
    W_maf <- data.table::fread(args$uncertainty, check.names = TRUE) %>%
      dplyr::filter(ids %in% all_ids) %>% arrange(ids) %>% select(-ids)
    W <- 1/matrix(apply(W_maf, 2, function(x) 2*x*(1-x)), nrow = nrow(W_maf), ncol = ncol(W_maf))
    X <- effects
  } else
  {
    message("No form selected. Please try again.")
    quit()
  }
  return(list("X" = X, "W" = W))
}


#' Function to order the columns by matching phenotype names
#'.
#' @param query.mat the matrix to check the order on
#' @param ref the reference order of phenotypes
#' @param force.ref the matrix doesn't yet have names, so just force the names in the order given in ref
#'
#' @return an ordered matrix
#' @export
orderColumnsByName <- function(query.mat, ref, force.ref = "")
{
  ret.mat <- as.matrix(query.mat) #matrix to return with right names
  if(force.ref == "") #if we learned the names on the fly from the effects file. Assumes the names ARE present in the file..
  {
    ret.mat <- ret.mat[,setColOrder(colnames(query.mat), ref)]
  }

  colnames(ret.mat) <- ref
  ret.mat
  #return(query.mat[,..o])
}

#' Title
#'
#' @param query
#' @param ref
#'
#' @return
#' @export
#'
#' @examples
setColOrder <- function(query, ref)
{
  #The reference might not include all the phenotypes, huh?
  ret.list <- 1:length(query)
  #best case- same entries
  if(length(query) == length(ref) & all(query == ref))
  {
    return(ret.list)
  }else if(all(ref %in% query)) #next best- its all there, just need to rearrange
  {
    reorder.q <- order(factor(query, levels = ref))
    return(reorder.q)
  }else #bad case
  {
    message("Phenotype names in query file don't match the reference phenotypes. Check this")
    print(query)
    print("")
    print(ref)
    quit()
  }
}
#' Load the relevant datasets for analysis based on an argument vector
#'
#' @param args Specifies the paths to files and program settings to run. Includes *gwas_effects, weighting_scheme, uncertainty, genomic_correction,* and *covar_matrix arguments*
#'
#' @return A list with entries:
#' "X": sorted GWAS effect sizes,
#' "W": sorted GWAS uncertainty weights (1/SE),
#' "ids": sorted ID order, "trait_names" = names of the traits,
#' "C": the read in covariance matrix and "W_c": the (blockified) matrix to be used for decorrelation
#' @export
readInData <- function(args)
{
  `%>%` <- magrittr::`%>%`
  #Load the effect size data
  rg <- NULL
  effects <- readInBetas(args$gwas_effects, args)

    all_ids <- unlist(effects[,1])
    #Get the trait names out
  if(args$trait_names == "")
  {
    message("No trait names provided. Using the identifiers in the tabular effect data instead.")
    names <- unlist(names(effects)[-1]) %>% make.names(.)
  } else{
    message("Using the provided trait names, and assuming all files have columns in the correct order.")
    message("It is the user's responsibility to verify this.")
    names <- scan(args$trait_names, what = character(), quiet = TRUE) %>% make.names(.)
  }
  if(length(names) > (ncol(effects) - 1))
  {
    message("Error- passed in effect size file and list of phenotype names are of different lengths.")
    quit()
  }
  weighted.dat <- SpecifyWeightingScheme(effects, all_ids,names, args)
  X <- weighted.dat$X; W <- as.matrix(weighted.dat$W);
  if(args$rg_ref != "")
  {
    message("Using an LDSC-rg based V initialization")
    rg <- readInCovariance(args$rg_ref, names)
  }
  if(args$drop_phenos != "")
  {
    message("Dropping data corresponding to names: ", args$drop_phenos)
    drops <- unlist(strsplit(args$drop_phenos, split = ",")[[1]])
    drop.indices <- which(names %in% drops)
    if(length(drop.indices) > 0)
    {
      X <- X[,-drop.indices]; W <- W[,-drop.indices]
      rg <- rg[-drop.indices,-drop.indices]
      names <- names[-drop.indices]
    }
  }

  #remove NAs, extreme values.
  r <- matrixGWASQC(X,W,all_ids)
  X <- r$clean_X;  W <- r$clean_W; all_ids <- r$snp.list
  if(length(r$dropped_cols) > 1 || r$dropped_cols != 0)
  {
    names <- names[-(r$dropped_cols)]
  }

  #Consider moving this into GWASQC, limit studies with fewer than N samples?
  if(args$scale_n != "")
  {
    N <- as.matrix(data.table::fread(args$scale_n, check.names = TRUE) %>% dplyr::filter(!row_number() %in% r$drops) %>%
      dplyr::filter(unlist(.[,1]) %in% all_ids) %>% quickSort(.,args))
    if(args$drop_phenos != "")
    {
      drops <- unlist(strsplit(args$drop_phenos, split = ",")[[1]])
      drop.indices <- which(colnames(N) %in% drops)
      if(length(drop.indices) > 0)
      {
        N <- N[,-drop.indices];
      }

      print(N)
      dim(N)
    }
    if(!all(all_ids == N[,1]))
    {
      #This would occur if we are missing variants.
      message("Counts not provided for all variants. Using the average where variants missing")
      vars <- all_ids[!(all_ids %in% unlist(N[,1]))]
      m <- colMeans(N[,-1])
      ndiff <- nrow(X) - nrow(N)
      pre <- lapply(1:ndiff, function(x) unlist(m))

      first <- do.call("rbind", pre)
      new_rows <- cbind("SNP" = vars,first)
      N <- rbind(N, new_rows) %>% quickSort(.,args)
      stopifnot(all(all_ids == N[,1]))
    }
    N <- as.matrix(N[,-1] %>% apply(., 2, as.numeric)) %>% orderColumnsByName(query.mat = ., ref=names,force.ref = args$trait_names)

    #if(args$drop_phenos != "") {N <- N[,-c(as.integer(drop.indices))] } #Drop }

   if(length(r$dropped_cols) > 1 | r$dropped_cols != 0)
   {
     N <- N[,-r$dropped_cols]
   }
    print(dim(N))
    print(dim(W))
    W <- W * (sqrt(N)) #I have been doing this wrong the whole time....
  }
  if(args$genomic_correction != "")
  {
    message("Including genomic correction in factorization...")
    GC <- readInLambdaGC(args$genomic_correction,X, names) #Check this
    X <- X * (1/sqrt(GC)) # we don't weight by this, we actually adjust the effect sizes by this.
  }

   message("Reading in covariance structure from sample overlap...")
  covar.dat <- SampleOverlapCovarHandler(args, names, X)
  return(list("X" = X, "W" = W, "ids" = all_ids, "trait_names" = names, "C" = covar.dat$C,
              "W_c" = covar.dat$W_c, "rg"=rg, "C_block"=covar.dat$C_block))

}

#' Wrapper to get all the important data extracted and used for cohort overlap covariance adjjustment
#'
#' @param args argument vector from initializing or calling gleanr.
#' @param names list of trait names, in order corresponding to X
#' @param X NxM matrix of SNP effect sizes
#'
#' @return list containing Whitening matrix W_c inverse, Block versino of the C matrix and unblocked version of C matrix.
#' @export
#'
SampleOverlapCovarHandler <- function(args, names, X)
{
  #If its an empty file, just limit to identity
  if(args$covar_matrix == "")
  {
    message("No covariance matrix provided, identity matrix will be used.")
    id = diag(length(names))
    return(list("W_c" =  id, "C_block"=id,"C" = id))
  }

  #Just read in the file
  C <- readInCovariance(args$covar_matrix, names)
  if(any(abs(C) > 1))
  {
    message("Warning- some of the covariance overlap effect estimates > 1. These may severely influence factorization results, especially if you aren't scaling by sample SD")
    message("Consider threshold these to < 1")
  }
  #If we are scaling by the sample standard deviation, assuming we use LDSC input
  sd.scaling = 1
  if(args$sample_sd != "")
  {
    message("Note to user: Verify that the input trait order in the file containing estimates of SD for each study corresponds to the order of GWAS in B and S.")
    sd.df <- data.table::fread(args$sample_sd)
    sd.scaling <- 1/(as.matrix(sd.df$V2) %*% t(as.matrix(sd.df$V2)))
    diag(sd.scaling) <- 1 #make it correlation matrix
    C <- C*sd.scaling
    write.table(C, file = paste0(args$output, "_scaledCovarMatrix.txt"), quote = FALSE, row.names = FALSE)
  }
  blocks <- create_blocks(C,cor_thr=args$block_covar)
  covar <- blockifyCovarianceMatrix(blocks, C)
  if(toupper(args$WLgamma) == "STRIMMER")
  {
    	se.path = args$covar_se_matrix
	#se.path= gsub(args$covar_matrix, pattern="gcov_int.tab.csv", replacement = "gcov_int_se.tab.csv")
    C_se <- readInCovariance(se.path, diag_enforce = NA)
    adjusted.C <- strimmerCovShrinkage(args, covar,C_se, sd.scaling)
    userMessage(args$verbosity, paste0("Norm following adjustment: ", norm(adjusted.C, "F")))
  }else if(toupper(args$WLgamma) == "MLE")
  {
	  message("MLE version isn't implemented, don't use it.")
    adjusted.C <- covar #no change, adjust later
  }else
  {
    adjusted.C <- linearShrinkLWSimple(covar, as.numeric(args$WLgamma))
  }

  whitening.dat <- buildWhiteningMatrix(adjusted.C, ncol(X),blockify = -1)
  write.table(adjusted.C, file = paste0(args$output, "_scaledShrunkBlockCovarMatrix.txt"), quote = FALSE, row.names = FALSE)
  return(list("W_c" =  whitening.dat$W_c, "C_block"=whitening.dat$C_block,"C" = adjusted.C))
}



#' Select the K to initialize with. A few different options, made for testing.
#'
#' @param args a list with settings; key settings are K desired (either a number or a method)
#' @param X_ matrix, weighted and adjusted for covariance
#' @param evals (optional) the egienvalues to evaluate
#'
#' @return a number of K to initialize with
#' @export
#'
#' @examples
selectInitK <- function(args,X_, evals = NULL)
{
  #Optiosn now are: MAX (default), KAISER, GD
  #GD method: from https://arxiv.org/pdf/1305.5870.pdf
  #Kaiser: from ??? TODO
  #MAX: m-1
  #M/2: M/2

  if(is.null(evals) & args$K %in% c("KAISER", "GD"))
  {
    #We've already done spectral decomp, no need to do it again...
    svrun <- svd(X_)
    evals <- svrun$d^2
  }

  corr.based <- c("KAISER-1", "SCREE", "BENTLER", "R2")
  if(is.null(evals) & args$K %in% corr.based)
  {
    #We've already done spectral decomp, no need to do it again...
    cor_struct <- cor2(X_)
    evals <- eigen(cor_struct,symmetric=TRUE, only.values = TRUE)$values
  }

  #if(is.numeric(args$K) & args$K !=0)
  if(args$K %in% paste(1:ncol(X_)))
  {
    message("Using the specified number of factors for initialization, K=", args$K)
    return(as.numeric(args$K))
  } else
  {
    #check this
    k = switch(
      args$K,
      "MAX"= ncol(X_)-1,
      "KAISER"= sum(evals > mean(evals)), #described in https://wires.onlinelibrary.wiley.com/doi/epdf/10.1002/wics.101; better called avergae.
      "AVG"= sum(evals > mean(evals)), #same as tanigawa et al.
      "KAISER-1"=  sum(evals > 1),
      "SCREE"=nFactors::nSeScree(evals,cor=FALSE)$nFactors[1],
      "BENTLER"=nFactors::nBentler(evals,N=nrow(X_))$nFactors,
      "R2"=nFactors::nSeScree(evals,cor=FALSE)$nFactors[2],
      "K/2"=ceiling(ncol(X_)/2),
      "K-2"=ceiling(ncol(X_)/2),
      as.numeric(args$K)
    )
    message("Proceeding with initialized K of ", k)

  }
  if(k == ncol(X_))
  {
    warning("K initialized to the same number as the columns of X. This will result in MN = NK, which will cause some BIC implementations (such as sklearn) to fail")
  }
  return(k)
}

#' Take command-line arguments and convert them into the gleanr options object
#'
#' @param args, a list as you would get from a command line argument parser
#'
#' @return an options list object
#' @export
#'
#' @examples
readInSettings <- function(args)
{
  message("")
  message("------------------------------ INPUT FILE PROCESSING ------------------------------")
 option <- list()
	#careful with this environment passed by reference- we do object copying versions so need to be consisgtent.
	#larger changes required if you're going to use this-- need to change some of the BIc functions
	# option <- listenv::listenv()
  option[['K']] <- args$nfactors
  option[['iter']] <- args$niter
  option[['convF']] <- 0
 option[["nsplits"]] <- as.numeric(args$ncores)
 option[["ncores"]] <- as.numeric(args$ncores)
  option[['conv0']] <- args$converged_obj_change
  option[['ones']] <- FALSE
  option[["plots"]] <- args$overview_plots
  option[['disp']] <- FALSE
  #F matrix initialization
  option[['f_init']] <- args$init_V
  option[['epsilon']] <- as.numeric(args$epsilon)
  option[['u_init']] <- args$init_U
  option[["preinitialize"]] <- FALSE
  option[['carry_coeffs']] <- FALSE
  option[["glmnet"]] <- FALSE
  option[["parallel"]] <- FALSE
  option[["fastReg"]] <- FALSE
  option[["ridge_L"]] <- FALSE
  option[['debug']] <- FALSE #args$debug
  option[["subsample"]] <- args$subsample
  #option[["gls"]] <- ifelse(args$covar_matrix != "", TRUE, FALSE)
  option$gls <- FALSE
  option[["covar"]] <- args$covar_matrix
  option$sort <- args$sort
  if(args$simulation)
  {
    message("specifying minimum K with simulation.")
    option$Kmin <- args$nfactors
  }else
  {
    option$Kmin <- 0
  }
  #Experimental

  #message("Scaling is off by default")
  #had to turn off for simulations, at least for now...
  option$scale <- FALSE
  option$burn.in <- 0
  option$fix.alt.setting <- NA
  option$swap <- FALSE
  option$bic.var <- args$bic_var
  option$svd_init <- args$svd_init
  userMessage(args$verbosity, "Standardizing explanatory variables (W_c(W_s*B)^T) by default")
  option$std_y <- TRUE #Standardizing explanatory variables (W_c(W_s*B)^T) by default
  option$param_conv_criteria <- args$param_conv_criteria
  option$min.bicsearch.iter <- args$min.bic.search.iter
  #Internal use only:
  option$actively_calibrating_sparsity <- FALSE
  option$save_out <- TRUE #save most of the time.
  #This looks alot like object-oriented programming.
  #You should just have this be a proper R object with all the attributes and data you need....


  if(args$regression_method %in% c("penalized", "glmnet", "OLS"))
  {
	 option[["regression_method"]] = args$regression_method #push this through all initializations.
  }else{
	  message("This method isn't recognized. Try penalized or glmnet")
	  quit()
  }
  option[["posF"]] <- args$posF
  option$out <- args$output
  #option$logu <- args$output
  option$logu <- NULL #for not writing to log fie output.
  option[["MAP_autofit"]] <- as.numeric(args$MAP_autofit)
  option$intercept_ubiq <- FALSE
  option$traitSpecificVar <- FALSE
  option$verbosity <- args$verbosity
  option$calibrate_sparsity <- args$scaled_sparsity
  if(args$ncores > 1)
  {
    updateLog(paste("Running in parallel on", args$ncores, "cores"))
    option[["parallel"]] <- TRUE
  }
  option[["ncores"]] <- args$ncores
  option[["fixed_ubiq"]] <- args$fixed_first
  option$std_coef <- args$std_coef
  return(option)
}



#####
## Setting defaults helpful for running elsewhere

defaultSettings <- function(K=0, init.mat = "V", fixed_ubiq= TRUE, conv_objective = 0.001,min_bic_search_iter=5, is_sim=FALSE,verbosity=1, covar_shrinkage=-1 )
{
  args <- defaultInteractiveArgs()
  args$niter <- 200
  args$uncertainty <- ""
  args$gwas_effects <- ""
  args$nfactors <- K
  args$verbosity <- verbosity
  args$scale_n <- ""
  #args$output <- "/scratch16/abattle4/ashton/snp_networks/scratch/testing_gwasMF_code/matrix_simulations/RUN"
  args$output <- "./"
  opath <- "gleanr"
  args$simulation <- is_sim
  args$sort <- FALSE #b/c default for sims.
  args$converged_obj_change <- conv_objective
  args$std_coef <- FALSE
  args$std_y <- TRUE  #Updated 12/04
  args$min.bic.search.iter <- min_bic_search_iter
  if(init.mat == "U")
  {
    args$init_U <- "std"
  }
  args$fixed_first <- fixed_ubiq #trying to see if this help- IT DOENS'T really appear to matter very much.
  args$WLgamma <-covar_shrinkage #Added jan 28, for some reason this was missing.
  args
}

DefaultSeed2Args <- function()
{
  args <- list()
  args$std_coef <- FALSE
  args$gwas_effects <-"/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_GWAS_h2-0.1_rg-0.9/B.tsv"
  args$uncertainty <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_GWAS_h2-0.1_rg-0.9/SE.tsv"
  args$fixed_first <- TRUE
  args$genomic_correction <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_GWAS_h2-0.1_rg-0.9/Lambda_gc.tsv"
  args$nfactors <- 15
  args$trait_names = ""
  args$niter <- 10
  args$alphas <- ""
  args$lambdas <- ""
  args$autofit <- -1
  args$ncores <- 1
  args$weighting_scheme = "B_SE"
  args$output <- "/scratch16/abattle4/ashton/snp_networks/scratch/testing_gwasMF_code/model_selection/bic_autofit/"
  args$converged_obj_change <- 1
  args$scaled_sparsity <- TRUE
  args$posF <- FALSE
  args$init_V <- "ones_eigenvect"
  args$init_U <- ""
  args$epsilon <- 1e-8
  args$verbosity <- 1
  args$scale_n <- ""
  args$MAP_autofit <- -1
  args$regression_method = "penalized"
  args$converged_obj_change <- 0.001
  args$sort <- TRUE
  args
}
YuanSimEasy <- function()
{
  args <- list()
  args$std_coef <- FALSE
  args$covar_matrix = ""
  args$gwas_effects <-"/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/yuan_simulations/Input_tau100_seed1_X.txt"
  args$uncertainty <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/yuan_simulations/Input_tau100_seed1_W.txt"
  args$fixed_first <- TRUE
  args$genomic_correction <- ""
  args$overview_plots <- FALSE
  args$nfactors <- 5
  args$trait_names = ""
  args$niter <- 100
  args$simulation <- TRUE
  args$alphas <- ""
  args$lambdas <- ""
  args$autofit <- -1
  args$ncores <- 1
  args$weighting_scheme = "B_SE"
  args$output <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/yuan_simulations/"
  args$converged_obj_change <- 1
  args$scaled_sparsity <- TRUE
  args$posF <- FALSE
  args$init_V <- "ones_eigenvect"
  args$init_U <- ""
  args$epsilon <- 1e-8
  args$verbosity <- 1
  args$sort <- FALSE
  args$scale_n <- ""
  args$MAP_autofit <- -1
  args$regression_method = "penalized"
  args$converged_obj_change <- 0.05 #this is the percent change from one to the next.
  args$prefix <- ""
  args$bic_var <- "mle"
  args$svd_init <- FALSE
  args
}
defaultInteractiveArgs <- function() #This used to be to update data with the udler data directly.
{
  #Set first
  args <- list()
  args$covar_matrix = ""
  args$gwas_effects <-""
  args$uncertainty <- ""
  args$fixed_first <- TRUE
  args$genomic_correction <- ""
  args$overview_plots <- FALSE
  args$nfactors <- "GRID"
  args$trait_names = ""
  args$niter <- 200
  args$ncores <- 1
  args$weighting_scheme = "B_SE"
  args$output <- ""
  args$converged_obj_change <- 1
  args$epsilon <- 1e-8
  args$converged_obj_change <- 0.001 #this is the percent change from one to the next.
  args$bic_var <- "sklearn_eBIC"
  args$param_conv_criteria <- "BIC.change"

  #redundant with fillDefaultSettings:
  fillDefaultSettings(args)
}

#Old udler file paths:
#args$scale_n <- "/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/processed_data/sample_counts_matrix.tsv"
#args$gwas_effects <-"/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/processed_data/beta_signed_matrix.tsv"
#args$uncertainty <- "/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/processed_data/se_matrix.tsv"
#args$output <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/udler_original/bic_version/"

#' Filling in settings, all relics of a past run.
#' TODO: basically get rid of all this
#'
#' @param curr.args
#'
#' @return updated arguments
#' @export
#'
fillDefaultSettings <- function(curr.args)
{
  curr.args$std_coef <- FALSE
  curr.args$sort <- TRUE #what is this doing? sorts the snps I think...
  curr.args$alphas <- ""
  curr.args$lambdas <- ""
  curr.args$autofit <- -1
  #curr.args$cores <- 1
  curr.args$svd_init <- TRUE
  curr.args$scaled_sparsity <- TRUE
  curr.args$posF <- FALSE
  curr.args$init_V <- "ones_eigenvect"
  curr.args$init_U <- ""
  curr.args$MAP_autofit <- -1
  curr.args$regression_method = "glmnet"
  curr.args$prefix <- ""
  curr.args$scale_n <- ""
  curr.args$output <- curr.args$outdir
  curr.args$simulation <- FALSE
  curr.args$std_y <- TRUE #Updated from FALSE, 12/04
  #TODO: add some checks for missing or cnflicing parameters
  if(is.null(curr.args$rg_ref)) {curr.args$rg_ref <- ""}
  if(is.null(curr.args$ncores)) {curr.args$ncores <- 1}
  curr.args
}


#' Quick readout to report settings, useful for debugging
#'
#' @param argsin the current arguments in.
#'
#' @return
#' @export
#'
#' @examples
writeRunReport <- function(argsin)
{
  message("Running GLEANR, with settings as follows:")
  message("")
  message("------------------------------ INPUT FILES ------------------------------")
  message("Effect sizes: ", basename(argsin$gwas_effects))
  message("Uncertainty estimates: ", basename(argsin$uncertainty))
  message("Cohort overlap adjustment: ", basename(argsin$covar_matrix))
  message("       Block distance: ", as.character(argsin$block_covar))
  message("       Shrinkage factor: ", as.character(argsin$WLgamma))
  message("Genomic correction terms: ", basename(argsin$genomic_correction))
  message("Z-score sample standard deviation: ", basename(argsin$sample_sd))
  message("")

  message("------------------------------ INPUT SETTINGS ------------------------------")
  message("BIC convergence criteria: ", argsin$param_conv_criteria)
  message("BIC method: ", argsin$bic_var)
  message("K init: ", argsin$nfactors)


  message("------------------------------ OUTPUT SETTINGS ------------------------------")
  message("Output directory: ", argsin$output)
}


#' Helper for controling how much output users get when running gleanr
#'
#' @param verbosity level to use- 1 means little, greater than 1 means everything
#' @param message message to share
#' @param thresh threshold to share message. Some messges are essential and so will be shared (thresh = 0)
#'
#' @return
#' @export
#'
#' @examples
userMessage <- function(verbosity, message, thresh=1)
{
  if(verbosity >= thresh)
  {
    message(message)
  }
}
