#!/usr/bin/env Rscript
#This script has been repursposed as our write_out script. Should be called that, but for convenience staying with this...
#Manually extract runs from this directory to check for enrichment, see what we get.
#s(alphas, lambdas, names, run_stats,output)
writeFactorMatrices <- function(specs1, specs2, studies, snps, run_stats, out_prefix)
{
  entry <- c()
  for(i in specs1) { for(j in specs2) { entry <- c(entry, paste0(i, " ", j)) }}
  run_list <- entry
  for(r in run_list)
  {
    index <- which(entry == r)
    curr_run <-  run_stats[[index]]
    writeFactorMatrix(studies, snps, curr_run, out_prefix)

  }
}

writeFactorMatrix <- function(studies, snps, curr_run,identifier, out_prefix)
{
  f_dat <- data.frame("rownames" = studies, curr_run$V)
  #write_tsv(f_dat, paste0(out_prefix, gsub(" ", "_", r), ".factors.txt"))
  write.table(f_dat, file = paste0(out_prefix, gsub(" ", "_", identifier), ".factors.txt"), sep = "\t", quote = FALSE, row.names=FALSE)
  if(length(curr_run) > 1)
  {
    if(nrow(curr_run$U) != length(snps))
    {
      message("Unable to complete L matrix, unsure why")
      print(head(curr_run$U))
      #readline()
    } else{
    l_dat <- data.frame("snps" = snps,curr_run$U)
    write.csv(x = l_dat, file = paste0(out_prefix, gsub(" ", "_", identifier), ".loadings.txt"), quote = FALSE, row.names=FALSE)
    }

  }else{
    message("No L matrix generated for this run, F too sparse.")
  }
}

#' Debugging tool to print output as tracking progress
#'
#' @param f vector of sparsity parameters
#' @param s vector of sparsity parameters
#' @param t vector of sparsity parameters
#'
#' @return none, just prints out
#' @export
reportSimilarities <- function(f, s, t)
{
  print("Alphas:")
  print(c(f$alpha, s$alpha, t$alpha))
  print("Lambdas:")
  print(c(f$lambda, s$lambda, t$lambda))

  print("Ks:")
  print(c(f$K, s$K, t$K))
}



#' Report update in job and write to log
#'
#' @param text What to write out
#' @param options Option list. Key values are 'ofile', which specifies where to write the log (or not to write to it at all if NULL),
#' and "V", which indcates how verbose to me
#'
#' @return None
#' @export
#'
updateLog <- function(text, options = NULL)
{
  #added back in some key stuff.
  verbosity=options$V
  if(is.null(options))
  {
    ofile = NULL
  }else
  {
    ofile = options$logu
  }

  if(is.null(ofile))
  {
    message(text)
  }else if (verbosity == 1)
  {
    #logr::log_print(text,console = TRUE, file = ofile)
    message(text)
  } else{
    #logr::log_print(text,console = FALSE, file = ofile)
    message(text)
  }

}

