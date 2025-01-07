# GLEANR: GWAS latent embeddings accounting for noise and regularization
GLEANER is a GWAS matrix factorization tool to estimate sparse latent pleiotropic genetic factors. Factors map traits to a distribution of SNP effects that may capture biological pathways or mechanisms shared by these traits. This repo contains the `gleanr` R package (in development), which we recommend using in conjunction with the [gleanr_workflow](https://github.com/aomdahl/gleanr_workflow) repository.
The bioRxiv preprint describing the `gleanr` method in detail is avaialable here:

[**Sparse matrix factorization robust to sample sharing across GWAS reveals interpretable genetic components**](https://www.biorxiv.org/content/10.1101/2024.11.12.623313v2).


## Installing GLEANR
This can be done directly from github using the  `devtools` package as follows:
```
devtools::install_github("aomdahl/gleanr")
```
## GLEANR method:
This is an ongoing project to develop a flexible, interpretable, and sparse factorization framework to integrate GWAS data across studies and cohorts. We employ a basic alternating least-squares matrix factoriztion algorithm with sparse priors on learned matrices, while accounting for study uncertainty.
Our approach was inspired by work from Yuan He [here](https://github.com/heyuan7676/ts_eQTLs).

## Running GLEANR
Development of tutorials/vignettes for `gleanr` are ongoing. For a basic interactive use case in `R`, see the [vignette associated with this package](https://github.com/aomdahl/gleanr/blob/main/vignettes/gleanr-basic.Rmd). If you'd like to run `gleanr` directly from the command line (our recommended use), use the script `src/gleanr_run.R` available in the [gleanr_workflow repository](https://github.com/aomdahl/gleanr_workflow) after installing this package to run analysis directly on input matrices of summary statistics.
### GLEANR inputs:
To run GLEANR, a user must provide: 
  - a matrix $B$ of $N$ SNPs by $M$ studies of GWAS effect sizes (e.g. $\beta$'s) (required)
    - Each SNP and trait should have a label, as in the example file [here](https://github.com/aomdahl/gleanr/blob/main/inst/extdata/sim1.effect_sizes.txt)
  - an $N \times M$ matrix of GWAS standard error estimates, with the same order as $B$ (required, example file [here](https://github.com/aomdahl/gleanr/blob/main/inst/extdata/sim1.std_error.txt))
  - an $M \times M$ matrix of estimated correlation due to sample sharing ($C$); this may be estimated using LDSC and should have  (optional, example file [here](https://github.com/aomdahl/gleanr/blob/main/inst/extdata/sim1.c_matrix.txt))
  - an $N \times M$ matrix  of esitmation error correlation due to sample sharing; this will be used to regularize $C$ (optional, example file [here](https://github.com/aomdahl/gleanr/blob/main/inst/extdata/sim1.c_se_matrix.txt))
  - an $M \times 1$ list of trait names corresponding to $M$ (required). This can be used to specify cleaner names for columns in $B$. These should be unique.
  - an $M \times 1$ list of standard deviation estimates across trait Z-scores (optional; only provide if using XT- LDSC to estimate degree of sample sharing)

## Development versions of gleanr (preceeding Nov 2024)
To review development versions of gleanr prior to the reorgnization of this github in Nov. 2024, please see the `gleanr_source_backup` directory in the [gleanr_workflow repository](https://github.com/aomdahl/gleanr_workflow/tree/main/gleanr_source_backup).
