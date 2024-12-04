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
Development of tutorials/vignettes for `gleanr` are ongoing. For a basic interactive use case in `R`, see the vignette associated with this package. If you'd like to run `gleanr` directly from the command line (our recommended use), use the script `src/gleaner_run.R` available in the [gleanr_workflow repository](https://github.com/aomdahl/gleanr_workflow) after installing this package to run analysis directly on input matrices of summary statistics.

## Development versions of gleanr (preceeding Nov 2024)
To review development versions of gleanr prior to the reorgnization of this github in Nov. 2024, please see the `gleanr_source_backup` directory in the [gleanr_workflow repositor](https://github.com/aomdahl/gleanr_workflow/tree/main/gleanr_source_backup).
