# GLEANR: GWAS latent embeddings accounting for noise and regularization
GLEANER is a GWAS matrix factorization tool to estimate sparse latent pleiotropic genetic factors. Factors map traits to a distribution of SNP effects that may capture biological pathways or mechanisms shared by these traits. This repo contains the `gleanr` R package (in development), which we recommend using in conjunction with the [gleanr_workflow](https://github.com/aomdahl/gleanr_workflow) repository.
The bioRxiv preprint describing the `gleanr` method in detail is avaialable here:

[**Sparse matrix factorization of GWAS summary statistics robust to sample sharing improves detection and interpretation of factors with diverse genetic architectures**](https://www.biorxiv.org/content/10.1101/2024.11.12.623313v1).


## Installing GLEANR
This can be done directly from github using the  `devtools` package as follows:
```
devtools::install_github("aomdahl/gleanr")
```
## GLEANR method:
This is an ongoing project to develop a flexible, interpretable, and sparse factorization framework to integrate GWAS data across studies and cohorts. We employ a basic alternating least-squares matrix factoriztion algorithm with sparse priors on learned matrices, while accounting for study uncertainty.
Our approach was inspired by work from Yuan He [here](https://github.com/heyuan7676/ts_eQTLs).

## Running GLEANR
Tutorials and vignettes will be posted within the next month. If you'd like to try `gleanr` in the meantime, use the script `src/gleaner_run.R` available in the [gleanr_workflow repository](https://github.com/aomdahl/gleanr_workflow) after installing this package to run analysis directly on input matrices of summary statistics.
