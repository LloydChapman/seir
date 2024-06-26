# seir

This package demonstrates use of the [odin](https://mrc-ide.github.io/odin/), [dust](https://mrc-ide.github.io/dust/) and [mcstate](https://mrc-ide.github.io/mcstate/) R packages for running stochastic SIR and SEIR models and fitting them to simulated data [1]. See also the SIR model [odin.dust vignette](https://mrc-ide.github.io/odin.dust/articles/sir_models.html) and [mcstate vignette](https://mrc-ide.github.io/mcstate/articles/sir_models.html).

## Requirements
The package requires the following versions of the [eigen1](https://github.com/mrc-ide/eigen1/), odin, dust, [odin.dust](https://mrc-ide.github.io/odin.dust/) and mcstate packages to be installed:
```
eigen1 0.1.3
odin 1.5.10
dust 0.15.1
odin.dust 0.3.10 
mcstate 0.9.20
```
These can be installed in R with:
```r
remotes::install_github(c(
  "mrc-ide/eigen1@v0.1.3",
  "mrc-ide/odin@v1.5.10",
  "mrc-ide/dust@v0.15.1",
  "mrc-ide/odin.dust@v0.3.10",
  "mrc-ide/mcstate@v0.9.20"))
```
The vignettes use the [coda](https://cran.r-project.org/web/packages/coda/index.html) package, which can be installed from CRAN with:
```r
install.packages("coda")
```

## Installation
Install the package with:

```r
remotes::install_github("LloydChapman/seir", upgrade = FALSE, build_vignettes = TRUE)
```

## References
1. FitzJohn RG, Knock ES, Whittles LK, Perez-Guzman PN, Bhatia S, Guntoro F, Watson OJ, Whittaker C, Ferguson NM, Cori A, Baguelin M, and Lees JA. Reproducible parallel inference and simulation of stochastic state space models using odin, dust, and mcstate. Wellcome Open Research, 5:288, 12 2021. [doi:10.12688/wellcomeopenres.16466.2](https://doi.org/10.12688/wellcomeopenres.16466.2)
