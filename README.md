# seir

This package demonstrates use of the [odin](https://mrc-ide.github.io/odin/), [dust](https://mrc-ide.github.io/dust/) and [mcstate](https://mrc-ide.github.io/mcstate/) R packages for running stochastic SIR and SEIR models and fitting them to simulated data.

## Requirements
The package requires the following versions of the [eigen1](https://github.com/mrc-ide/eigen1/), odin, dust, [odin.dust](https://mrc-ide.github.io/odin.dust/) and mcstate packages to be installed:
```
eigen1 0.1.1
odin 1.3.2
dust 0.11.24
odin.dust 0.2.16 
mcstate 0.9.0
```
These can be installed in R with:
```r
remotes::install_github(c(
  "mrc-ide/eigen1@v0.1.1",
  "mrc-ide/odin@v1.3.2",
  "mrc-ide/dust@v0.11.24",
  "mrc-ide/odin.dust@v0.2.16",
  "mrc-ide/mcstate@v0.9.0"))
```

## Installation
Install the package with:

```r
remotes::install_github("LloydChapman/seir", upgrade = FALSE, build_vignettes = TRUE)
```
