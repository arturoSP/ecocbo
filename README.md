
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ecocbo

<!-- badges: start -->

[![R-CMD-check](https://github.com/arturoSP/ecocbo/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/arturoSP/ecocbo/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of ecocbo is to help scientists to calculate an optimum
sampling effort for community ecology projects, based on data from a
pilot study only. This package is based on the principles developed on
[SSP](https://github.com/edlinguerra/SSP), an R package that simulates
ecological communities by extracting and using parameters that control
the simulation. The simulated communities are then compared with
PERMANOVA, to estimate its componrnts of variation and consequently the
optimal sampling effort depending on wether the relevant issue is an
economic budget or a required level of precision.

## Installation

You can install the development version of ecocbo from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("edlinguerra/ecocbo")
```

## Example

This is a basic example which shows you how to use the different
functions in the package:

``` r
# Prepare the data.
# Load data and adjust it.
data(epiDat)

epiH0 <- epiDat
epiH0[,"site"] <- as.factor("T0")
epiHa <- epiDat
epiHa[,"site"] <- as.factor(epiHa[,"site"])

# Calculate simulation parameters.
parH0 <- SSP::assempar(data = epiH0, type = "counts", Sest.method = "average")
parHa <- SSP::assempar(data = epiHa, type = "counts", Sest.method = "average")

# Simulation.
simH0Dat <- SSP::simdata(parH0, cases = 3, N = 1000, sites = 1)
simHaDat <- SSP::simdata(parHa, cases = 3, N = 100, sites = 10)
```

``` r
# Calculate statistical power.
betaResult <- sim_beta(simH0 = simH0Dat, simHa = simHaDat, 
                       n = 10, m = 3, k = 20, alpha = 0.05)
betaResult
#> Power at different sampling efforts (m x n):
#>       n = 2 n = 3 n = 4 n = 5 n = 6 n = 7 n = 8 n = 9 n = 10
#> m = 2 0.167 0.350 0.717 0.750 0.800 0.983 0.983 0.983      1
#> m = 3 0.433 0.817 0.850 0.967 0.983 1.000 1.000 1.000      1
```

``` r
# Plot the power progression as sampling increases.
plot_power(data = betaResult, n = NULL, m = 3, method = "power")
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

``` r
#Calculate components of variation.
compVar <- scompvar(data = betaResult)
compVar
#>   compVarA compVarR
#> 1 0.070026 0.332595
```

``` r
# Determine optimal sampling effort
cboResult <- sim_cbo(comp.var = compVar, ct = 20000, ck = 100, cj = 2500)
cboResult
#>   nOpt bOpt
#> 1   10    5
```
