
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ecocbo

<!-- badges: start -->

[![R-CMD-check](https://github.com/arturoSP/ecocbo/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/arturoSP/ecocbo/actions/workflows/R-CMD-check.yaml)
[![License](https://img.shields.io/badge/License-GPL3-blue.svg)](https://github.com/arturoSP/ecocbo/blob/master/LICENSE.md)

<!-- badges: end -->

## Calculating Optimum Sampling Effort in Community Ecology

**ecocbo** is an R package that helps scientists calculate the optimum
sampling effort for community ecology projects. The package is based on
the principles developed in the
[SSP](https://github.com/edlinguerra/SSP) package, which simulates
ecological communities by extracting and using parameters that control
the simulation. The simulated communities are then compared with
PERMANOVA to estimate their components of variation and consequently the
optimal sampling effort.

**ecocbo** is a valuable tool for scientists who need to design
efficient sampling plans. The package can help scientists to save time
and money by ensuring that they collect the minimum amount of data
necessary to achieve their research goals.

## Installation

``` r
# You can easilly obtain 'ecocbo' from CRAN:
install.packages("ecocbo")

# Alternatively, you can install the development version of ecocbo from [GitHub](https://github.com/):
# install.packages("devtools")
devtools::install_github("arturoSP/ecocbo")
```

## Example

This is a basic example which shows you how to use the different
functions in the package:

### Prepare the data

``` r
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

### Calculate statistical power

``` r
betaResult <- sim_beta(simH0 = simH0Dat, simHa = simHaDat, 
                       n = 10, m = 3, k = 20, alpha = 0.05,
                       useParallel = F)
betaResult
#> Power at different sampling efforts (m x n):
#>       n = 2 n = 3 n = 4 n = 5 n = 6 n = 7 n = 8 n = 9 n = 10
#> m = 2  0.20  0.55  0.72  0.82  0.93  0.97  1.00  0.98      1
#> m = 3  0.43  0.50  0.98  0.92  1.00  1.00  0.98  1.00      1
```

### Plot the power progression as sampling increases.

``` r
plot_power(data = betaResult, n = NULL, m = 3, method = "power")
```

<figure>
<img src="man/figures/plotm3n4.png"
alt="This plot will look different in every simulation" />
<figcaption aria-hidden="true">This plot will look different in every
simulation</figcaption>
</figure>

### Calculate components of variation.

``` r
compVar <- scompvar(data = betaResult)
compVar
#>    compVarA  compVarR
#> 1 0.0717776 0.3303087
```

### Determine optimal sampling effort

The sampling effort can be evaluated depending on an economic budget
(ct) or desired precision level (multSE), depending on the proposed
parameter, the function will calculate optimal values for number of
treatments (bOpt) and replicates (nOpt).

``` r
cboCost <- sim_cbo(comp.var = compVar, ct = 20000, ck = 100, cj = 2500)
cboPrecision <- sim_cbo(comp.var = compVar, multSE = 0.10, ck = 100, cj = 2500)
cboCost
#>   nOpt mOpt
#> 1   10    5
cboPrecision
#>   nOpt mOpt
#> 1   10   10
```

## R packages required for running ecocbo

- Required: ggplot2, ggpubr, sampling, stats, vegan
- Suggested: SSP, knitr, rmarkdown, testthat

## Sponsored by

<img src="man/figures/logoCONACYT.png" height="121" />
<img src="man/figures/logoENES.png" height="121" />
