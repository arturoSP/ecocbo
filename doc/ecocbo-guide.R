## ----setup, include = FALSE---------------------------------------------------
library(ecocbo)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.retina=2,
  fig.align='center',
  fig.width = 7,
  fig.height = 5,
  warning = FALSE,
  message = FALSE
)

## ----eval=FALSE---------------------------------------------------------------
#  # Load data and adjust it.
#  data(epiDat)
#
#  epiH0 <- epiDat
#  epiH0[,"site"] <- as.factor("T0")
#  epiHa <- epiDat
#  epiHa[,"site"] <- as.factor(epiHa[,"site"])
#
#  # Calculate simulation parameters.
#  parH0 <- SSP::assempar(data = epiH0, type = "counts", Sest.method = "average")
#  parHa <- SSP::assempar(data = epiHa, type = "counts", Sest.method = "average")
#
#  # Simulation.
#  simH0Dat <- SSP::simdata(parH0, cases = 3, N = 1000, sites = 1)
#  simHaDat <- SSP::simdata(parHa, cases = 3, N = 100, sites = 10)
#

## -----------------------------------------------------------------------------
simResults <- prep_data(data = epiDat, type = "counts", Sest.method = "average",
                        cases = 5, N = 100, M = 10,
                        n = 5, m = 5, k = 30,
                        transformation = "none", method = "bray",
                        dummy = FALSE, useParallel = TRUE,
                        model = "single.factor")

betaResult <- sim_beta(simResults, alpha = 0.05)
betaResult


## -----------------------------------------------------------------------------
plot_power(data = betaResult, method = "power")


## -----------------------------------------------------------------------------
compVar <- scompvar(data = simResults)
compVar


## -----------------------------------------------------------------------------
cboCost <- Underwood_cbo(comp.var = compVar, budget = 20000, a = 3, ca = 1200, cn = 100, cm = 2500)
cboPrecision <- Underwood_cbo(comp.var = compVar, multSE = 0.10, cn = 100, cm = 2500)
cboCost
cboPrecision


