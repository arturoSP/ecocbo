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
betaResult <- sim_beta(simH0 = simH0Dat, simHa = simHaDat, 
                       n = 10, m = 3, k = 20, alpha = 0.05,
                       transformation = "none", method = "bray", dummy = FALSE,
                       useParallel = F)
betaResult


## -----------------------------------------------------------------------------
plot_power(data = betaResult, n = NULL, m = 3, method = "both")


## -----------------------------------------------------------------------------
compVar <- scompvar(data = betaResult)
compVar


## -----------------------------------------------------------------------------
cboCost <- sim_cbo(comp.var = compVar, ct = 20000, ck = 100, cj = 2500)
cboPrecision <- sim_cbo(comp.var = compVar, multSE = 0.10, ck = 100, cj = 2500)
cboCost
cboPrecision


