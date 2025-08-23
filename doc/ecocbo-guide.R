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

## ----step0, eval=FALSE--------------------------------------------------------
# # Load data and adjust it.
# data(epiDat)
# 
# simResults <- prep_data(data = epiDat, type = "counts", Sest.method = "average",
#                         cases = 5, N = 100, M = 10,
#                         n = 5, k = 30,
#                         transformation = "none", method = "bray",
#                         dummy = FALSE, useParallel = TRUE,
#                         model = "single.factor")

## ----step1--------------------------------------------------------------------
compVar <- scompvar(data = simResults)
compVar


## ----step2--------------------------------------------------------------------
betaResult <- sim_beta(simResults, alpha = 0.05)
betaResult


## ----step31-------------------------------------------------------------------
cboCost <- sim_cbo(betaResult, cn = 75)
cboCost

## ----step4--------------------------------------------------------------------
plot_power(data = betaResult, n = NULL, method = "power")


