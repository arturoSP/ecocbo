library(SSP)
library(ecocbo)

test_that("sim_beta returns an object with certain characteristics",{
  data("epiDat")
  epiH0 <- epiDat
  epiH0[,"site"] <- as.factor("T0")
  epiHa <- epiDat
  epiHa[,"site"] <- as.factor(epiHa[,"site"])
  parH0 <- SSP::assempar(data = epiH0, type = "counts", Sest.method = "average")
  parHa <- SSP::assempar(data = epiHa, type = "counts", Sest.method = "average")
  simH0 <- SSP::simdata(parH0, cases = 3, N = 1000, sites = 1)
  simHa <- SSP::simdata(parHa, cases = 3, N = 100, sites = 10)

  n = 10
  m = 3
  k = 20
  alpha = 0.05

  results <- sim_beta(simH0, simHa, n, m, k, alpha)

  expect_s3_class(results, "ecocbo_beta")
  expect_equal(nrow(results$Power), (m-1) * (n-1))
  expect_equal(nrow(results$Results), length(simHa) * k * (m-1) * (n-1))
  expect_error(sim_beta(simH0, simHa, n = 1, m, k, alpha), "larger")
  expect_error(sim_beta(simH0, simHa, n, m = 1, k, alpha), "larger")
  expect_error(sim_beta(simH0, simHa, n, m, k, alpha = 5), "smaller")
})
