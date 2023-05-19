test_that("sim_beta returns an object with certain characteristics",{
  n = 10
  m = 3
  k = 20
  alpha = 0.05
  results <- sim_beta(simH0Dat, simHaDat, n, m, k, alpha,
                      transformation = "none", method = "bray", dummy = FALSE,
                      useParallel = FALSE)

  expect_s3_class(results, "ecocbo_beta")
  expect_equal(nrow(results$Power), (m-1) * (n-1))
  expect_error(sim_beta(simH0Dat, simHaDat, n = 1, m, k, alpha), "larger")
  expect_error(sim_beta(simH0Dat, simHaDat, n, m = 1, k, alpha), "larger")
  expect_error(sim_beta(simH0Dat, simHaDat, n, m, k, alpha = 5), "smaller")
})
