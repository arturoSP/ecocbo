# chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
#
# if (nzchar(chk) && chk == "TRUE") {
#   # use 2 cores in CRAN/Travis/AppVeyor
#   num_workers <- 2L
# } else {
#   # use all cores in devtools::test()
#   num_workers <- parallel::detectCores()
# }

test_that("sim_beta returns an object with certain characteristics",{
  n = 10
  m = 3
  k = 20
  alpha = 0.05
  results <- sim_beta(simH0Dat, simHaDat, n, m, k, alpha,
                      transformation = "none", method = "bray", dummy = FALSE)

  expect_s3_class(results, "ecocbo_beta")
  expect_equal(nrow(results$Power), (m-1) * (n-1))
  expect_error(sim_beta(simH0Dat, simHaDat, n = 1, m, k, alpha), "larger")
  expect_error(sim_beta(simH0Dat, simHaDat, n, m = 1, k, alpha), "larger")
  expect_error(sim_beta(simH0Dat, simHaDat, n, m, k, alpha = 5), "smaller")
})
