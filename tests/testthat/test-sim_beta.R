data("simResults")

test_that("sim_beta returns an object with certain characteristics",{
  alpha = 0.05
  results <- sim_beta(simResults, alpha)

  testthat::expect_s3_class(results, "ecocbo_beta")
  testthat::expect_equal(nrow(results$Power), 4)
  testthat::expect_error(sim_beta(simResults, alpha = 2), "smaller")
})

test_that("no error in documentation", {
  expect_silent(help("sim_beta"))
})
