data("simResults")

test_that("sim_beta returns an object with certain characteristics",{
  alpha = 0.05
  results <- sim_beta(simResults, alpha)

  testthat::expect_s3_class(results, "ecocbo_beta")
  testthat::expect_equal(nrow(results$Power),
                         (max(simResults$Results[,3])-1) * (max(simResults$Results[,4])-1))
  testthat::expect_error(sim_beta(simResults, alpha = 2), "smaller")
})
