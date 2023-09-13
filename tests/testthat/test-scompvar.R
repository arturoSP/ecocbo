data("simResults")

test_that("function works", {
  testthat::expect_no_condition(scompvar(simResults))
  testthat::expect_error(scompvar(simResults, n = 1), "larger")
  testthat::expect_error(scompvar(simResults, m = 8), "larger")
})

test_that("results are reproducible", {
  testthat::expect_equal(round(scompvar(simResults, n = 4, m = 4),6),
               data.frame(compVarA = 0.072905, compVarR = 0.330303))
})
