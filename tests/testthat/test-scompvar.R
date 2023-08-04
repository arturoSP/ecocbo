data("epiBetaR")
scompvar(epiBetaR)

test_that("function works", {
  testthat::expect_no_condition(scompvar(epiBetaR))
  testthat::expect_error(scompvar(epiBetaR, n = 1), "larger")
  testthat::expect_error(scompvar(epiBetaR, m = 8), "larger")
})

test_that("results are reproducible", {
  testthat::expect_equal(scompvar(epiBetaR, n = 4, m = 4),
               data.frame(compVarA = 0.07410969, compVarR = 0.3307416))
})
