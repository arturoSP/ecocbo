data("epiBetaR")
scompvar(epiBetaR)

test_that("function works", {
  expect_no_condition(scompvar(epiBetaR))
  expect_error(scompvar(epiBetaR, n = 1), "larger")
  expect_error(scompvar(epiBetaR, m = 8), "larger")
})

test_that("results are reproducible", {
  expect_equal(scompvar(epiBetaR, n = 4, m = 4),
               data.frame(compVarA = 0.075459167, compVarR = 0.32864833))
})
