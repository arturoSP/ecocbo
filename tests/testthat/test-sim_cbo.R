data("epiBetaR")
data("betaNested")

test_that("function works in its two modes", {
  testthat::expect_no_condition(sim_cbo(epiBetaR, cn = 75))
  testthat::expect_no_condition(sim_cbo(betaNested, cn = 75, cm = 125))

  testthat::expect_error(sim_cbo(betaNested, cn = 75), "missing")

  testthat::expect_equal(dim(sim_cbo(betaNested, cn = 75, cm = 125)), c(126,7))
})

test_that("no error in documentation", {
  expect_silent(help("sim_cbo"))
})

