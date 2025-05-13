data("simResults")
compVar <- ecocbo::scompvar(simResults)

test_that("function works in its two modes", {
  testthat::expect_no_condition(Underwood_cbo(compVar, budget = 20000, a = 3, ca = 100, cn = 2500))
  testthat::expect_no_condition(Underwood_cbo(compVar, multSE = 0.14, cn = 100, cm = 2500))

  testthat::expect_error(Underwood_cbo(compVar, cn = 100, cm = 2500), "necessary")

  testthat::expect_equal(dim(Underwood_cbo(compVar, budget = 20000, a = 3, ca = 100, cn = 2500)), c(1,1))
})

test_that("no error in documentation", {
  expect_silent(help("sim_cbo"))
})

