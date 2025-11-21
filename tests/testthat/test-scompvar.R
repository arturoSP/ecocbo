data("simResults")

test_that("function works", {
  testthat::expect_condition(scompvar(simResults))
  testthat::expect_error(scompvar(simResults, n = 1), "larger")
})

test_that("results are reproducible", {
  testthat::expect_equal(round(scompvar(simResults, n = 4, m = 4)[,2],6),
               c(0.333171))
})

test_that("no error in documentation", {
  expect_silent(help("scompvar"))
})
