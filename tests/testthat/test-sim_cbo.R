data("simResults")
compVar <- ecocbo::scompvar(simResults)

test_that("function works in its two modes", {
  testthat::expect_no_condition(sim_cbo(compVar, ct = 20000, ck = 100, cj = 2500))
  testthat::expect_no_condition(sim_cbo(compVar, multSE = 0.14, ck = 100, cj = 2500))

  testthat::expect_error(sim_cbo(compVar, ck = 100, cj = 2500), "necessary")

  testthat::expect_equal(dim(sim_cbo(compVar, ct = 21000, ck = 150, cj = 2500)), c(1,1))
})

test_that("no error in documentation", {
  expect_silent(help("sim_cbo"))
})
