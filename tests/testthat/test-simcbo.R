data("epiBetaR")
compVar <- ecocbo::scompvar(epiBetaR)

test_that("function works in its two modes", {
  expect_no_condition(simcbo(compVar, ct = 20000, ck = 100, cj = 2500))
  expect_no_condition(simcbo(compVar, multSE = 0.14, ck = 100, cj = 2500))

  expect_error(simcbo(compVar, ck = 100, cj = 2500), "necessary")

  expect_equal(dim(simcbo(compVar, ct = 21000, ck = 150, cj = 2500)), c(1,2))
})


