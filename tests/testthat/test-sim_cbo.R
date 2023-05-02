data("epiBetaR")
compVar <- ecocbo::scompvar(epiBetaR)

test_that("function works in its two modes", {
  expect_no_condition(sim_cbo(compVar, ct = 20000, ck = 100, cj = 2500))
  expect_no_condition(sim_cbo(compVar, multSE = 0.14, ck = 100, cj = 2500))

  expect_error(sim_cbo(compVar, ck = 100, cj = 2500), "necessary")

  expect_equal(dim(sim_cbo(compVar, ct = 21000, ck = 150, cj = 2500)), c(1,2))
})
