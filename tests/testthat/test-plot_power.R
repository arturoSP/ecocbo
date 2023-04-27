library(ecocbo)
data("epiBetaR")

test_that("plots are plotted", {
  expect_type(plot_power(epiBetaR, n = 4, m = 4, method = "power"), "list")
  expect_no_condition(plot_power(epiBetaR, m = 4, method = "both"))
  expect_no_error(plot_power(epiBetaR, n = 3, m = 4))
  expect_error(plot_power(epiBetaR, n = 1, m = 4, method = "both"), "larger")
})

test_that("calculation of n does not affect graph result", {
  p1 <- plot_power(epiBetaR, m = 2, method = "power")
  p2 <- plot_power(epiBetaR, m = 3, method = "density")

  expect_equal(plot_power(epiBetaR, n = 5, m = 2, method = "power"), p1)
  expect_equal(plot_power(epiBetaR, n = 5, m = 3, method = "density"), p2)
})
