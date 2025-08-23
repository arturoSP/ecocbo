data("epiBetaR")
data("betaNested")

test_that("plots are plotted", {
  testthat::expect_true(ggplot2::is_ggplot(plot_power(epiBetaR, n = 4, m = 4, method = "power")))
  testthat::expect_condition(plot_power(betaNested, m = 4, method = "both"))
  testthat::expect_no_error(plot_power(epiBetaR, n = 3, m = 4))
  testthat::expect_error(plot_power(epiBetaR, n = 1, m = 4, method = "both"), "larger")
})

test_that("calculation of n does not affect graph result", {
  p1 <- plot_power(epiBetaR, m = 2, method = "power")
  p2 <- plot_power(betaNested, m = 3, method = "density")

  expect_equal(length(plot_power(epiBetaR, n = 5, m = 2, method = "power")$layers),
               length(p1$layers))
  expect_equal(length(plot_power(betaNested, n = 10, m = 3, method = "density")$layers),
               length(p2$layers))
})

test_that("no error in documentation", {
  expect_silent(help("plot_power"))
})
