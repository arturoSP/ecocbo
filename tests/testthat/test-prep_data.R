data("epiDat")
type = "counts"
Sest.method = "average"
cases = 5
N = 50
sites = 10
n = 7
m = 5
k = 20
transformation = "none"
method = "bray"
dummy = FALSE
useParallel = FALSE
model = "single.factor"

result <- prep_data(data = epiDat, type, Sest.method,
                    cases, N, sites,
                    n, m, k,
                    transformation, method,
                    dummy, useParallel,
                    model)

test_that("prep_data() returns an object with certain characteristics",{
  testthat::expect_s3_class(result, "ecocbo_data")
  testthat::expect_equal(ncol(result$Results), 8)
  testthat::expect_error(prep_data(epiDat, n = 1, m, k, useParallel = FALSE), "larger")
  testthat::expect_error(prep_data(epiDat, n = n, m = 1, k, useParallel = FALSE), "larger")
  testthat::expect_equal(max(result$Results[,2]), k)
})

test_that("no error in documentation", {
  expect_silent(help("prep_data"))
})
