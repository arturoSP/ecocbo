data("epiDat")

test_that("calc_dist uses centroid distance when there are exactly 2 groups", {
  set.seed(123)
  toy <- data.frame(
    grp = factor(rep(c("A", "B"), each = 6)),
    sp1 = c(rpois(6, 8), rpois(6, 2)),
    sp2 = c(rpois(6, 1), rpois(6, 7)),
    sp3 = c(rpois(6, 3), rpois(6, 5))
  )

  out <- ecocbo:::calc_dist(toy, return = "both")
  expect_equal(out$n_groups, 2)
  expect_equal(out$ecological_effect, as.numeric(out$centroid_dist_matrix[1, 2]))
})

test_that("calc_dist uses weighted global RMS when groups are > 2", {
  set.seed(123)
  toy <- data.frame(
    grp = factor(rep(c("A", "B", "C"), each = 5)),
    sp1 = c(rpois(5, 8), rpois(5, 2), rpois(5, 4)),
    sp2 = c(rpois(5, 1), rpois(5, 7), rpois(5, 4)),
    sp3 = c(rpois(5, 3), rpois(5, 5), rpois(5, 6))
  )

  out <- ecocbo:::calc_dist(toy, return = "both")
  cent <- as.matrix(out$centroids[, -1, drop = FALSE])
  wg <- as.numeric(table(toy$grp))
  cbar <- colMeans(cent)
  sqd <- rowSums((cent - matrix(cbar, nrow = nrow(cent), ncol = ncol(cent), byrow = TRUE))^2)
  expected <- sqrt(sum(wg * sqd) / sum(wg))

  expect_equal(out$n_groups, 3)
  expect_equal(out$ecological_effect, expected, tolerance = 1e-10)
  expect_equal(dim(out$centroid_dist_matrix), c(3, 3))
})

test_that("dbmanova_oneway exposes inferential components and formulas", {
  toy <- data.frame(
    sp1 = c(8, 7, 9, 2, 1, 3),
    sp2 = c(1, 2, 0, 8, 7, 9),
    sp3 = c(3, 4, 2, 5, 6, 7)
  )
  grp <- factor(rep(c("A", "B"), each = 3))

  out <- ecocbo:::dbmanova_oneway(toy, grp, return = "list")
  expect_equal(out$R2, out$SS_between / out$SS_total)
  expect_equal(
    out$omega2,
    (out$SS_between - out$df_between * out$MS_residual) /
      (out$SS_total + out$MS_residual)
  )
})

test_that("sim_ES returns redesigned columns and stores centroid distances", {
  set.seed(12)
  es <- sim_ES(
    data = epiDat,
    steps = 2,
    type = "counts",
    Sest.method = "average",
    cases = 2,
    N = 20,
    n = 4,
    k = 2,
    transformation = "none",
    method = "bray",
    dummy = FALSE,
    useParallel = FALSE,
    model = "single.factor",
    jitter.base = 0
  )

  required_cols <- c(
    "reduction_level", "k", "m", "n",
    "ecological_effect", "omega2", "R2", "pseudoF",
    "centroid_dist_matrix"
  )
  expect_true(all(required_cols %in% colnames(es)))
  expect_true(all(vapply(es$centroid_dist_matrix, is.matrix, logical(1))))
  expect_true(inherits(es, "effect_size_data"))
  pcoa <- attr(es, "pcoa")
  expect_true(is.data.frame(pcoa))
  expect_true(all(c(
    "reduction_level", "step", "dat_sim", "k", "m", "n",
    "group", "Axis1", "Axis2"
  ) %in% names(pcoa)))
})

test_that("sim_ES blocks nested model in redesigned pipeline", {
  expect_error(
    sim_ES(
      data = epiDat,
      steps = 1,
      cases = 1,
      N = 10,
      n = 3,
      k = 1,
      useParallel = FALSE,
      model = "nested.symmetric"
    ),
    "not implemented"
  )
})

test_that("balanced_sampling legacy signature is preserved for prep_data compatibility", {
  expected <- c(
    "i", "Y", "mm", "nn", "YPU", "H0Sim", "HaSim", "resultsHa",
    "transformation", "method"
  )
  expect_equal(names(formals(ecocbo:::balanced_sampling)), expected)
})

test_that("effect_size_data supports ordination plotting using saved PCoA", {
  set.seed(7)
  es <- sim_ES(
    data = epiDat,
    steps = 1,
    type = "counts",
    Sest.method = "average",
    cases = 1,
    N = 16,
    n = 4,
    k = 1,
    transformation = "none",
    method = "bray",
    dummy = FALSE,
    useParallel = FALSE,
    model = "single.factor",
    jitter.base = 0
  )

  p <- plot(es, type = "ordination")
  expect_s3_class(p, "ggplot")
})
