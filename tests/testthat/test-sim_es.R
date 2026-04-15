test_that("nested ecological helper computes hierarchical outputs", {
  skip_if_not_installed("vegan")

  comm <- rbind(
    c(8, 2, 0, 1),
    c(7, 3, 0, 1),
    c(1, 8, 2, 0),
    c(0, 9, 1, 0),
    c(6, 1, 0, 2),
    c(5, 2, 0, 3)
  )

  main <- factor(c("A", "A", "B", "B", "A", "A"))
  nested <- factor(c("s1", "s1", "s2", "s2", "s3", "s3"))

  out <- ecocbo:::.nested_ecological_effect(
    comm = comm,
    main_factor = main,
    nested_factor = nested,
    method = "bray"
  )

  expect_true(is.list(out))
  expect_true(all(
    c(
      "effect_ecol_main",
      "centroids_nested",
      "centroids_main_factor",
      "dist_main_pairwise",
      "disp_nested_within_main",
      "disp_nested_global",
      "pcoa_points"
    ) %in% names(out)
  ))
  expect_true(is.numeric(out$effect_ecol_main))
  expect_true(out$effect_ecol_main >= 0)
  expect_true(nrow(out$centroids_nested) >= nrow(out$centroids_main_factor))
})

test_that("sim_ES single.factor remains backward compatible in core columns", {
  skip_if_not_installed("SSP")
  skip_if_not_installed("vegan")

  data("epiDat", package = "ecocbo")

  res <- sim_ES(
    data = epiDat,
    steps = 1,
    type = "counts",
    Sest.method = "average",
    cases = 1,
    N = 10,
    n = 3,
    k = 1,
    transformation = "none",
    method = "bray",
    dummy = FALSE,
    useParallel = FALSE,
    model = "single.factor",
    jitter.base = 0
  )

  expect_s3_class(res, "effect_size_data")
  expect_true(all(c("ecological_effect", "omega2", "R2", "pseudoF") %in% names(res)))
  expect_false("effect_ecol_main" %in% names(res))
})

test_that("sim_ES nested.symmetric returns hierarchical ecological outputs by step", {
  skip_if_not_installed("SSP")
  skip_if_not_installed("vegan")

  data("dataFish", package = "ecocbo")

  dat <- dataFish[, c(1, 2, 3:8)]
  names(dat)[1:2] <- c("main", "nested")

  # Build a symmetric subset: 2 main levels x 2 nested levels x 3 replicates
  main_levels <- unique(dat$main)[1:2]
  dat <- dat[dat$main %in% main_levels, , drop = FALSE]

  key <- interaction(dat$main, dat$nested, drop = TRUE)
  key_levels <- unique(key)
  key_levels <- key_levels[seq_len(min(4, length(key_levels)))]
  dat <- dat[key %in% key_levels, , drop = FALSE]

  groups <- split(dat, interaction(dat$main, dat$nested, drop = TRUE))
  min_n <- min(vapply(groups, nrow, numeric(1)))
  dat_sym <- do.call(rbind, lapply(groups, function(x) x[seq_len(min(3, min_n)), , drop = FALSE]))
  rownames(dat_sym) <- NULL

  res <- sim_ES(
    data = dat_sym,
    steps = 2,
    type = "cover",
    Sest.method = "average",
    cases = 1,
    N = 6,
    n = 3,
    m = 2,
    k = 1,
    transformation = "none",
    method = "bray",
    dummy = FALSE,
    useParallel = FALSE,
    model = "nested.symmetric",
    jitter.base = 0
  )

  expect_s3_class(res, "effect_size_data")
  expect_true(all(c("effect_ecol_main", "disp_nested_within_main") %in% names(res)))
  expect_true(is.list(attr(res, "nested_ecology")))
  expect_true(length(unique(res$step)) == 3) # steps 0,1,2
  expect_true(all(res$ecological_effect == res$effect_ecol_main))
})
