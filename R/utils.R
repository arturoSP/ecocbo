#' Sum of squares using Huygen Theorem
#'
#' Calculates sum of squares using Huygen theorem as implemented by Anderson (2014).
#'
#' @param d distance matrix from which the sum of squares will be calculated.
#'
#' @return A numeric vector containing the dimension for the distance matrix, and
#' the value for the sum of squares for the matrix.
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @references Anderson, M. J. (2014). Permutational multivariate analysis of
#' variance (PERMANOVA). Wiley statsref: statistics reference online, 1-15.
#'
#'
#' @keywords internal
#' @noRd
#'

SS <- function(d) {
  ss <- numeric(2)
  ss[1] <- dim(as.matrix(d))[1]
  ss[2] <- sum(d^2) / ss[1]
  return(ss)
}

#' dbMANOVA / PERMANOVA one way (single factor)
#' Calculates observed F and mean squares for the residuals and among sites. This
#' function is a helper for [prep_data()].
#'
#' @param x        Matrix/data.frame with the community data (observations x species).
#' @param factEnv  data.frame or vector with the environmental factor (first column = factor A).
#' @param transformation  Transformation process to reduce dominance. Options include: "none", "square root", "fourth root", "Log (X+1)".
#' @param method   Distance or dissimilarity method for using with vegan::vegdist. Defaults to "bray".
#' @param model    Label for the model: "single.factor".
#' @param return   "table", "list" o "both". What format does the user require for the output? (ANOVA-like table, list of SS/df, or both).
#'
#' @return data.frame showing an ANOVA table or a list with SS, df, MS, pseudoF.
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @references Underwood, A. J. (1997). Experiments in ecology: their logical
#' design and interpretation using analysis of variance. Cambridge university
#' press.
#' @references Anderson, M. J. (2014). Permutational multivariate analysis of
#' variance (PERMANOVA). Wiley statsref: statistics reference online, 1-15.
#'
#' @importFrom vegan vegdist
#' @importFrom stats model.matrix
#'
#' @keywords internal
#' @noRd
#'

## PERMANOVA ----
dbmanova_oneway <- function(
  x,
  factEnv,
  transformation = c("none", "square root", "fourth root", "Log (X+1)"),
  method = "bray",
  model = c("single.factor"),
  return = c("table", "list", "both")
) {
  # --- Validating inputs ---
  transformation <- match.arg(transformation)
  method <- match.arg(method)
  model <- match.arg(model)
  return <- match.arg(return)

  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("`x` must be matrix/data.frame.")
  }
  Xcomm <- as.matrix(x)
  if (any(!is.finite(Xcomm))) {
    stop("`x` contains non finite values (NA/NaN/Inf).")
  }

  # factEnv: function can take either a vector or data.frame; first column is the factor
  if (is.null(dim(factEnv))) {
    fac <- factor(factEnv, exclude = NULL)
    nameA <- deparse(substitute(factEnv))
  } else {
    factEnv <- as.data.frame(factEnv)
    if (ncol(factEnv) < 1L) {
      stop("`factEnv` debe tener al menos 1 columna (factor).")
    }
    fac <- factor(factEnv[[1]], exclude = NULL)
    nameA <- if (!is.null(colnames(factEnv))) colnames(factEnv)[1] else "A"
  }

  if (nrow(Xcomm) != length(fac)) {
    stop("`x` y `factEnv` must have the same number of rows.")
  }

  # --- Transformations to reduce dominance ---
  transform_comm <- function(M, how = "square root") {
    how <- match.arg(how, c("none", "square root", "fourth root", "Log (X+1)"))
    if (how == "none") {
      return(M)
    }
    if (how == "square root") {
      return(sqrt(M))
    }
    if (how == "fourth root") {
      return(M^(1 / 4))
    }
    if (how == "Log (X+1)") return(log1p(M))
  }
  Xtr <- transform_comm(Xcomm, transformation)

  # --- Distance and PCoA (Gower) ---
  d <- vegan::vegdist(Xtr, method = method)
  Dm <- as.matrix(d)

  A <- -0.5 * (Dm * Dm)

  n <- nrow(A)
  rbar <- rowMeans(A)
  cbar <- colMeans(A)
  abar <- mean(A)
  # G_ij = A_ij - rbar_i - cbar_j + abar
  G <- A
  G <- sweep(G, 1, rbar, `-`)
  G <- sweep(G, 2, cbar, `-`)
  G <- G + abar

  # Group statistics
  a <- nlevels(fac)
  ng <- as.integer(table(fac)) # n_g
  # Total sum and trace
  sum_all <- sum(G)
  trG <- sum(diag(G))

  # Sum of block groups (i,j in g)
  # trick: use logical indexes per group (don't create P_A)
  split_idx <- split(seq_len(n), fac)
  sum_in_g <- vapply(
    split_idx,
    function(idx) sum(G[idx, idx, drop = FALSE]),
    numeric(1)
  )

  # trace(P_A G) = sum_g (1/n_g) * sum_{i,j in g} G_ij
  trace_PAG <- sum(sum_in_g / ng)

  # SS y DF
  SSA <- trace_PAG - (1 / n) * sum_all
  SSR <- trG - trace_PAG
  SST <- trG - (1 / n) * sum_all

  # --- Degrees of freedom ---
  dfA <- a - 1L
  dfR <- n - a
  dfT <- n - 1L

  # --- MS and pseudo-F ---
  MS_A <- SSA / dfA
  MS_R <- SSR / dfR
  F_A <- MS_A / MS_R

  anova_tbl <- data.frame(
    Source = c(nameA, "Residuals", "Total"),
    df = c(dfA, dfR, dfT),
    SS = c(SSA, SSR, SST),
    MS = c(MS_A, MS_R, NA_real_),
    pseudoF = c(F_A, NA_real_, NA_real_),
    check.names = FALSE
  )

  out_list <- list(
    SS = c(SSA = SSA, SSR = SSR, SST = SST),
    df = c(dfA = dfA, dfR = dfR, dfT = dfT),
    MS = c(MS_A = MS_A, MS_R = MS_R),
    F = c(F_A = F_A),
    SS_between = SSA,
    SS_total = SST,
    df_between = dfA,
    MS_residual = MS_R,
    pseudoF = F_A,
    R2 = SSA / SST,
    omega2 = (SSA - dfA * MS_R) / (SST + MS_R),
    method = method,
    transformation = transformation,
    model = model
  )

  switch(
    return,
    table = anova_tbl,
    list = out_list,
    both = list(table = anova_tbl, stats = out_list)
  )
}

#' Balanced sampling for effect size simulations
#'
#' Builds one balanced sample and computes ecological and inferential effect sizes.
#'
#' @param i Integer. Iteration index.
#' @param Y Integer vector/matrix with row indices to sample from.
#' @param mm Integer vector with number of groups in each iteration.
#' @param nn Integer vector with total sample size in each iteration.
#' @param YPU Integer vector labeling PSUs for balanced sampling.
#' @param HaSim 3D array of simulated Ha communities.
#' @param resultsHa Helper matrix with labels (at least dat.sim, k, m, n).
#' @param transformation Character. Data transformation.
#' @param method Character. Dissimilarity metric.
#'
#' @return A named list with one-row metrics for the selected simulation/sample.
#'
#' @keywords internal
#' @noRd
balanced_sampling_es <- function(
  i,
  Y,
  mm,
  nn,
  YPU,
  HaSim,
  resultsHa,
  transformation,
  method
) {
  m <- mm[i]
  n <- nn[i]

  df_idx <- tibble::tibble(
    idx = seq_along(Y),
    PU = YPU
  )

  psu_sel <- df_idx |>
    dplyr::distinct(PU) |>
    dplyr::slice_sample(n = m) |>
    dplyr::pull(PU)

  sel_rows <- df_idx |>
    dplyr::filter(PU %in% psu_sel) |>
    dplyr::group_by(PU) |>
    dplyr::slice_sample(n = n / m) |>
    dplyr::ungroup() |>
    dplyr::pull(idx)

  sel <- matrix(0L, nrow = nrow(df_idx), ncol = 1)
  sel[sel_rows, 1] <- 1L
  ones <- which(sel[, 1] %in% 1)

  ya <- HaSim[ones, , resultsHa[i, 1]]
  yHa <- dim(ya)[2] - 2

  comm <- ya[, 1:yHa, drop = FALSE]
  grp <- ya[, yHa + 2]

  infer_out <- dbmanova_oneway(
    x = comm,
    factEnv = grp,
    transformation = transformation,
    method = method,
    return = "table"
  )

  rs <- rowSums(comm)

  if (any(rs == 0, na.rm = TRUE)) {
    warning(
      "If you see this message several times, it is suggested you change dummy to TRUE",
      call. = FALSE
    )

    levs <- unique(grp)
    n_groups <- length(levs)

    eco <- list(
      ecological_effect = NA_real_,
      centroid_dist_matrix = matrix(
        NA_real_,
        nrow = n_groups,
        ncol = n_groups,
        dimnames = list(as.character(levs), as.character(levs))
      ),
      n_groups = n_groups,
      pcoa_points = data.frame(
        sample_id = seq_len(nrow(comm)),
        group = as.character(grp),
        Axis1 = NA_real_,
        Axis2 = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  } else {
    eco <- calc_dist(
      datHa = cbind(group = grp, comm),
      method = method,
      return = "both"
    )
  }

  # pseudoF = infer_out$table[1, "pseudoF"]
  R2 = infer_out[1, "SS"] / infer_out[3, "SS"]
  omega2 = (infer_out[1, "SS"] - infer_out[1, "df"] * infer_out[2, "MS"]) /
    (infer_out[3, "SS"] + infer_out[2, "MS"])

  pcoa_points <- eco$pcoa_points
  pcoa_points$dat_sim <- resultsHa[i, "dat.sim"]
  pcoa_points$k <- resultsHa[i, "k"]
  pcoa_points$m <- resultsHa[i, "m"]
  pcoa_points$n <- resultsHa[i, "n"]

  list(
    pseudoF = infer_out[1, "pseudoF"],
    omega2 = omega2,
    R2 = R2,
    SS_between = infer_out[1, "SS"],
    SS_total = infer_out[3, "SS"],
    df_between = infer_out[1, "df"],
    MS_residual = infer_out[2, "MS"],
    ecological_effect = eco$ecological_effect,
    centroid_dist_matrix = eco$centroid_dist_matrix,
    n_groups = eco$n_groups,
    infer_table = infer_out,
    pcoa_points = pcoa_points
  )
}

## Balanced sampling

#' Balanced sampling
#'
#' Develops the experimental design based on the provided conditions
#'
#' @param i pointer to the index in the list of experimental designs to try.
#' @param Y index to the data.frame the function will work with.
#' @param mm number of site the function is working with in each iteration.
#' @param nn number of samples to consider in each iteration.
#' @param YPU label for the sites in each iteration, as used by
#' [sampling::balancedtwostage()]
#' @param H0Sim simulated community from \code{SSP::simdata()} in which H0 is
#' true.
#' @param HaSim simulated community from \code{SSP::simdata()} in which H0 is
#' false.
#' @param resultsHa helper matrix that stores labels and later the results.
#' @param method appropriate distance/dissimilarity metric (e.g. Gower,
#' Bray–Curtis, Jaccard, etc).
#' @param transformation Mathematical function to reduce the weight of very
#' dominant species.
#'
#' @return a data frame with values for observed F (for H0 and Ha), and the Ha mean
#' squares for residuals and variation among sites.
#'
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @references Underwood, A. J. (1997). Experiments in ecology: their logical
#' design and interpretation using analysis of variance. Cambridge university
#' press.
#'
#' @importFrom tibble tibble
#' @importFrom dplyr filter distinct slice_sample pull group_by ungroup arrange
#'
#' @keywords internal
#' @noRd
#'

balanced_sampling <- function(
  i,
  Y,
  mm,
  nn,
  YPU,
  H0Sim,
  HaSim,
  resultsHa,
  transformation,
  method
) {
  # m and n equivalent to mm[i] y nn[i]
  m <- mm[i]
  n <- nn[i]

  # Constructing a light df with index and id for PSU
  df_idx <- tibble::tibble(
    idx = seq_along(Y),
    PU = YPU
  )

  # 1) Select PSUs (distinct m)
  psu_sel <- df_idx |>
    dplyr::distinct(PU) |>
    dplyr::slice_sample(n = m) |>
    dplyr::pull(PU)

  # 2) Select SSU within each PSU (n per PSU)
  sel_rows <- df_idx |>
    dplyr::filter(PU %in% psu_sel) |>
    dplyr::group_by(PU) |>
    dplyr::slice_sample(n = n / m) |>
    dplyr::ungroup() |>
    dplyr::pull(idx)

  # 3) Build a matrix with 0/1 to use as filter pointer
  sel <- matrix(0L, nrow = nrow(df_idx), ncol = 1)
  sel[sel_rows, 1] <- 1L

  ones <- which(sel[, 1] %in% 1)
  y0 <- H0Sim[ones, , resultsHa[i, 1]]
  ya <- HaSim[ones, , resultsHa[i, 1]]

  # Validate dimensions of H0 and Ha matrices
  if (!all(dim(y0) == dim(ya))) {
    stop("The dimensions of H0 and Ha data matrices do not match.")
  }

  yHa <- dim(y0)[2] - 2

  # Apply PERMANOVA to get pseudoF and mean squares
  result1 <- dbmanova_oneway(
    x = y0[, 1:yHa],
    factEnv = y0[, yHa + 2],
    transformation = transformation,
    method = method
  )
  result2 <- dbmanova_oneway(
    x = ya[, 1:yHa],
    factEnv = y0[, (yHa + 2)],
    transformation = transformation,
    method = method
  )

  # Create result matrix
  result0 <- matrix(nrow = 1, ncol = 7)
  colnames(result0) <- c("FobsH0", "FobsHa", "MSA", "MSR", "SSf", "SSr", "SSt")

  # Gather the results and return
  result0[, 1] <- result1[1, 5]
  result0[, 2] <- result2[1, 5]
  result0[, 3] <- result2[1, 4]
  result0[, 4] <- result2[2, 4]
  result0[, 5] <- result2[1, 3]
  result0[, 6] <- result2[2, 3]
  result0[, 7] <- result2[3, 3]
  return(result0)
}

#' dbMANOVA for nested design B(A) (2-way PERMANOVA)
#'
#' @param x        Matrix/data.frame with the community data (observations x species).
#' @param factEnv  data.frame with environmental factors in the first two columns:
#'                 - A: main factor (e.g. sector)
#'                 - B: nested factor (e.g. site)
#' @param transformation  Transformation process to reduce dominance. Options include: "none", "square root", "fourth root", "Log (X+1)".
#' @param method   Distance or dissimilarity method for using with vegan::vegdist. Defaults to "bray".
#' @param model    Label for the model: "nested.symmetric"
#' @param return   "table", "list" o "both". What format does the user require for the output? (ANOVA-like table, list of SS/df, or both).
#'
#' @return data.frame or list with SS, df, MS and pseudo-F:
#'         F_A = MS_A / MS_B(A), F_B(A) = MS_B(A) / MS_R.
#' @importFrom stats dist rchisq
#' @importFrom vegan vegdist
#'
#' @keywords internal
#' @noRd
#'

## PERMANOVA Two factors ----

dbmanova_nested <- function(
  x,
  factEnv,
  transformation = c("square root", "none", "fourth root", "Log (X+1)"),
  method = "bray",
  model = c("nested.symmetric"),
  return = c("table", "list", "both")
) {
  # --- Validating inputs ---
  transformation <- match.arg(transformation)
  method <- match.arg(method)
  model <- match.arg(model)
  return <- match.arg(return)

  if (!is.data.frame(factEnv)) {
    factEnv <- as.data.frame(factEnv)
  }
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("`x` must be matrix/data.frame (comunidad).")
  }

  Xcomm <- as.matrix(x)
  if (any(!is.finite(Xcomm))) {
    stop("`x` contains non-finite values (NA/NaN/Inf).")
  }

  if (nrow(Xcomm) != nrow(factEnv)) {
    stop("`x` y `factEnv` must have the same number of rows (samples).")
  }

  # Coerce factors
  facA <- factor(factEnv[[1]], exclude = NULL)
  facB <- factor(factEnv[[2]], exclude = NULL)
  facBA <- interaction(facA, facB, drop = TRUE)

  # --- Transformations to reduce dominance ---
  transform_comm <- function(M, how) {
    switch(
      how,
      "none" = M,
      "square root" = sqrt(M),
      "fourth root" = M^(1 / 4),
      "Log (X+1)" = log1p(M)
    )
  }

  Xtr <- transform_comm(Xcomm, transformation)

  # --- Distance and PCoA (Gower) ---
  d <- vegan::vegdist(Xtr, method = method)
  Dm <- as.matrix(d)

  A <- -0.5 * (Dm * Dm)

  n <- nrow(A)
  rbar <- rowMeans(A)
  cbar <- colMeans(A)
  abar <- mean(A)
  # G_ij = A_ij - rbar_i - cbar_j + abar
  G <- A
  G <- sweep(G, 1, rbar, `-`)
  G <- sweep(G, 2, cbar, `-`)
  G <- G + abar

  # --- basic statistics ---
  trG <- sum(diag(G))
  sum_all <- sum(G)

  # --- Traces for P_A G and P_BA G done with sums per block ---
  split_A <- split(seq_len(n), facA)
  split_BA <- split(seq_len(n), facBA)

  # trace(P_A G) = sum_g (1/n_g) * sum_{i,j in g} G_ij
  ng_A <- vapply(split_A, length, integer(1))
  sum_in_A <- vapply(
    split_A,
    function(idx) sum(G[idx, idx, drop = FALSE]),
    numeric(1)
  )
  trace_PAG <- sum(sum_in_A / ng_A)

  # trace(P_BA G) = sum_h (1/n_h) * sum_{i,j in h} G_ij
  ng_BA <- vapply(split_BA, length, integer(1))
  sum_in_BA <- vapply(
    split_BA,
    function(idx) sum(G[idx, idx, drop = FALSE]),
    numeric(1)
  )
  trace_PBAG <- sum(sum_in_BA / ng_BA)

  # --- SS y DF ---
  # SST and SSA use (1/n) * sum_all = trace(P_T G)
  SSA <- trace_PAG - (1 / n) * sum_all
  SSBA <- trace_PBAG - trace_PAG
  SSR <- trG - trace_PBAG
  SST <- trG - (1 / n) * sum_all

  a <- nlevels(facA)
  b_tot <- nlevels(facBA)

  # --- Degrees of freedom ---
  dfA <- a - 1L
  dfBA <- b_tot - a
  dfR <- n - b_tot
  dfT <- n - 1L

  # --- Mean squares and pseudo-F ---
  MS_A <- SSA / dfA
  MS_BA <- SSBA / dfBA
  MS_R <- SSR / dfR

  F_A <- MS_A / MS_BA
  F_BA <- MS_BA / MS_R

  anova_tbl <- data.frame(
    Source = c(
      names(factEnv)[1],
      paste0(names(factEnv)[2], "(", names(factEnv)[1], ")"),
      "Residuals",
      "Total"
    ),
    df = c(dfA, dfBA, dfR, dfT),
    SS = c(SSA, SSBA, SSR, SST),
    MS = c(MS_A, MS_BA, MS_R, NA_real_),
    `pseudoF` = c(F_A, F_BA, NA_real_, NA_real_),
    check.names = FALSE
  )

  out_list <- list(
    SS = c(SSA = SSA, SSBA = SSBA, SSR = SSR, SST = SST),
    df = c(dfA = dfA, dfBA = dfBA, dfR = dfR, dfT = dfT),
    MS = c(MS_A = MS_A, MS_BA = MS_BA, MS_R = MS_R),
    F = c(F_A = F_A, F_BA = F_BA),
    method = method,
    transformation = transformation,
    model = model
  )

  switch(
    return,
    table = anova_tbl,
    list = out_list,
    both = list(table = anova_tbl, stats = out_list)
  )
}


#' balanced_sampling_es_nested
#'
#' @returns a list
#'
#' @keywords internal
#' @noRd

balanced_sampling_es_nested <- function(
  i,
  Y1,
  mn,
  nSect,
  M,
  N,
  HaSim,
  resultsHa,
  factEnv,
  transformation,
  method,
  model
) {
  # Determine index for sampling units
  m_psu <- mn[i, 1] # número de PSUs a seleccionar
  n_ssu <- mn[i, 2] / m_psu # número de SSUs por PSU seleccionada

  Y1_df <- as.data.frame(Y1)
  row_id <- seq_len(nrow(Y1_df))
  PU_vec <- Y1_df[[2]] # equivalente a 'PU = Y1[, 2]' en balancedtwostage()

  # 1) Selección de PSUs
  psu_sel <- tibble::tibble(PU = PU_vec) |>
    dplyr::distinct(PU) |>
    dplyr::slice_sample(n = m_psu) |>
    dplyr::pull(PU)

  # 2) Selección de SSUs dentro de cada PSU seleccionada
  indice <- tibble::tibble(row_id = row_id, PU = PU_vec) |>
    dplyr::filter(PU %in% psu_sel) |>
    dplyr::group_by(PU) |>
    dplyr::slice_sample(n = n_ssu) |>
    dplyr::ungroup() |>
    dplyr::arrange(row_id) |>
    dplyr::pull(row_id)

  ones_n <- rep(indice, nSect)
  ones_s <- rep(c(0:(nSect - 1)) * M * N, each = length(indice))
  ones <- ones_n + ones_s

  rm(ones_n, ones_s, indice)

  # Extract samples from the datasets and evaluate with PERMANOVA
  ya <- HaSim[ones, , resultsHa[i, 1]]
  rownames(ya) <- ones
  factEnvX <- factEnv[ones, ]

  infer_out <- dbmanova_nested(
    x = ya,
    factEnv = factEnvX,
    transformation = transformation,
    method = method,
    model = model,
    return = "table"
  )

  rs <- rowSums(ya)

  if (any(rs == 0, na.rm = TRUE)) {
    warning(
      "If you see this message several times, it is suggested you change dummy to TRUE",
      call. = FALSE
    )

    levs_A <- unique(as.character(factEnvX[, 1]))
    levs_BA <- unique(interaction(factEnvX[, 1], factEnvX[, 2], drop = TRUE))

    n_groups_A <- length(levs_A)
    n_groups_BA <- length(levs_BA)

    eco <- list(
      ecological_effect_A = NA_real_,
      ecological_effect_BA = NA_real_,
      centroid_dist_matrix_A = matrix(
        NA_real_,
        nrow = n_groups_A,
        ncol = n_groups_A,
        dimnames = list(levs_A, levs_A)
      ),
      n_groups_A = n_groups_A,
      n_groups_BA = n_groups_BA,
      pcoa_points = data.frame(
        sample_id = seq_len(nrow(ya)),
        group = as.character(factEnvX[, 1]),
        sector = as.character(factEnvX[, 1]),
        site = as.character(factEnvX[, 2]),
        sector_site = as.character(interaction(
          factEnvX[, 1],
          factEnvX[, 2],
          drop = TRUE
        )),
        Axis1 = NA_real_,
        Axis2 = NA_real_,
        stringsAsFactors = FALSE
      ),
      site_sector_table = data.frame(
        sector = as.character(factEnvX[, 1]),
        site = as.character(factEnvX[, 2]),
        sector_site = as.character(interaction(
          factEnvX[, 1],
          factEnvX[, 2],
          drop = TRUE
        )),
        n_site = NA_real_,
        dist_to_sector = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  } else {
    eco <- calc_dist_nested(
      x = ya,
      factEnv = factEnvX,
      method = method,
      return = "both"
    )
  }

  # pseudoF = infer_out$table[1, "pseudoF"]
  R2_A = infer_out[1, "SS"] / infer_out[4, "SS"]
  R2_BA = infer_out[2, "SS"] / infer_out[4, "SS"]
  omega2_A = (infer_out[1, "SS"] - infer_out[1, "df"] * infer_out[2, "MS"]) /
    (infer_out[4, "SS"] + infer_out[2, "MS"])
  omega2_BA = (infer_out[2, "SS"] - infer_out[2, "df"] * infer_out[3, "MS"]) /
    (infer_out[4, "SS"] + infer_out[3, "MS"])

  pcoa_points <- eco$pcoa_points
  pcoa_points$dat_sim <- resultsHa[i, "dat.sim"]
  pcoa_points$k <- resultsHa[i, "k"]
  pcoa_points$m <- resultsHa[i, "m"]
  pcoa_points$n <- resultsHa[i, "n"]

  list(
    pseudoF_A = infer_out[1, "pseudoF"],
    pseudoF_BA = infer_out[2, "pseudoF"],
    omega2_A = omega2_A,
    omega2_BA = omega2_BA,
    R2_A = R2_A,
    R2_BA = R2_BA,
    SS_A = infer_out[1, "SS"],
    SS_BA = infer_out[2, "SS"],
    SS_total = infer_out[4, "SS"],
    df_A = infer_out[1, "df"],
    df_BA = infer_out[2, "df"],
    MS_BA = infer_out[2, "MS"],
    MS_residual = infer_out[3, "MS"],
    ecological_effect_A = eco$ecological_effect_A,
    ecological_effect_BA = eco$ecological_effect_BA,
    centroid_dist_matrix_A = eco$centroid_dist_matrix_A,
    n_groups_A = eco$n_groups_A,
    n_groups_BA = eco$n_groups_BA,
    infer_table = infer_out,
    pcoa_points = pcoa_points,
    site_sector_table = eco$site_sector_table
  )
}


#' Balanced sampling 2
#'
#' Develops the experimental design based on the provided conditions
#'
#' @param i pointer to the index in the list of experimental designs to try.
#' @param NN Total number of iterations that the experiment will consider.
#' @param Y1 A data frame with two columns, one indicates the auxiliary variables
#' on which the sample must be balanced and the other contains the vector of integers
#' that defines the primary sampling units. This is used by
#' \code{sampling::balancedtwostage()}
#' @param mn A data frame with two columns, one indicates the number of primary
#' sampling units to be selected and the other the number of second-stage sampling
#' units to be selected in the iteration. This is used by
#' \code{sampling::balancedtwostage()}
#' @param nSect Total number of sectors to be simulated in each data set.
#' @param M Total number of replicates to be simulated in each data set.
#' @param N Total number of samples to be simulated in each site.
#' @param H0Sim simulated community from \code{SSP::simdata()} in which H0 is true.
#' @param HaSim simulated community from \code{SSP::simdata()} in which H0 is false.
#' @param resultsHa helper matrix that stores labels and later the results.
#' @param factEnv a data frame for indicating the treatment, replicate and sampling
#' unit lables in each experiment.
#' @param transformation Mathematical function to reduce the weight of very
#' dominant species.
#' @param method appropriate distance/dissimilarity metric (e.g. Gower,
#' Bray–Curtis, Jaccard, etc).
#' @param model which algorithm to use for the calculation? At the moment, the only
#' option is "nested.symmetric".
#'
#' @return a data frame with values for observed F (for H0 and Ha), and the Ha mean
#' squares for residuals (MS_R) and variation among sites (MS_B(A)).
#'
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @references Underwood, A. J. (1997). Experiments in ecology: their logical
#' design and interpretation using analysis of variance. Cambridge university
#' press.
#'
#' @importFrom tibble tibble
#' @importFrom dplyr filter distinct slice_sample pull group_by ungroup arrange
#'
#' @keywords internal
#' @noRd
#'

balanced_sampling2 <- function(
  i,
  NN,
  Y1,
  mn,
  nSect,
  M,
  N,
  H0Sim,
  HaSim,
  resultsHa,
  factEnv,
  transformation,
  method,
  model
) {
  # Determine index for sampling units
  m_psu <- mn[i, 1] # número de PSUs a seleccionar
  n_ssu <- mn[i, 2] / m_psu # número de SSUs por PSU seleccionada

  Y1_df <- as.data.frame(Y1)
  row_id <- seq_len(nrow(Y1_df))
  PU_vec <- Y1_df[[2]] # equivalente a 'PU = Y1[, 2]' en balancedtwostage()

  # 1) Selección de PSUs
  psu_sel <- tibble::tibble(PU = PU_vec) |>
    dplyr::distinct(PU) |>
    dplyr::slice_sample(n = m_psu) |>
    dplyr::pull(PU)

  # 2) Selección de SSUs dentro de cada PSU seleccionada
  indice <- tibble::tibble(row_id = row_id, PU = PU_vec) |>
    dplyr::filter(PU %in% psu_sel) |>
    dplyr::group_by(PU) |>
    dplyr::slice_sample(n = n_ssu) |>
    dplyr::ungroup() |>
    dplyr::arrange(row_id) |>
    dplyr::pull(row_id)

  ones_n <- rep(indice, nSect)
  ones_s <- rep(c(0:(nSect - 1)) * M * N, each = length(indice))
  ones <- ones_n + ones_s

  rm(ones_n, ones_s, indice)

  # Extract samples from the datasets and evaluate with PERMANOVA
  y0 <- H0Sim[ones, , resultsHa[i, 1]]
  rownames(y0) <- ones
  ya <- HaSim[ones, , resultsHa[i, 1]]
  rownames(ya) <- ones
  factEnvX <- factEnv[ones, ]

  result_0 <- dbmanova_nested(
    x = y0,
    factEnv = factEnvX,
    transformation = transformation,
    method = method,
    model = model
  )
  result_a <- dbmanova_nested(
    x = ya,
    factEnv = factEnvX,
    transformation = transformation,
    method = method,
    model = model
  )

  # Assemble results
  result1 <- c(result_0[1, 5], result_a[1, 5], result_a[2, 4], result_a[3, 4])

  return(result1)
}

## Other helper functions ====

#' mimimun permutations
#'
#' Determine if sampling effort allows for at least 100 permutations
#'
#' @param model which algorithm to use for the calculation? At the moment, the only
#' option is "nested.symmetric".
#' @param a Integer. Levels for the treatment factor.
#' @param m Integer. Levels for site within treatment. Only used in Nested Symmetric
#' experiments.
#' @param n Integer. Replicates in the experiment (either per treatment or site).
#' @param perm Integer. Minimum number of permutations needed to reject the null
#' hypothesis. Defaults to 100, as it would allow for rejecting with alpha = 0.05,
#' the user can change this value to make the testing more strict (e.g. 200 for
#' testing alpha = 0.01 or 5000 for testing alpha = 0.001).
#'
#' @return Logical. TRUE if the required number of permutations are guaranteed.
#'
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @keywords internal
#' @noRd
#'

minimum_cbo <- function(model, a, n, m = NULL, perm) {
  # thr <- c(log(100))
  thr <- perm

  if (model == "single.factor") {
    permA <- factorial(a * n) / (factorial(a) * (factorial(n)^a)) #num. permutaciones factor A
    # permA <- lgamma(a * n+1) - lgamma(a + 1) - a * lgamma(n + 1)

    return(permA > thr)
  } else {
    if (is.null(m)) {
      stop("m is required for the nested model")
    }

    permA <- factorial(a * m) / (factorial(a) * (factorial(m)^a)) #num. permutaciones factor A
    permBA <- factorial(n)^(a * m) #num permutaciones factor B(A)
    # permA <- lgamma(a * m + 1) - lgamma(a + 1) - a * lgamma(m + 1)
    # permBA <- a * (lgamma(m * n + 1) - lgamma(m + 1) - m * lgamma(n + 1))

    return((permA > thr) & (permBA > thr))
  }
}

#' PERMANOVA by exchanging labels
#'
#' Classic Permutational Multivariate Analysis of Variance
#'
#' @param data Data frame where columns represent species names and rows correspond
#' to samples.
#' @param factEnv Data frame containing the environmental factors for the data.
#' @param method Character. Dissimilarity metric used [vegan::vegdist()]. Common
#' options include: "Gower", "Bray–Curtis", "Jaccard", etc.
#' @param transformation Character. Transformation applied to reduce the weight
#' of dominant species: "square root", "fourth root", "Log (X+1)", "P/A", "none".
#' @param dummy Character. Transformation applied to reduce the weight
#' of dominant species: "square root", "fourth root", "Log (X+1)", "P/A", "none".
#' @param model Character. Select the model to use. Options, so far, are
#' @param k Integer. Number of resampling iterations. Defaults to 50.
#'
#' @return PERMANOVA table
#'
#' @importFrom vegan vegdist betadisper
#'
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @keywords internal
#' @noRd
#'

label_permanova <- function(
  dataP,
  factEnvP,
  method,
  transformation,
  dummy,
  model
) {
  if (dummy) {
    dataP$dummy = 1
  }

  # Apply transformation and calculate distance matrix
  if (transformation == "square root") {
    x.t <- sqrt(dataP)
    d <- vegan::vegdist(x.t, method = method)
  } else if (transformation == "fourth root") {
    x.t <- sqrt(sqrt(dataP))
    d <- vegan::vegdist(x.t, method = method)
  } else if (transformation == "Log (X+1)") {
    x.t <- log(dataP + 1)
    d <- vegan::vegdist(x.t, method = method)
  } else if (transformation == "P/A") {
    x.t <- 1 * (dataP > 0)
    d <- vegan::vegdist(x.t, method = method, binary = TRUE)
  } else {
    x.t <- dataP
    d <- vegan::vegdist(x.t, method = method)
  }
  rm(dataP)

  # Compute the number of permutations available to the experiment,
  # then compare it with the given k
  a = nlevels(as.factor(factEnvP[, 1])) # number of treatments (A)
  b = length(unique(as.factor(factEnvP[, 2]))) # number of replicates (B)
  factEnvP["secsit"] <- paste0(factEnvP[, 1], factEnvP[, 2]) # intersections AB
  nBA = nlevels(as.factor(factEnvP$secsit)) # number of intersections AB
  nRep = dim(factEnvP)[1] / nBA # number of times we're repeating each intersection
  nNm = unique(factEnvP[, 1]) # unique values for the sectors
  nScSt = unique(factEnvP$secsit) # unique values for the intersections site-sector

  permutaciones_rep <- replicate(
    999,
    expr = factEnvP[sample(nrow(factEnvP)), ],
    simplify = FALSE
  )

  # calculates SS for all
  SST <- SS(d)[2]

  permList <- vector("list", 999 + 1)
  # degrees of freedom
  DoFA <- a - 1
  DoFBA <- a * (b - 1)
  DoFR <- a * b * (nRep - 1)
  DoFT <- (a * b * nRep) - 1

  for (i in c(1:999)) {
    # dataframe in this iteration
    currentPerm <- permutaciones_rep[[i]]

    # calculates SS within replicates
    # secsite_groups <- split(rownames(permutaciones_rep[[i]]), factEnvP$secsit)
    secsite_groups <- split(rownames(x.t), currentPerm$secsit)
    listR <- sapply(
      secsite_groups,
      function(rw) {
        SS(vegan::vegdist(x.t[rw, ], method = method))
      },
      simplify = "array"
    )
    SSR <- sum(listR[2, ])
    # calculates SS_B(A)
    # sector_groups <- split(rownames(permutaciones_rep[[i]]), factEnvP[,1])
    sector_groups <- split(rownames(x.t), currentPerm$Locality)

    listBA <- sapply(
      sector_groups,
      function(rw) {
        dBA <- vegan::vegdist(x.t[rw, ], method = "bray")
        tCentroid <- vegan::betadisper(
          dBA,
          group = currentPerm[rw, "secsit"],
          type = "centroid",
          bias.adjust = FALSE
        )
        Eig <- which(tCentroid$eig > 0)
        SS(vegan::vegdist(
          tCentroid$centroids[, Eig, drop = FALSE],
          method = "euclidean"
        ))
      },
      simplify = "array"
    )
    SSBA <- sum(listBA[2, ]) * nRep # * nRep is added so that when calculating
    # the variation between nested groups it is weighted by the size of the subgroup
    # rm(sector_groups, listBA)

    # calculates SSA
    SSA <- SST - SSBA - SSR

    # fill the permanova table
    # mean squares
    MSA <- SSA / DoFA
    MSBA <- SSBA / DoFBA
    MSR <- SSR / DoFR

    # observed pseudoF
    FobsA <- MSA / MSBA
    FobsBA <- MSBA / MSR

    Fobs <- as.data.frame(matrix(nrow = 4, ncol = 4))
    colnames(Fobs) <- c("SS", "DoF", "MS", "F")
    rownames(Fobs) <- c("A", "B(A)", "R", "T")
    Fobs[1, 1] <- SSA
    Fobs[1, 2] <- DoFA
    Fobs[1, 3] <- MSA
    Fobs[1, 4] <- FobsA
    Fobs[2, 1] <- SSBA
    Fobs[2, 2] <- DoFBA
    Fobs[2, 3] <- MSBA
    Fobs[2, 4] <- FobsBA
    Fobs[3, 1] <- SSR
    Fobs[3, 2] <- DoFR
    Fobs[3, 3] <- MSR
    Fobs[4, 1] <- SST
    Fobs[4, 2] <- DoFT

    permList[[i]] <- Fobs
  }

  # Compute ANOVA for the original data
  secsite_groups <- split(rownames(factEnvP), factEnvP$secsit)
  listR <- sapply(
    secsite_groups,
    function(rw) {
      SS(vegan::vegdist(x.t[rw, ], method = method))
    },
    simplify = "array"
  )
  SSR <- sum(listR[2, ])
  # calculates SS_B(A)
  sector_groups <- split(rownames(factEnvP), factEnvP[, 1])

  listBA <- sapply(
    sector_groups,
    function(rw) {
      dBA <- vegan::vegdist(x.t[rw, ], method = "bray")
      tCentroid <- vegan::betadisper(
        dBA,
        group = factEnvP[rw, "secsit"],
        type = "centroid",
        bias.adjust = FALSE
      )
      Eig <- which(tCentroid$eig > 0)
      SS(vegan::vegdist(
        tCentroid$centroids[, Eig, drop = FALSE],
        method = "euclidean"
      ))
    },
    simplify = "array"
  )
  SSBA <- sum(listBA[2, ]) * nRep

  # calculates SSA
  SSA <- SST - SSBA - SSR

  # fill the permanova table
  # mean squares
  MSA <- SSA / DoFA
  MSBA <- SSBA / DoFBA
  MSR <- SSR / DoFR

  # observed pseudoF
  FobsA <- MSA / MSBA
  FobsBA <- MSBA / MSR

  Fobs <- as.data.frame(matrix(nrow = 4, ncol = 4))
  colnames(Fobs) <- c("SS", "DoF", "MS", "F")
  rownames(Fobs) <- c("A", "B(A)", "R", "T")
  Fobs[1, 1] <- SSA
  Fobs[1, 2] <- DoFA
  Fobs[1, 3] <- MSA
  Fobs[1, 4] <- FobsA
  Fobs[2, 1] <- SSBA
  Fobs[2, 2] <- DoFBA
  Fobs[2, 3] <- MSBA
  Fobs[2, 4] <- FobsBA
  Fobs[3, 1] <- SSR
  Fobs[3, 2] <- DoFR
  Fobs[3, 3] <- MSR
  Fobs[4, 1] <- SST
  Fobs[4, 2] <- DoFT

  # Only necessary for full PERMANOVA
  permList[[1000]] <- Fobs

  # If this were full PERMANOVA, it's necessary to change the return to permList
  return(permList)
}

#' Weighted average of Similarity Percentages
#'
#' simper() from 'vegan' is adapted so that the between-group dissimilarities are condensed
#' into a single number that represents the per-species/per-site contribution.
#'
#' @param datHa Data frame where columns represent species names and rows correspond
#' to samples.
#'
#' @return Weighted average of species contribution to average between-group dissimilarity.
#'
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @importFrom vegan simper
#' @importFrom dplyr arrange mutate
#' @importFrom stats aggregate
#'
#' @keywords internal
#' @noRd
#'

use_simper <- function(datHa) {
  # Se llama a simper para evaluar contribución de las diferentes especies
  weight_spp <- vegan::simper(
    comm = datHa[, c(2:(ncol(datHa)))],
    group = datHa[, 1],
    permutations = 99
  )

  # Se guardan datos de apoyo
  pares <- names(weight_spp)
  flat_weight <- vector("list", length = length(weight_spp))
  names(flat_weight) <- pares

  # Se calcula promedio de contribución ponderado por par de sitios
  for (i in pares) {
    flat_weight[[i]] <- data.frame(
      spp = weight_spp[[i]]["species"],
      average = weight_spp[[i]]["average"]
    ) |>
      dplyr::mutate(
        norm = average / sum(average)
      )
  }

  # Se identifica especes con mayor contribución promedio
  relevant_spp <- do.call(rbind, flat_weight) |>
    stats::aggregate(norm ~ species, FUN = mean) |>
    dplyr::arrange(desc(norm))

  return(relevant_spp)
}

#' Weighted average of Similarity Percentages for nested experiments
#'
#' simper() from 'vegan' is adapted so that the between-group dissimilarities are condensed
#' into a single number that represents the per-species/per-site contribution.
#'
#' @param datHa Data frame where columns represent species names and rows correspond
#' to samples.
#'
#' @return Weighted average of species contribution to average between-group dissimilarity.
#'
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @importFrom vegan simper
#' @importFrom dplyr arrange mutate
#' @importFrom stats aggregate
#'
#' @keywords internal
#' @noRd
#'

use_simper_nested <- function(datHa) {
  # Se extraen los datos de diferentes comunidades
  pooled <- dplyr::bind_rows(
    lapply(datHa, function(x) {
      as.data.frame(x, stringsAsFactors = FALSE)
    }),
    .id = ".case"
  )

  # Indexes for creating the comm and group dataframes
  meta_cols <- c(".case", "sector", "sites", "N")
  species_cols <- setdiff(names(pooled), meta_cols)

  # Force community matrix to numeric
  pooled[species_cols] <- lapply(
    pooled[species_cols],
    function(z) {
      as.numeric(as.character(z))
    }
  )
  pooled$sector <- as.factor(pooled$sector)
  pooled$sites <- as.factor(pooled$sites)

  # Collapse observations to site level
  site_level <- pooled |>
    dplyr::group_by(sector, sites) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(species_cols),
        ~ mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )

  # Remove species with zero total abundance after collapse
  keep_spp <- colSums(site_level[, species_cols, drop = FALSE], na.rm = TRUE) >
    0
  species_cols_use <- species_cols[keep_spp]

  if (length(species_cols_use) == 0) {
    stop("All community columns are zero after site-level collapse.")
  }

  if (nlevels(site_level$sector) < 2) {
    stop("`sector` must have at least two levels for SIMPER.")
  }

  # Run simper at sector level
  weight_spp <- vegan::simper(
    comm = site_level[, species_cols_use, drop = FALSE],
    group = site_level$sector,
    permutations = 99
  )

  # Se guardan datos de apoyo
  pares <- names(weight_spp)
  flat_weight <- vector("list", length = length(weight_spp))
  names(flat_weight) <- pares

  # Se calcula promedio de contribución ponderado por par de sitios
  for (i in pares) {
    spp_i <- weight_spp[[i]][["species"]]
    avg_i <- weight_spp[[i]][["average"]]

    denom_i <- sum(avg_i, na.rm = TRUE)

    norm_i <- if (is.finite(denom_i) && denom_i > 0) {
      avg_i / denom_i
    } else {
      rep(0, length(avg_i))
    }

    flat_weight[[i]] <- data.frame(
      species = spp_i,
      norm = norm_i
    )
  }

  # Se identifica especes con mayor contribución promedio
  relevant_spp <- dplyr::bind_rows(flat_weight) |>
    dplyr::group_by(species) |>
    dplyr::summarise(norm = mean(norm, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(norm))

  return(relevant_spp)
}

#' Calculate distances to estimate ecological effect size
#'
#' Auxiliary function that calculates distances matrix for a simulated community.
#'
#' @param datHa Data frame where columns represent species names and rows correspond
#' to samples.
#'
#' @return A distance matrix for the groups in the ecological simulated community.
#'
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @importFrom vegan vegdist wcmdscale
#' @importFrom dplyr arrange mutate
#' @importFrom stats aggregate dist
#'
#' @keywords internal
#' @noRd
#'

calc_dist <- function(datHa, method = "bray", return = c("matrix", "both")) {
  return <- match.arg(return)

  datHa_site <- as.factor(datHa[, 1])
  datHa_ <- datHa[, c(2:ncol(datHa)), drop = FALSE]

  DistBC <- vegan::vegdist(datHa_, method = method)
  ord <- vegan::wcmdscale(DistBC, eig = TRUE, add = "lingoes")

  eig <- ord$eig
  eig_pos <- which(eig > 0)
  n_axes <- length(eig_pos)

  # Analytical coordinates.
  # Important: keep numeric PCoA coordinates separate from grouping metadata.
  if (n_axes == 0) {
    coords_num <- data.frame(
      Axis1 = rep(0, nrow(datHa_))
    )
    axis_cols <- "Axis1"
  } else {
    coords_num <- as.data.frame(ord$points[, eig_pos, drop = FALSE])
    names(coords_num) <- paste0("Axis", seq_len(ncol(coords_num)))
    axis_cols <- names(coords_num)
  }

  coords <- coords_num
  coords$site <- datHa_site

  # 2-axis view used for ordination plotting.
  # Analytical ES still uses all positive axes.
  if (n_axes == 0) {
    pcoa_points <- data.frame(
      sample_id = seq_len(nrow(datHa_)),
      group = datHa_site,
      Axis1 = rep(0, nrow(datHa_)),
      Axis2 = rep(0, nrow(datHa_)),
      stringsAsFactors = FALSE
    )
  } else if (n_axes == 1) {
    pcoa_points <- data.frame(
      sample_id = seq_len(nrow(datHa_)),
      group = datHa_site,
      Axis1 = coords_num[[1]],
      Axis2 = rep(0, nrow(datHa_)),
      stringsAsFactors = FALSE
    )
  } else {
    pcoa_points <- data.frame(
      sample_id = seq_len(nrow(datHa_)),
      group = datHa_site,
      Axis1 = coords_num[[1]],
      Axis2 = coords_num[[2]],
      stringsAsFactors = FALSE
    )
  }

  pcoa_points$Axis1 <- as.numeric(pcoa_points$Axis1)
  pcoa_points$Axis2 <- as.numeric(pcoa_points$Axis2)

  centroides <- stats::aggregate(
    coords[, axis_cols, drop = FALSE],
    by = list(site = coords$site),
    FUN = mean
  )

  if (nrow(centroides) <= 1) {
    centroid_dist <- matrix(
      0,
      nrow = 1,
      ncol = 1,
      dimnames = list(
        as.character(centroides$site),
        as.character(centroides$site)
      )
    )
  } else {
    centroid_dist <- stats::dist(centroides[, -1, drop = FALSE]) |>
      as.matrix()
    rownames(centroid_dist) <- centroides$site
    colnames(centroid_dist) <- centroides$site
  }

  n_groups <- nrow(centroides)
  group_size <- as.numeric(table(datHa_site)[as.character(centroides$site)])
  centroid_mat <- as.matrix(centroides[, -1, drop = FALSE])

  ecological_effect <- if (n_groups <= 1) {
    0
  } else if (n_groups == 2) {
    as.numeric(centroid_dist[1, 2])
  } else {
    c_bar <- colMeans(centroid_mat)
    sq_dist <- rowSums(
      (centroid_mat -
         matrix(
           c_bar,
           nrow = nrow(centroid_mat),
           ncol = ncol(centroid_mat),
           byrow = TRUE
         ))^2
    )
    sqrt(sum(group_size * sq_dist) / sum(group_size))
  }

  if (return == "matrix") {
    return(centroid_dist)
  }

  list(
    ecological_effect = ecological_effect,
    centroid_dist_matrix = centroid_dist,
    centroids = centroides,
    n_groups = n_groups,
    positive_axes = eig_pos,
    pcoa_points = pcoa_points
  )
}

#' Calculate distances to estimate ecological effect size for nested experiments
#'
#' Auxiliary function that calculates distances matrix for simulated communities.
#'
#' @param datHa Data frame where columns represent species names and rows correspond
#' to samples.
#'
#' @return A distance matrix for the groups in the ecological simulated community.
#'
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @importFrom vegan vegdist wcmdscale
#' @importFrom dplyr arrange mutate
#' @importFrom stats aggregate dist
#'
#' @keywords internal
#' @noRd
#'

calc_dist_nested <- function(
    x,
    factEnv,
    method = "bray",
    return = c("matrix", "both")
) {
  return <- match.arg(return)

  facA <- as.factor(factEnv[, 1])
  facB <- as.factor(factEnv[, 2])
  facBA <- interaction(facA, facB, drop = TRUE)

  DistBC <- vegan::vegdist(x, method = method)
  ord <- vegan::wcmdscale(DistBC, eig = TRUE, add = "lingoes")

  eig <- ord$eig
  eig_pos <- which(eig > 0)
  n_axes <- length(eig_pos)

  # Analytical coordinates.
  # Keep PCoA axes separate from sector/site metadata.
  if (n_axes == 0) {
    coords_num <- data.frame(
      Axis1 = rep(0, nrow(x))
    )
    axis_cols <- "Axis1"
  } else {
    coords_num <- as.data.frame(ord$points[, eig_pos, drop = FALSE])
    names(coords_num) <- paste0("Axis", seq_len(ncol(coords_num)))
    axis_cols <- names(coords_num)
  }

  coords <- coords_num
  coords$sector <- facA
  coords$site <- facB
  coords$sector_site <- facBA

  # 2-axis view used for ordination plotting.
  # Analytical ES still uses all positive axes.
  if (n_axes == 0) {
    pcoa_points <- data.frame(
      sample_id = seq_len(nrow(x)),
      group = facA,
      sector = facA,
      site = facB,
      sector_site = facBA,
      Axis1 = rep(0, nrow(x)),
      Axis2 = rep(0, nrow(x)),
      stringsAsFactors = FALSE
    )
  } else if (n_axes == 1) {
    pcoa_points <- data.frame(
      sample_id = seq_len(nrow(x)),
      group = facA,
      sector = facA,
      site = facB,
      sector_site = facBA,
      Axis1 = coords_num[[1]],
      Axis2 = rep(0, nrow(x)),
      stringsAsFactors = FALSE
    )
  } else {
    pcoa_points <- data.frame(
      sample_id = seq_len(nrow(x)),
      group = facA,
      sector = facA,
      site = facB,
      sector_site = facBA,
      Axis1 = coords_num[[1]],
      Axis2 = coords_num[[2]],
      stringsAsFactors = FALSE
    )
  }

  pcoa_points$Axis1 <- as.numeric(pcoa_points$Axis1)
  pcoa_points$Axis2 <- as.numeric(pcoa_points$Axis2)

  # Centroids for A
  centroides_A <- stats::aggregate(
    coords[, axis_cols, drop = FALSE],
    by = list(sector = coords$sector),
    FUN = mean
  )

  if (nrow(centroides_A) <= 1) {
    centroid_dist_A <- matrix(
      0,
      nrow = 1,
      ncol = 1,
      dimnames = list(
        as.character(centroides_A$sector),
        as.character(centroides_A$sector)
      )
    )
  } else {
    centroid_dist_A <- stats::dist(centroides_A[, -1, drop = FALSE]) |>
      as.matrix()
    rownames(centroid_dist_A) <- centroides_A$sector
    colnames(centroid_dist_A) <- centroides_A$sector
  }

  n_groups_A <- nrow(centroides_A)
  group_size_A <- as.numeric(table(facA)[as.character(centroides_A$sector)])
  centroid_mat_A <- as.matrix(centroides_A[, -1, drop = FALSE])

  ecological_effect_A <- if (n_groups_A <= 1) {
    0
  } else if (n_groups_A == 2) {
    as.numeric(centroid_dist_A[1, 2])
  } else {
    c_bar <- colMeans(centroid_mat_A)
    sq_dist <- rowSums(
      (centroid_mat_A -
         matrix(
           c_bar,
           nrow = nrow(centroid_mat_A),
           ncol = ncol(centroid_mat_A),
           byrow = TRUE
         ))^2
    )
    sqrt(sum(group_size_A * sq_dist) / sum(group_size_A))
  }

  # Centroids for B(A)
  centroides_BA <- stats::aggregate(
    coords[, axis_cols, drop = FALSE],
    by = list(
      sector = coords$sector,
      site = coords$site,
      sector_site = coords$sector_site
    ),
    FUN = mean
  )

  site_size <- as.data.frame(table(coords$sector_site))
  names(site_size) <- c("sector_site", "n_site")
  site_size$sector_site <- as.character(site_size$sector_site)

  site_lookup <- unique(coords[, c("sector", "site", "sector_site")])
  site_lookup$sector_site <- as.character(site_lookup$sector_site)

  site_size <- dplyr::left_join(
    site_size,
    site_lookup,
    by = "sector_site"
  )

  # Many-to-one join:
  # each site centroid receives the centroid of its parent sector.
  centroides_BA_join <- dplyr::left_join(
    centroides_BA,
    centroides_A,
    by = "sector",
    suffix = c("_site", "_sector")
  )

  axis_cols_site <- paste0(axis_cols, "_site")
  axis_cols_sector <- paste0(axis_cols, "_sector")

  dist_to_sector <- sqrt(
    rowSums(
      (
        as.matrix(centroides_BA_join[, axis_cols_site, drop = FALSE]) -
          as.matrix(centroides_BA_join[, axis_cols_sector, drop = FALSE])
      )^2
    )
  )

  site_sector_table <- dplyr::left_join(
    centroides_BA_join[
      ,
      c("sector", "site", "sector_site", axis_cols_site, axis_cols_sector),
      drop = FALSE
    ],
    site_size,
    by = c("sector", "site", "sector_site")
  )

  site_sector_table$dist_to_sector <- as.numeric(dist_to_sector)

  ecological_effect_BA <- if (nrow(site_sector_table) == 0) {
    0
  } else {
    sqrt(
      sum(site_sector_table$n_site * site_sector_table$dist_to_sector^2) /
        sum(site_sector_table$n_site)
    )
  }

  if (return == "matrix") {
    return(centroid_dist_A)
  }

  list(
    ecological_effect_A = ecological_effect_A,
    ecological_effect_BA = ecological_effect_BA,
    centroid_dist_matrix_A = centroid_dist_A,
    centroids_A = centroides_A,
    centroids_BA = centroides_BA,
    n_groups_A = n_groups_A,
    n_groups_BA = nrow(centroides_BA),
    positive_axes = eig_pos,
    pcoa_points = pcoa_points,
    site_sector_table = site_sector_table
  )
}

#' pooling fw between sectors
#' @keywords internal
#' @noRd
#'

pool_fw_across_sectors <- function(ListParam0, sector_weights = NULL) {
  # Validate input
  if (!is.list(ListParam0) || length(ListParam0) == 0) {
    stop("`ListParam0` must be a non-empty list.")
  }

  if (is.null(names(ListParam0))) {
    stop("`ListParam0` must be a named list, with sector labels as names.")
  }

  # Default weights: equal weights across sectors
  if (is.null(sector_weights)) {
    sector_weights <- rep(1, length(ListParam0))
    names(sector_weights) <- names(ListParam0)
  }

  # Validate weights
  if (is.null(names(sector_weights))) {
    stop("`sector_weights` must be a named numeric vector.")
  }

  if (!all(names(ListParam0) %in% names(sector_weights))) {
    stop(
      "All sector names in `ListParam0` must be present in `sector_weights`."
    )
  }

  sector_weights <- sector_weights[names(ListParam0)]

  # 1. Bind all sector-specific parameter tables in long format
  param_long <- purrr::imap_dfr(ListParam0, function(df, sec) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)

    if (!all(c("Species", "fw") %in% names(df))) {
      stop(
        "Each element of `ListParam0` must contain columns `Species` and `fw`."
      )
    }

    df$sector <- sec
    df[, c("Species", "sector", "fw")]
  })

  # 2. Complete all Species x sector combinations
  all_species <- sort(unique(param_long$Species))
  all_sectors <- names(ListParam0)

  param_long <- tidyr::complete(
    param_long,
    Species = all_species,
    sector = all_sectors,
    fill = list(fw = 0)
  )

  # 3. Add sector weights
  weights_df <- data.frame(
    sector = names(sector_weights),
    weight = as.numeric(sector_weights),
    stringsAsFactors = FALSE
  )

  param_long <- dplyr::left_join(param_long, weights_df, by = "sector")

  # 4. Weighted average across sectors
  param_pool <- param_long |>
    dplyr::group_by(Species) |>
    dplyr::summarise(
      fw_pool = sum(fw * weight, na.rm = TRUE) / sum(weight, na.rm = TRUE),
      n_sectors_present = sum(fw > 0, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(Species)

  return(param_pool)
}

# helper para blindar balanced_sampling
#' @keywords internal
#' @noRd
#'

safe_balanced_sampling_es_nested <- function(
  i,
  Y1,
  mn,
  nSect,
  M,
  N,
  HaSim,
  resultsHa,
  factEnv,
  transformation,
  method,
  model
) {
  tryCatch(
    {
      out <- balanced_sampling_es_nested(
        i = i,
        Y1 = Y1,
        mn = mn,
        nSect = nSect,
        M = M,
        N = N,
        HaSim = HaSim,
        resultsHa = resultsHa,
        factEnv = factEnv,
        transformation = transformation,
        method = method,
        model = model
      )

      out$ok <- TRUE
      out$error_message <- NA_character_

      out
    },
    error = function(e) {
      empty_pcoa <- data.frame(
        sample_id = integer(0),
        group = character(0),
        sector = character(0),
        site = character(0),
        sector_site = character(0),
        Axis1 = numeric(0),
        Axis2 = numeric(0),
        stringsAsFactors = FALSE
      )

      empty_site_sector <- data.frame(
        sector = character(0),
        site = character(0),
        sector_site = character(0),
        n_site = numeric(0),
        dist_to_sector = numeric(0),
        stringsAsFactors = FALSE
      )

      list(
        pseudoF_A = NA_real_,
        pseudoF_BA = NA_real_,
        omega2_A = NA_real_,
        omega2_BA = NA_real_,
        R2_A = NA_real_,
        R2_BA = NA_real_,
        SS_A = NA_real_,
        SS_BA = NA_real_,
        SS_total = NA_real_,
        df_A = NA_real_,
        df_BA = NA_real_,
        MS_residual = NA_real_,
        ecological_effect_A = NA_real_,
        ecological_effect_BA = NA_real_,
        centroid_dist_matrix_A = matrix(NA_real_, nrow = 0, ncol = 0),
        n_groups_A = NA_real_,
        n_groups_BA = NA_real_,
        infer_table = NULL,
        pcoa_points = empty_pcoa,
        site_sector_table = empty_site_sector,
        ok = FALSE,
        error_message = conditionMessage(e)
      )
    }
  )
}

#' @keywords internal
#' @noRD
#'

rank_fw_contribution <- function(ListParamA, ParamPoolA, weights = NULL) {
  if (!is.list(ListParamA) || length(ListParamA) == 0) {
    stop("`ListParamA` must be a non-empty named list.", call. = FALSE)
  }

  if (is.null(names(ListParamA))) {
    stop("`ListParamA` must be named by group or sector.", call. = FALSE)
  }

  if (!all(c("Species", "fw_pool") %in% names(ParamPoolA))) {
    stop("`ParamPoolA` must contain columns `Species` and `fw_pool`.", call. = FALSE)
  }

  groups <- names(ListParamA)

  if (is.null(weights)) {
    weights <- rep(1, length(groups))
    names(weights) <- groups
  }

  if (is.null(names(weights))) {
    stop("`weights` must be a named numeric vector.", call. = FALSE)
  }

  weights <- weights[groups]

  param_A <- purrr::imap_dfr(ListParamA, function(obj, group_i) {
    par_i <- as.data.frame(obj$par)

    if (!all(c("Species", "fw") %in% names(par_i))) {
      stop(
        "Each element of `ListParamA` must contain `$par` with columns ",
        "`Species` and `fw`.",
        call. = FALSE
      )
    }

    data.frame(
      Species = par_i$Species,
      group = group_i,
      fw = par_i$fw,
      stringsAsFactors = FALSE
    )
  })

  all_species <- sort(unique(c(param_A$Species, ParamPoolA$Species)))

  param_A <- tidyr::complete(
    param_A,
    Species = all_species,
    group = groups,
    fill = list(fw = 0)
  )

  weights_df <- data.frame(
    group = names(weights),
    weight = as.numeric(weights),
    stringsAsFactors = FALSE
  )

  param_A <- dplyr::left_join(param_A, weights_df, by = "group")

  param_A <- dplyr::left_join(
    param_A,
    ParamPoolA[, c("Species", "fw_pool")],
    by = "Species"
  )

  out <- param_A |>
    dplyr::group_by(Species) |>
    dplyr::summarise(
      score = sqrt(
        sum(weight * (fw - fw_pool)^2, na.rm = TRUE) /
          sum(weight, na.rm = TRUE)
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      norm = if (sum(score, na.rm = TRUE) > 0) {
        score / sum(score, na.rm = TRUE)
      } else {
        0
      }
    ) |>
    dplyr::transmute(
      species = Species,
      norm = norm
    ) |>
    dplyr::arrange(dplyr::desc(norm))

  out
}
