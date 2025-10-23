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
  result0 <- matrix(nrow = 1, ncol = 4)
  colnames(result0) <- c("FobsH0", "FobsHa", "MSA", "MSR")

  # Gather the results and return
  result0[, 1] <- result1[1, 5]
  result0[, 2] <- result2[1, 5]
  result0[, 3] <- result2[1, 4]
  result0[, 4] <- result2[2, 4]
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
    stop("`x` y `factEnv` must have the same number of rows (muestras).")
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

## Other helper functions

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
