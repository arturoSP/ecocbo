#' Calculate Simulated Effect Sizes
#'
#' Descripcion
#'
#' @param data Data frame where columns represent species names and rows correspond
#' to samples.
#'   - For `"single.factor"` analysis: The first column should indicate the replicate
#'   to which the sample belongs.
#'   - For `"nested.symmetric"` analysis: The first column should indicate the
#'   treatment, and the second column should indicate the replicate.
#' @param steps Number of steps needed to calculate the relevance depletion.
#' @param type Character. Nature of the data to be processed. It may be presence
#'  / absence ("P/A"), counts of individuals ("counts"), or coverage ("cover").
#' @param Sest.method Character Method for estimating species richness using
#' [vegan::specpool()]. Available methods are the incidence-based Chao ("chao"),
#' first order jackknife ("jack1"), second order jackknife ("jack2") and Bootstrap
#' ("boot"). By default, the average ("average") of the four estimates is used.
#' @param cases Integer. Number of simulated datasets.
#' @param N Integer. Total number of samples simulated per site.
#' @param M Integer. Total number of replicates simulated per dataset. Not needed
#' for single factor experiments.
#' @param n Integer. Maximum number of samples to consider (must be `<= N`).
#' @param m Integer. Number of replicates to consider. (must be `<=M`). Not needed
#' for single factor experiments.
#' @param k Integer. Number of resampling iterations. Defaults to 50.
#' @param transformation Character. Transformation applied to reduce the weight
#' of dominant species: "square root", "fourth root", "Log (X+1)", "P/A", "none".
#' @param method Character. Dissimilarity metric used [vegan::vegdist()]. Common
#' options include: "Gower", "Bray–Curtis", "Jaccard", etc.
#' @param dummy Logical. If `TRUE`, adds a small constant to empty observations.
#' @param useParallel Logical.  If `TRUE`, enables parallel computation. Defaults
#' to `TRUE`.
#' @param model Character. Select the model to use. Options are `"single.factor"`
#' and `"nested.symmetric"`.
#' @param jitter.base Numeric. Standard deviation multiplier used to add Gaussian
#' jitter to \code{fs} and \code{fw}. Defaults to 0.5.
#'
#' @details
#' The input dataset should have:
#' - One or two leading columns for treatment/replicate labels.
#' - Subsequent columns representing species presence/absence, counts, or coverage.
#' - `"single.factor"` requires a single column for replicates.
#' - `"nested.symmetric"` requires two columns: treatment and replicate in that
#' order.
#'
#' @return \code{prep_data()} returns an object of class "ecocbo_data".
#'
#' An object of class "ecocbo_data" is a list containing:
#'   - \code{$Results}, a data frame that lists the estimates of pseudoF for
#'   \code{simH0} and \code{simHa}, useful for statistical power analysis. It also
#'   includes mean squares for variance component estimation.
#'   - \code{$model}, a label for keeping track of the model that is being used
#'   in the analysis.
#'   - \code{$a}, an integer for the number of treatments recorded from the original
#'   data.
#'
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @references
#' - Underwood, A. J. (1997). Experiments in ecology: their logical
#' design and interpretation using analysis of variance. Cambridge university
#' press.
#' - Underwood, A. J., & Chapman, M. G. (2003). Power, precaution,
#' Type II error and sampling design in assessment of environmental impacts.
#' Journal of Experimental Marine Biology and Ecology, 296(1), 49-70.
#'
#' @seealso
#' [prep_data()]
#' [sim_beta()]
#' [vegan::simper()]
#'
#' @aliases simes
#'
#' @export
#' @importFrom SSP assempar simdata
#'
#' @examples
#' \donttest{
#' ES1 <- sim_es(data = epiDat, steps = 10, type = "counts", Sest.method = "average",
#'                         cases = 5, N = 100, M = 10,
#'                         n = 5, m = 5, k = 30,
#'                         transformation = "none", method = "bray",
#'                         dummy = FALSE, useParallel = FALSE,
#'                         model = "single.factor",
#'                         jitter.base = 0)
#' }
#' ES1
#'

.nested_ecological_effect <- function(
  comm,
  main_factor,
  nested_factor,
  method = "bray"
) {
  DistBC <- vegan::vegdist(comm, method = method)
  ord <- vegan::wcmdscale(DistBC, eig = TRUE, add = "lingoes")

  eig <- ord$eig
  eig_pos <- which(eig > 0)
  coords <- as.data.frame(ord$points[, eig_pos, drop = FALSE])
  main_factor <- as.factor(main_factor)
  nested_factor <- as.factor(nested_factor)

  n_obs <- nrow(comm)
  if (ncol(coords) == 0) {
    pcoa_points <- data.frame(
      sample_id = seq_len(n_obs),
      group = main_factor,
      nested = nested_factor,
      Axis1 = rep(0, n_obs),
      Axis2 = rep(0, n_obs),
      stringsAsFactors = FALSE
    )
  } else if (ncol(coords) == 1) {
    pcoa_points <- data.frame(
      sample_id = seq_len(n_obs),
      group = main_factor,
      nested = nested_factor,
      Axis1 = coords[[1]],
      Axis2 = rep(0, n_obs),
      stringsAsFactors = FALSE
    )
  } else {
    pcoa_points <- data.frame(
      sample_id = seq_len(n_obs),
      group = main_factor,
      nested = nested_factor,
      Axis1 = coords[[1]],
      Axis2 = coords[[2]],
      stringsAsFactors = FALSE
    )
  }

  coords_nested <- cbind(
    data.frame(main = main_factor, nested = nested_factor),
    coords
  )

  centroids_nested <- stats::aggregate(. ~ main + nested, data = coords_nested, FUN = mean)
  centroids_main <- stats::aggregate(. ~ main, data = centroids_nested, FUN = mean)

  centroid_main_mat <- as.matrix(centroids_main[, -(1), drop = FALSE])
  dist_main_pairwise <- stats::dist(centroid_main_mat) |> as.matrix()
  rownames(dist_main_pairwise) <- centroids_main$main
  colnames(dist_main_pairwise) <- centroids_main$main

  A <- nrow(centroids_main)
  if (A >= 2) {
    upper <- dist_main_pairwise[upper.tri(dist_main_pairwise)]
    effect_ecol_main <- sqrt(mean(upper^2))
  } else {
    effect_ecol_main <- 0
  }

  centroids_with_main <- merge(
    centroids_nested,
    centroids_main,
    by = "main",
    suffixes = c("_nested", "_main")
  )

  axis_nested <- grep("_nested$", names(centroids_with_main), value = TRUE)
  axis_main <- sub("_nested$", "_main", axis_nested)

  d_to_main <- sqrt(rowSums(
    (as.matrix(centroids_with_main[, axis_nested, drop = FALSE]) -
      as.matrix(centroids_with_main[, axis_main, drop = FALSE]))^2
  ))

  disp_nested <- data.frame(
    main = centroids_with_main$main,
    nested = centroids_with_main$nested,
    dist_to_main = d_to_main,
    stringsAsFactors = FALSE
  )

  disp_nested_within_main <- stats::aggregate(
    dist_to_main ~ main,
    data = disp_nested,
    FUN = mean
  )

  list(
    ecological_effect = effect_ecol_main,
    effect_ecol_main = effect_ecol_main,
    centroids_nested = centroids_nested,
    centroids_main_factor = centroids_main,
    dist_main_pairwise = dist_main_pairwise,
    disp_nested_within_main = disp_nested_within_main,
    disp_nested_global = mean(disp_nested$dist_to_main),
    pcoa_points = pcoa_points,
    positive_axes = eig_pos
  )
}

.balanced_sampling_es_nested <- function(
  i,
  Y,
  mm,
  nn,
  YPU,
  HaSim,
  resultsHa,
  factEnv,
  n_species,
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
  ya_comm <- as.matrix(ya[, seq_len(n_species), drop = FALSE])
  fact_x <- factEnv[ones, , drop = FALSE]

  infer_out <- dbmanova_nested(
    x = ya_comm,
    factEnv = fact_x,
    transformation = transformation,
    method = method,
    model = "nested.symmetric",
    return = "table"
  )

  eco <- .nested_ecological_effect(
    comm = ya_comm,
    main_factor = fact_x[, 1],
    nested_factor = interaction(fact_x[, 1], fact_x[, 2], drop = TRUE),
    method = method
  )

  R2 <- infer_out[1, "SS"] / infer_out[4, "SS"]
  omega2 <- (infer_out[1, "SS"] - infer_out[1, "df"] * infer_out[2, "MS"]) /
    (infer_out[4, "SS"] + infer_out[2, "MS"])

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
    SS_total = infer_out[4, "SS"],
    df_between = infer_out[1, "df"],
    MS_residual = infer_out[3, "MS"],
    ecological_effect = eco$ecological_effect,
    effect_ecol_main = eco$effect_ecol_main,
    centroids_nested = eco$centroids_nested,
    centroids_main_factor = eco$centroids_main_factor,
    dist_main_pairwise = eco$dist_main_pairwise,
    disp_nested_within_main = eco$disp_nested_within_main,
    disp_nested_global = eco$disp_nested_global,
    n_groups = nrow(eco$centroids_main_factor),
    infer_table = infer_out,
    pcoa_points = pcoa_points
  )
}

sim_ES <- function(
  data,
  steps = 10,
  type = "counts",
  Sest.method = "average",
  cases = 5,
  N = 100,
  M = NULL,
  n,
  m = NULL,
  k = 50,
  transformation = "none",
  method = "bray",
  dummy = FALSE,
  useParallel = TRUE,
  model = "single.factor",
  jitter.base = 0.5
) {
  if (!model %in% c("single.factor", "nested.symmetric")) {
    stop(
      "`model` must be either `single.factor` or `nested.symmetric`."
    )
  }

  # read data and store it in two objects, one for H0 and one for Ha ----
  datH0 <- data
  datHa <- data

  if (model == "single.factor") {
    datH0[, 1] <- as.factor("zero")
    datHa[, 1] <- as.factor(data[, 1])
    a <- nlevels(datHa[, 1])
    M <- a
    n_species <- ncol(data) - 1
    factEnv <- NULL
  } else {
    datH0[, 1] <- as.factor("zero")
    datH0[, 2] <- as.factor(data[, 2])
    datHa[, 1] <- as.factor(data[, 1])
    datHa[, 2] <- as.factor(data[, 2])
    trt_nested <- unique(data.frame(
      main = datHa[, 1],
      nested = datHa[, 2]
    ))
    a <- nlevels(as.factor(datHa[, 1]))
    M <- nrow(trt_nested)
    n_species <- ncol(data) - 2
    factEnv <- trt_nested[rep(seq_len(nrow(trt_nested)), each = N), , drop = FALSE]
  }

  if (!is.null(m) && m > M) {
    stop("`m` must be <= total number of nested units available in simulation (`M`).")
  }

  ## Helper matrix to store labels ----
  NN <- cases * k * (n - 1)
  resultsHa <- matrix(NA_real_, nrow = NN, ncol = 6)
  colnames(resultsHa) <- c(
    "dat.sim",
    "k",
    "m",
    "n",
    "FobsHa",
    "FobsH0"
  )

  resultsHa[, 1] <- rep(seq(cases), times = 1, each = (k * (n - 1)))
  resultsHa[, 2] <- rep(1:k, times = (n - 1) * cases)
  resultsHa[, 3] <- if (model == "nested.symmetric" && !is.null(m)) m else M
  resultsHa[, 4] <- rep(seq(2, n), times = 1, each = k)

  Y <- cbind(1:(N * M))
  YPU <- as.numeric(gl(M, N))
  mm <- resultsHa[, 3]
  nn <- resultsHa[, 3] * resultsHa[, 4]

  # Calculate H0 simulation parameters
  parH0 <- SSP::assempar(data = datH0, type = type, Sest.method = Sest.method)
  parHa <- SSP::assempar(data = datHa, type = type, Sest.method = Sest.method)

  ## Calculate species similarity percentage from Ha ----
  simper_data <- if (model == "single.factor") {
    datHa
  } else {
    cbind(datHa[, 1, drop = FALSE], datHa[, -(1:2), drop = FALSE])
  }
  sppContribution <- use_simper(simper_data)

  ## Output container ----
  resultOut <- vector("list", length = NN * (steps + 1))
  pcoaOut <- vector("list", length = NN * (steps + 1))
  nestedDetailOut <- vector("list", length = NN * (steps + 1))

  simH0 <- SSP::simdata(
    parH0,
    cases = cases,
    N = N,
    sites = M,
    jitter.base = jitter.base
  )

  xH0 <- dim(simH0[[1]])[1]
  yH0 <- dim(simH0[[1]])[2]
  casesHa <- length(simH0)

  H0Sim <- simH0

  if (dummy == TRUE) {
    yH0 <- yH0 + 1
    for (i in seq_len(casesHa)) {
      H0Sim[[i]] <- cbind(simH0[[i]], dummy = 1)
      H0Sim[[i]] <- H0Sim[[i]][, c(1:(yH0 - 3), (yH0), (yH0 - 2):(yH0 - 1))]
    }
  }

  rm(simH0)
  H0Sim <- array(unlist(H0Sim), dim = c(xH0, yH0, casesHa))

  # Loop for species contribution attenuation ----
  n_spp <- dim(sppContribution)[1]
  propSpp <- 1 / steps

  for (st in 0:steps) {
    # Calculate number of species to alter
    removing <- st * propSpp
    adjustN <- ceiling(nrow(sppContribution) * removing)

    if (adjustN != 0) {
      adjusting <- sppContribution[1:adjustN, 1]

      # Obtain H0 parameters from null hypothesis
      zeroParameter <- parH0$par[
        parH0$par$Species %in% adjusting,
        c(1, 3)
      ]

      # Alter the Ha parameters with the corresponding null hypothesis version
      parHaTemp <- parHa
      parHaTemp$par[
        parHaTemp$par$Species %in% adjusting,
        c(3)
      ] <- zeroParameter[, 2]
    } else {
      parHaTemp <- parHa
    }

    simHa <- SSP::simdata(
      parHaTemp,
      cases = cases,
      N = N,
      sites = M,
      jitter.base = jitter.base
    )

    # Simulation arguments ---
    xH0 <- dim(simHa[[1]])[1]
    yH0 <- dim(simHa[[1]])[2]
    casesHa <- length(simHa)

    HaSim <- simHa

    if (dummy == TRUE) {
      yH0 <- yH0 + 1
      for (i in seq_len(casesHa)) {
        HaSim[[i]] <- cbind(simHa[[i]], dummy = 1)
        HaSim[[i]] <- HaSim[[i]][, c(1:(yH0 - 3), (yH0), (yH0 - 2):(yH0 - 1))]
      }
    }

    rm(simHa)

    HaSim <- array(unlist(HaSim), dim = c(xH0, yH0, casesHa))

    # Loop to calculate ecological and inferential effect sizes ----

    if (useParallel) {
      # Registering the cluster of workers with parabar
      parabar::configure_bar(type = "basic", style = 3)
      cl <- parabar::start_backend(
        cores = parallelly::availableCores() / 2,
        cluster_type = "psock",
        backend_type = "async"
      )

      # Exporing functions needed for the parallel iterations
      parabar::export(
        cl,
        variables = c(
          if (model == "single.factor") "balanced_sampling_es" else ".balanced_sampling_es_nested",
          "dbmanova_oneway",
          "dbmanova_nested",
          ".nested_ecological_effect",
          "calc_dist"
        ),
        environment = asNamespace("ecocbo")
      )

      # Executing the loop in parallel
      result1 <- if (model == "single.factor") {
        parabar::par_lapply(
          cl,
          x = 1:NN,
          fun = balanced_sampling_es,
          Y,
          mm,
          nn,
          YPU,
          HaSim,
          resultsHa,
          transformation,
          method
        )
      } else {
        parabar::par_lapply(
          cl,
          x = 1:NN,
          fun = .balanced_sampling_es_nested,
          Y,
          mm,
          nn,
          YPU,
          HaSim,
          resultsHa,
          factEnv,
          n_species,
          transformation,
          method
        )
      }

      parabar::stop_backend(cl)
    } else {
      pb <- txtProgressBar(max = NN, style = 3)
      result1 <- vector("list", length = NN)

      for (i in seq_len(NN)) {
        # Performs the operation iteratively in a for loop
        result1[[i]] <- if (model == "single.factor") {
          balanced_sampling_es(
            i,
            Y,
            mm,
            nn,
            YPU,
            HaSim,
            resultsHa,
            transformation,
            method
          )
        } else {
          .balanced_sampling_es_nested(
            i,
            Y,
            mm,
            nn,
            YPU,
            HaSim,
            resultsHa,
            factEnv,
            n_species,
            transformation,
            method
          )
        }

        # Updating the progress bar
        setTxtProgressBar(pb, i)
      }
      close(pb)
    }

    idx <- (st * NN + 1):((st + 1) * NN)

    resultOut[idx] <- lapply(seq_len(NN), function(i) {
      cur <- result1[[i]]

      out <- data.frame(
        reduction_level = st / steps,
        step = st,
        dat_sim = resultsHa[i, "dat.sim"],
        k = resultsHa[i, "k"],
        m = resultsHa[i, "m"],
        n = resultsHa[i, "n"],
        ecological_effect = cur$ecological_effect,
        omega2 = cur$omega2,
        R2 = cur$R2,
        pseudoF = cur$pseudoF,
        SS_between = cur$SS_between,
        SS_total = cur$SS_total,
        df_between = cur$df_between,
        MS_residual = cur$MS_residual,
        n_groups = cur$n_groups,
        stringsAsFactors = FALSE
      )

      if (model == "nested.symmetric") {
        out$effect_ecol_main <- cur$effect_ecol_main
        out$disp_nested_within_main <- cur$disp_nested_global
      }

      out
    })

    pcoaOut[idx] <- lapply(seq_len(NN), function(i) {
      cur <- result1[[i]]
      pts <- cur$pcoa_points
      pts$reduction_level <- st / steps
      pts$step <- st
      pts
    })

    if (model == "nested.symmetric") {
      nestedDetailOut[idx] <- lapply(seq_len(NN), function(i) {
        cur <- result1[[i]]
        list(
          step = st,
          reduction_level = st / steps,
          dat_sim = resultsHa[i, "dat.sim"],
          k = resultsHa[i, "k"],
          m = resultsHa[i, "m"],
          n = resultsHa[i, "n"],
          centroids_nested = cur$centroids_nested,
          centroids_main_factor = cur$centroids_main_factor,
          dist_main_pairwise = cur$dist_main_pairwise,
          disp_nested_within_main = cur$disp_nested_within_main
        )
      })
    }

    cat("Step ", st, ": Done! - ", steps - st, " remaining")
  }

  resultOut <- dplyr::bind_rows(resultOut)
  pcoaOut <- dplyr::bind_rows(pcoaOut)

  es_obj <- new_effect_size_data(
    data = resultOut,
    pcoa = pcoaOut,
    call = match.call()
  )

  if (model == "nested.symmetric") {
    attr(es_obj, "nested_ecology") <- nestedDetailOut
  }

  return(es_obj)
}

# =========================================================
# effect_size_data S3 class
# =========================================================

#' Create an effect_size_data object
#'
#' @param data A data frame with at least `reduction_level`
#'   and `ecological_effect`. Inferential columns such as
#'   `omega2`, `R2`, or `pseudoF` are optional but recommended.
#' @param centroid_distances Optional list or object storing
#'   centroid distance matrices.
#' @param pcoa Optional list or object storing PCoA results.
#' @param call Optional matched call.
#'
#' @return An object of class `effect_size_data`.
#' @export
new_effect_size_data <- function(
  data,
  centroid_distances = NULL,
  pcoa = NULL,
  call = NULL
) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  validate_effect_size_data(data)

  class(data) <- unique(c("effect_size_data", class(data)))
  attr(data, "centroid_distances") <- centroid_distances
  attr(data, "pcoa") <- pcoa
  attr(data, "call") <- call
  data
}

#' Validate an effect_size_data object
#'
#' @param x A data frame or `effect_size_data` object.
#'
#' @return Invisibly returns `TRUE` if validation passes.
#' @keywords internal
validate_effect_size_data <- function(x) {
  if (!is.data.frame(x)) {
    stop("`x` must be a data.frame or effect_size_data object.", call. = FALSE)
  }

  required_cols <- c("reduction_level", "ecological_effect")
  missing_cols <- setdiff(required_cols, names(x))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(x$ecological_effect)) {
    stop("`ecological_effect` must be numeric.", call. = FALSE)
  }

  inferential_candidates <- c("omega2", "R2", "pseudoF")
  inferential_present <- intersect(inferential_candidates, names(x))

  if (length(inferential_present) == 0) {
    warning(
      "No inferential columns found. Consider including one of: ",
      paste(inferential_candidates, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' @export
as.data.frame.effect_size_data <- function(x, ...) {
  class(x) <- setdiff(class(x), "effect_size_data")
  x
}

# =========================================================
# print() method
# =========================================================

#' Print method for effect_size_data objects
#'
#' @param x An object of class `effect_size_data`.
#' @param digits Number of digits to print.
#' @param preview_n Number of rows to preview.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.effect_size_data <- function(
  x,
  digits = max(3L, getOption("digits") - 2L),
  preview_n = 6L,
  ...
) {
  validate_effect_size_data(x)

  df <- as.data.frame(x)
  n_rows <- nrow(df)
  n_cols <- ncol(df)

  inferential_candidates <- c("omega2", "R2", "pseudoF")
  inferential_present <- intersect(inferential_candidates, names(df))

  reduction_n <- length(unique(df$reduction_level))
  reduction_range <- tryCatch(
    range(df$reduction_level, na.rm = TRUE),
    error = function(e) NULL
  )

  has_pcoa <- !is.null(attr(x, "pcoa"))

  cat("<effect_size_data>\n")
  cat("Rows:", n_rows, "\n")
  cat("Columns:", n_cols, "\n")
  cat("Reduction levels:", reduction_n, "\n")
  cat("PCoA data available:", if (has_pcoa) "yes" else "no", "\n")

  if (!is.null(reduction_range) && length(reduction_range) == 2) {
    cat(
      "Reduction range:",
      paste(signif(reduction_range, digits), collapse = " to "),
      "\n"
    )
  }

  if ("k" %in% names(df)) {
    cat("Communities (k):", length(unique(df$k)), "\n")
  }
  if ("m" %in% names(df)) {
    cat("Group labels (m):", length(unique(df$m)), "\n")
  }
  if ("n" %in% names(df)) {
    cat("Replicates (n):", length(unique(df$n)), "\n")
  }

  cat(
    "Inferential metrics available:",
    if (length(inferential_present) > 0) {
      paste(inferential_present, collapse = ", ")
    } else {
      "none"
    },
    "\n"
  )

  preview_cols <- intersect(
    c(
      "reduction_level",
      "k",
      "m",
      "n",
      "ecological_effect",
      "effect_ecol_main",
      "disp_nested_within_main",
      "omega2",
      "R2",
      "pseudoF"
    ),
    names(df)
  )

  if (length(preview_cols) > 0) {
    cat("\nPreview:\n")
    print(
      utils::head(df[, preview_cols, drop = FALSE], n = preview_n),
      digits = digits,
      row.names = FALSE
    )
  }

  invisible(x)
}

# =========================================================
# plot() helpers
# =========================================================

#' @keywords internal
.es_center <- function(z, summary_stat = c("median", "mean")) {
  summary_stat <- match.arg(summary_stat)
  if (summary_stat == "mean") {
    mean(z, na.rm = TRUE)
  } else {
    stats::median(z, na.rm = TRUE)
  }
}

#' @keywords internal
.es_low <- function(
  z,
  summary_stat = c("median", "mean"),
  interval = c("iqr", "sd")
) {
  summary_stat <- match.arg(summary_stat)
  interval <- match.arg(interval)

  if (interval == "sd") {
    center <- .es_center(
      z,
      summary_stat = if (summary_stat == "mean") "mean" else "median"
    )
    center - stats::sd(z, na.rm = TRUE)
  } else {
    stats::quantile(z, probs = 0.25, na.rm = TRUE, names = FALSE)
  }
}

#' @keywords internal
.es_high <- function(
  z,
  summary_stat = c("median", "mean"),
  interval = c("iqr", "sd")
) {
  summary_stat <- match.arg(summary_stat)
  interval <- match.arg(interval)

  if (interval == "sd") {
    center <- .es_center(
      z,
      summary_stat = if (summary_stat == "mean") "mean" else "median"
    )
    center + stats::sd(z, na.rm = TRUE)
  } else {
    stats::quantile(z, probs = 0.75, na.rm = TRUE, names = FALSE)
  }
}

#' @keywords internal
.summarise_effect_size_data <- function(
  x,
  inferential_var = c("omega2", "R2", "pseudoF"),
  summary_stat = c("median", "mean"),
  interval = c("iqr", "sd")
) {
  inferential_var <- match.arg(inferential_var)
  summary_stat <- match.arg(summary_stat)
  interval <- match.arg(interval)

  df <- as.data.frame(x)

  if (!inferential_var %in% names(df)) {
    stop(
      "Column `",
      inferential_var,
      "` is not available in the object.",
      call. = FALSE
    )
  }

  dplyr::group_by(df, reduction_level) |>
    dplyr::summarise(
      eco_center = .es_center(ecological_effect, summary_stat),
      eco_low = .es_low(ecological_effect, summary_stat, interval),
      eco_high = .es_high(ecological_effect, summary_stat, interval),
      inf_center = .es_center(.data[[inferential_var]], summary_stat),
      inf_low = .es_low(.data[[inferential_var]], summary_stat, interval),
      inf_high = .es_high(.data[[inferential_var]], summary_stat, interval),
      .groups = "drop"
    )
}

#' @keywords internal
.plot_effect_size_stacked <- function(
  x,
  inferential_var = c("omega2", "R2", "pseudoF"),
  summary_stat = c("median", "mean"),
  interval = c("iqr", "sd"),
  ribbon_alpha = 0.25,
  line_size = 0.4,
  point_size = 1.5
) {
  inferential_var <- match.arg(inferential_var)
  summary_stat <- match.arg(summary_stat)
  interval <- match.arg(interval)

  df_sum <- .summarise_effect_size_data(
    x = x,
    inferential_var = inferential_var,
    summary_stat = summary_stat,
    interval = interval
  )

  y_lab_inf <- switch(
    inferential_var,
    omega2 = expression(omega^2),
    R2 = expression(R^2),
    pseudoF = "pseudo-F"
  )

  base_theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(linewidth = 0.4),
      axis.ticks = ggplot2::element_line(linewidth = 0.2)
    )

  p_eco <- ggplot2::ggplot(df_sum, ggplot2::aes(x = reduction_level)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = eco_low, ymax = eco_high),
      alpha = ribbon_alpha
    ) +
    ggplot2::geom_line(ggplot2::aes(y = eco_center), linewidth = line_size) +
    ggplot2::geom_point(ggplot2::aes(y = eco_center), size = point_size) +
    ggplot2::labs(
      x = NULL,
      y = "Ecological effect size"
    ) +
    base_theme

  p_inf <- ggplot2::ggplot(df_sum, ggplot2::aes(x = reduction_level)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = inf_low, ymax = inf_high),
      alpha = ribbon_alpha
    ) +
    ggplot2::geom_line(ggplot2::aes(y = inf_center), linewidth = line_size) +
    ggplot2::geom_point(ggplot2::aes(y = inf_center), size = point_size) +
    ggplot2::labs(
      x = "% removed SIMPER contribution",
      y = y_lab_inf
    ) +
    base_theme

  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(patchwork::wrap_plots(p_eco, p_inf, ncol = 1))
  }

  message(
    "Package `patchwork` is not installed. Returning a faceted fallback plot."
  )

  df_long <- rbind(
    data.frame(
      reduction_level = df_sum$reduction_level,
      center = df_sum$eco_center,
      low = df_sum$eco_low,
      high = df_sum$eco_high,
      metric = "Ecological effect size"
    ),
    data.frame(
      reduction_level = df_sum$reduction_level,
      center = df_sum$inf_center,
      low = df_sum$inf_low,
      high = df_sum$inf_high,
      metric = inferential_var
    )
  )

  ggplot2::ggplot(df_long, ggplot2::aes(x = reduction_level, y = center)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = low, ymax = high),
      alpha = ribbon_alpha
    ) +
    ggplot2::geom_line(linewidth = line_size) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::facet_wrap(~metric, ncol = 1, scales = "free_y") +
    ggplot2::labs(
      x = "% removed SIMPER contribution",
      y = NULL
    ) +
    base_theme
}

#' @keywords internal
.plot_effect_size_scatter <- function(
  x,
  inferential_var = c("omega2", "R2", "pseudoF"),
  colour_by = c("reduction_level", "none"),
  scatter_data = c("raw", "summary"),
  summary_stat = c("median", "mean"),
  point_alpha = 0.7,
  point_size = 1.8,
  add_smooth = FALSE
) {
  inferential_var <- match.arg(inferential_var)
  colour_by <- match.arg(colour_by)
  scatter_data <- match.arg(scatter_data)
  summary_stat <- match.arg(summary_stat)

  df <- as.data.frame(x)

  if (!inferential_var %in% names(df)) {
    stop(
      "Column `",
      inferential_var,
      "` is not available in the object.",
      call. = FALSE
    )
  }

  if (scatter_data == "summary") {
    df <- dplyr::group_by(df, reduction_level) |>
      dplyr::summarise(
        ecological_effect = .es_center(ecological_effect, summary_stat),
        inferential_value = .es_center(.data[[inferential_var]], summary_stat),
        .groups = "drop"
      )
  } else {
    df$inferential_value <- df[[inferential_var]]
  }

  y_lab_inf <- switch(
    inferential_var,
    omega2 = expression(omega^2),
    R2 = expression(R^2),
    pseudoF = "pseudo-F"
  )

  mapping <- if (colour_by == "reduction_level") {
    ggplot2::aes(
      x = inferential_value,
      y = ecological_effect,
      colour = reduction_level
    )
  } else {
    ggplot2::aes(
      x = inferential_value,
      y = ecological_effect
    )
  }

  p <- ggplot2::ggplot(df, mapping) +
    ggplot2::geom_point(alpha = point_alpha, size = point_size) +
    ggplot2::labs(
      x = y_lab_inf,
      y = "Ecological effect size",
      colour = "% removed SIMPER contribution"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(linewidth = 0.4),
      axis.ticks = ggplot2::element_line(linewidth = 0.2)
    )

  if (add_smooth) {
    p <- p + ggplot2::geom_smooth(method = "lm", se = FALSE)
  }

  p
}

#' @keywords internal
.plot_effect_size_ordination <- function(
  x,
  reduction_levels = NULL,
  selection = c("median", "first", "min", "max"),
  show_centroids = TRUE,
  label_centroids = TRUE,
  point_alpha = 0.7,
  point_size = 1.8,
  centroid_size = 2.8,
  facet_ncol = NULL
) {
  selection <- match.arg(selection)

  pcoa_df <- .extract_pcoa_data(x)

  plot_df <- .select_representative_pcoa(
    x = x,
    pcoa_df = pcoa_df,
    reduction_levels = reduction_levels,
    selection = selection
  )

  centroids <- .compute_pcoa_centroids(plot_df)

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = Axis1, y = Axis2, colour = group)
  ) +
    ggplot2::geom_point(alpha = point_alpha, size = point_size)

  if (show_centroids) {
    p <- p +
      ggplot2::geom_point(
        data = centroids,
        ggplot2::aes(x = Axis1, y = Axis2, colour = group),
        inherit.aes = FALSE,
        size = centroid_size,
        shape = 4
      )
  }

  if (label_centroids) {
    p <- p +
      ggplot2::geom_text(
        data = centroids,
        ggplot2::aes(x = Axis1, y = Axis2, label = group, colour = group),
        inherit.aes = FALSE,
        vjust = -0.7,
        show.legend = FALSE
      )
  }

  p +
    ggplot2::facet_wrap(~panel_label, ncol = facet_ncol) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = "PCoA 1",
      y = "PCoA 2",
      colour = "Group"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(linewidth = 0.4),
      axis.ticks = ggplot2::element_line(linewidth = 0.2)
    )
}

# =========================================================
# plot() method
# =========================================================

#' Plot method for effect_size_data objects
#'
#' @param x An object of class `effect_size_data`.
#' @param type Type of plot. `"stacked"` returns the two-panel
#'   summary plot; `"scatter"` returns the inferential vs ecological plot.
#' @param inferential_var Which inferential metric to use.
#'   One of `"omega2"`, `"R2"`, or `"pseudoF"`.
#' @param summary_stat Summary statistic used in grouped summaries.
#'   One of `"median"` or `"mean"`.
#' @param interval Interval used in the stacked summary plot.
#'   One of `"iqr"` or `"sd"`.
#' @param colour_by Whether to colour scatter points by
#'   `reduction_level` or not.
#' @param scatter_data Whether scatter should show raw observations
#'   or summaries by `reduction_level`.
#' @param ribbon_alpha Alpha transparency for ribbons in stacked plots.
#' @param point_alpha Alpha transparency for points in scatter plots.
#' @param line_size Line width in stacked plots.
#' @param point_size Point size.
#' @param add_smooth Logical; add linear smooth to scatter plot.
#' @param ... Unused, reserved for future extensions.
#'
#' @return A ggplot object or a patchwork object.
#' @export
plot.effect_size_data <- function(
  x,
  type = c("stacked", "scatter", "ordination"),
  inferential_var = c("omega2", "R2", "pseudoF"),
  summary_stat = c("median", "mean"),
  interval = c("iqr", "sd"),
  colour_by = c("reduction_level", "none"),
  scatter_data = c("summary", "raw"),
  reduction_levels = NULL,
  selection = c("median", "first", "min", "max"),
  show_centroids = TRUE,
  label_centroids = TRUE,
  ribbon_alpha = 0.25,
  point_alpha = 0.7,
  line_size = 0.4,
  point_size = 1.8,
  centroid_size = 2.8,
  facet_ncol = NULL,
  add_smooth = FALSE,
  ...
) {
  validate_effect_size_data(x)

  type <- match.arg(type)
  inferential_var <- match.arg(inferential_var)
  summary_stat <- match.arg(summary_stat)
  interval <- match.arg(interval)
  colour_by <- match.arg(colour_by)
  scatter_data <- match.arg(scatter_data)
  selection <- match.arg(selection)

  if (type %in% c("stacked", "scatter") && !inferential_var %in% names(x)) {
    stop(
      "Column `",
      inferential_var,
      "` is not available in the object.",
      call. = FALSE
    )
  }

  switch(
    type,
    stacked = .plot_effect_size_stacked(
      x = x,
      inferential_var = inferential_var,
      summary_stat = summary_stat,
      interval = interval,
      ribbon_alpha = ribbon_alpha,
      line_size = line_size,
      point_size = point_size
    ),
    scatter = .plot_effect_size_scatter(
      x = x,
      inferential_var = inferential_var,
      colour_by = colour_by,
      scatter_data = scatter_data,
      summary_stat = summary_stat,
      point_alpha = point_alpha,
      point_size = point_size,
      add_smooth = add_smooth
    ),
    ordination = .plot_effect_size_ordination(
      x = x,
      reduction_levels = reduction_levels,
      selection = selection,
      show_centroids = show_centroids,
      label_centroids = label_centroids,
      point_alpha = point_alpha,
      point_size = point_size,
      centroid_size = centroid_size,
      facet_ncol = facet_ncol
    )
  )
}

# =========================================================
# summary() helpers
# =========================================================

#' @keywords internal
.es_quantile_name <- function(p) {
  paste0("q", formatC(100 * p, format = "fg", width = 1, digits = 0))
}

#' @keywords internal
.summarise_vector <- function(
  z,
  center = c("median", "mean"),
  probs = c(0.25, 0.5, 0.75)
) {
  center <- match.arg(center)

  z <- z[!is.na(z)]
  if (length(z) == 0) {
    out <- c(
      n = 0,
      mean = NA_real_,
      median = NA_real_,
      sd = NA_real_,
      min = NA_real_,
      max = NA_real_
    )
    qs <- stats::setNames(
      rep(NA_real_, length(probs)),
      .es_quantile_name(probs)
    )
    return(c(out, qs))
  }

  out <- c(
    n = length(z),
    mean = mean(z),
    median = stats::median(z),
    sd = stats::sd(z),
    min = min(z),
    max = max(z)
  )

  qs <- stats::quantile(z, probs = probs, na.rm = TRUE, names = FALSE)
  names(qs) <- .es_quantile_name(probs)

  c(out, qs)
}

#' @keywords internal
.summarise_by_reduction <- function(df, inferential_var, center, probs) {
  dplyr::group_by(df, reduction_level) |>
    dplyr::group_modify(function(.x, .y) {
      eco <- .summarise_vector(
        .x$ecological_effect,
        center = center,
        probs = probs
      )
      inf <- .summarise_vector(
        .x[[inferential_var]],
        center = center,
        probs = probs
      )

      out <- data.frame(
        n_obs = nrow(.x),
        eco_mean = eco["mean"],
        eco_median = eco["median"],
        eco_sd = eco["sd"],
        eco_min = eco["min"],
        eco_max = eco["max"],
        inf_mean = inf["mean"],
        inf_median = inf["median"],
        inf_sd = inf["sd"],
        inf_min = inf["min"],
        inf_max = inf["max"]
      )

      for (p in probs) {
        nm <- .es_quantile_name(p)
        out[[paste0("eco_", nm)]] <- eco[[nm]]
        out[[paste0("inf_", nm)]] <- inf[[nm]]
      }

      out
    }) |>
    dplyr::ungroup()
}

# =========================================================
# summary() method
# =========================================================

#' Summary method for effect_size_data objects
#'
#' @param object An object of class `effect_size_data`.
#' @param inferential_var Which inferential metric to summarise.
#'   One of `"omega2"`, `"R2"`, or `"pseudoF"`.
#' @param by Whether to summarise by `reduction_level` or only globally.
#' @param center Central tendency target. Currently informative only for printing.
#' @param probs Probabilities for quantile summaries.
#' @param association Association method between ecological and inferential effects.
#'   One of `"spearman"`, `"pearson"`, or `"none"`.
#' @param ... Unused, reserved for future extensions.
#'
#' @return An object of class `summary_effect_size_data`.
#' @export
summary.effect_size_data <- function(
  object,
  inferential_var = c("omega2", "R2", "pseudoF"),
  by = c("reduction_level", "none"),
  center = c("median", "mean"),
  probs = c(0.25, 0.5, 0.75),
  association = c("spearman", "pearson", "none"),
  ...
) {
  validate_effect_size_data(object)

  inferential_var <- match.arg(inferential_var)
  by <- match.arg(by)
  center <- match.arg(center)
  association <- match.arg(association)

  df <- as.data.frame(object)

  if (!inferential_var %in% names(df)) {
    stop(
      "Column `",
      inferential_var,
      "` is not available in the object.",
      call. = FALSE
    )
  }

  meta <- list(
    n_rows = nrow(df),
    n_cols = ncol(df),
    reduction_levels = sort(unique(df$reduction_level)),
    n_reduction_levels = length(unique(df$reduction_level)),
    inferential_var = inferential_var,
    available_inferential = intersect(c("omega2", "R2", "pseudoF"), names(df)),
    has_k = "k" %in% names(df),
    has_m = "m" %in% names(df),
    has_n = "n" %in% names(df)
  )

  missing <- data.frame(
    variable = intersect(
      c("ecological_effect", "omega2", "R2", "pseudoF"),
      names(df)
    ),
    n_missing = vapply(
      intersect(c("ecological_effect", "omega2", "R2", "pseudoF"), names(df)),
      function(v) sum(is.na(df[[v]])),
      numeric(1)
    ),
    row.names = NULL
  )

  overall <- data.frame(
    metric = c("ecological_effect", inferential_var),
    rbind(
      .summarise_vector(df$ecological_effect, center = center, probs = probs),
      .summarise_vector(df[[inferential_var]], center = center, probs = probs)
    ),
    row.names = NULL,
    check.names = FALSE
  )

  by_reduction <- NULL
  if (by == "reduction_level") {
    by_reduction <- .summarise_by_reduction(
      df = df,
      inferential_var = inferential_var,
      center = center,
      probs = probs
    )
  }

  association_out <- NULL
  if (association != "none") {
    keep <- stats::complete.cases(df[, c("ecological_effect", inferential_var)])
    if (sum(keep) >= 3) {
      ct <- stats::cor.test(
        x = df$ecological_effect[keep],
        y = df[[inferential_var]][keep],
        method = association
      )

      association_out <- data.frame(
        method = association,
        estimate = unname(ct$estimate),
        p_value = ct$p.value,
        n = sum(keep),
        row.names = NULL
      )
    } else {
      association_out <- data.frame(
        method = association,
        estimate = NA_real_,
        p_value = NA_real_,
        n = sum(keep),
        row.names = NULL
      )
    }
  }

  out <- list(
    meta = meta,
    missing = missing,
    overall = overall,
    by_reduction = by_reduction,
    association = association_out,
    call = match.call()
  )

  class(out) <- c("summary_effect_size_data", "list")
  out
}

# =========================================================
# print.summary() method
# =========================================================

#' Print method for summary_effect_size_data objects
#'
#' @param x An object of class `summary_effect_size_data`.
#' @param digits Number of digits to print.
#' @param preview_n Number of rows to show for grouped summaries.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.summary_effect_size_data <- function(
  x,
  digits = max(3L, getOption("digits") - 2L),
  preview_n = 6L,
  ...
) {
  cat("<summary_effect_size_data>\n\n")

  cat("Inferential metric:", x$meta$inferential_var, "\n")
  cat("Observations:", x$meta$n_rows, "\n")
  cat("Reduction levels:", x$meta$n_reduction_levels, "\n")
  cat(
    "Available inferential metrics:",
    paste(x$meta$available_inferential, collapse = ", "),
    "\n"
  )

  cat("\nMissing values:\n")
  print(x$missing, row.names = FALSE)

  cat("\nOverall summary:\n")
  print(x$overall, digits = digits, row.names = FALSE)

  if (!is.null(x$association)) {
    cat("\nAssociation between ecological and inferential effects:\n")
    print(x$association, digits = digits, row.names = FALSE)
  }

  if (!is.null(x$by_reduction)) {
    cat("\nBy reduction level (preview):\n")
    print(
      utils::head(x$by_reduction, n = preview_n),
      digits = digits,
      row.names = FALSE
    )
  }

  invisible(x)
}

#' @keywords internal
.extract_pcoa_data <- function(x) {
  pcoa <- attr(x, "pcoa")

  if (is.null(pcoa)) {
    stop("This object does not contain stored PCoA data.", call. = FALSE)
  }

  if (!is.data.frame(pcoa)) {
    stop("`attr(x, 'pcoa')` must be a data.frame.", call. = FALSE)
  }

  required_cols <- c(
    "reduction_level",
    "step",
    "dat_sim",
    "k",
    "m",
    "n",
    "group",
    "Axis1",
    "Axis2"
  )

  missing_cols <- setdiff(required_cols, names(pcoa))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in stored PCoA data: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  pcoa
}

#' @keywords internal
.select_representative_pcoa <- function(
  x,
  pcoa_df,
  reduction_levels = NULL,
  selection = c("median", "first", "min", "max")
) {
  selection <- match.arg(selection)

  df <- as.data.frame(x)

  if (!is.null(reduction_levels)) {
    df <- df[df$reduction_level %in% reduction_levels, , drop = FALSE]
    pcoa_df <- pcoa_df[
      pcoa_df$reduction_level %in% reduction_levels,
      ,
      drop = FALSE
    ]
  }

  id_cols <- c("reduction_level", "step", "dat_sim", "k", "m", "n")

  selected <- switch(
    selection,
    median = {
      targets <- df |>
        dplyr::group_by(reduction_level) |>
        dplyr::summarise(
          target_eco = stats::median(ecological_effect, na.rm = TRUE),
          .groups = "drop"
        )

      df |>
        dplyr::inner_join(targets, by = "reduction_level") |>
        dplyr::mutate(.dist_to_target = abs(ecological_effect - target_eco)) |>
        dplyr::group_by(reduction_level) |>
        dplyr::slice_min(.dist_to_target, n = 1, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::select(dplyr::all_of(id_cols))
    },
    first = {
      df |>
        dplyr::group_by(reduction_level) |>
        dplyr::slice_head(n = 1) |>
        dplyr::ungroup() |>
        dplyr::select(dplyr::all_of(id_cols))
    },
    min = {
      df |>
        dplyr::group_by(reduction_level) |>
        dplyr::slice_min(ecological_effect, n = 1, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::select(dplyr::all_of(id_cols))
    },
    max = {
      df |>
        dplyr::group_by(reduction_level) |>
        dplyr::slice_max(ecological_effect, n = 1, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::select(dplyr::all_of(id_cols))
    }
  )

  out <- dplyr::inner_join(pcoa_df, selected, by = id_cols)

  out$panel_label <- factor(
    paste0("Reduction = ", out$reduction_level),
    levels = unique(paste0("Reduction = ", selected$reduction_level))
  )

  out
}

#' @keywords internal
.compute_pcoa_centroids <- function(df) {
  df |>
    dplyr::group_by(panel_label, reduction_level, group) |>
    dplyr::summarise(
      Axis1 = mean(Axis1, na.rm = TRUE),
      Axis2 = mean(Axis2, na.rm = TRUE),
      .groups = "drop"
    )
}
