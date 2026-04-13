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
#' @param steps =====PENDING====
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
  if (model != "single.factor") {
    stop(
      "`sim_ES()` currently supports only `model = 'single.factor'`. Nested designs are not implemented yet."
    )
  }

  # read data and store it in two objects, one for H0 and one for Ha ----
  datH0 <- data
  datH0[, 1] <- as.factor("zero")
  datHa <- data
  datHa[, 1] <- as.factor(data[, 1])
  a <- nlevels(datHa[, 1])
  M <- a

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
  resultsHa[, 3] <- M
  resultsHa[, 4] <- rep(seq(2, n), times = 1, each = k)

  Y <- cbind(1:(N * M))
  YPU <- as.numeric(gl(M, N))
  mm <- resultsHa[, 3]
  nn <- resultsHa[, 3] * resultsHa[, 4]

  # Calculate H0 simulation parameters
  parH0 <- SSP::assempar(data = datH0, type = type, Sest.method = Sest.method)
  parHa <- SSP::assempar(data = datHa, type = type, Sest.method = Sest.method)

  ## Calculate species similarity percentage from Ha ----
  sppContribution <- use_simper(datHa)

  ## Output container ----
  resultOut <- vector("list", length = NN * (steps + 1))

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

  # Here we start the loop for species contribution attenuation ----
  n_spp <- dim(sppContribution)[1]
  propSpp <- 1 / steps

  for (st in 0:steps) {
    # Calculate number of species to alter
    removing <- st * propSpp

    adjustN <- ceiling(nrow(sppContribution) * removing)

    if (adjustN != 0) {
      adjusting <- sppContribution[1:adjustN, 1]

      # Obtain H0 parameters from null hypothesis
      zeroParameter <- parH0$par
      zeroParameter <- zeroParameter[
        zeroParameter$Species %in% adjusting,
        c(1, 3)
      ] # CONFIRMAR CON EDLIN QUE SE USA fw

      # Alter the Ha parameters with the corresponding null hypothesis version
      parHaTemp <- parHa
      parHaTemp$par[
        parHaTemp$par$Species %in% adjusting,
        c(3)
      ] <- zeroParameter[,
        2
      ]
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
          "balanced_sampling_es",
          "dbmanova_oneway",
          "calc_dist"
        ),
        environment = asNamespace("ecocbo")
      )

      # Executing the loop in parallel
      result1 <- parabar::par_lapply(
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

      parabar::stop_backend(cl)
    } else {
      pb <- txtProgressBar(max = NN, style = 3)
      result1 <- vector("list", length = NN)

      for (i in seq_len(NN)) {
        # Performs the operation iteratively in a for loop
        result1[[i]] <- balanced_sampling_es(
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

        # Updating the progress bar
        setTxtProgressBar(pb, i)
      }
      close(pb)
    }

    idx <- (st * NN + 1):((st + 1) * NN)
    resultOut[idx] <- lapply(seq_len(NN), function(i) {
      cur <- result1[[i]]
      data.frame(
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
      ) #|>
      # dplyr::mutate(
      #   centroid_dist_matrix = list(cur$centroid_dist_matrix),
      #   infer_table = list(cur$infer_table)
      # )
    })

    cat("Step ", st, ": Done! - ", steps - st, " remaining")
  }

  resultOut <- dplyr::bind_rows(resultOut)
  # class(resultOut) <- c("effect_size_data", class(resultOut))
  new_effect_size_data(resultOut)

  return(resultOut)
}


# #-------------------------------------------
# ## S3Methods plot()
# #-------------------------------------------
# #' Plot effect_size_data
# #'
# #' @method plot effect_size_data
# #'
# #' @param x An object of class \code{effect_size_data}.
# #' @param ... Additional arguments, ignored.
# #'
# #' @return A ggplot representation of the \code{effect_size_data} object.
# #' @export

# plot.effect_size_data <- function(x, type = c("stacked", "scatter"), ...) {
#   if (!all(c("reduction_level", "ecological_effect") %in% colnames(x))) {
#     stop("`plot.effect_size_data()` expects redesigned output from `sim_ES()`.")
#   }

#   p1 <- x |>
#     as.data.frame() |>
#     dplyr::group_by(reduction_level) |>
#     dplyr::summarise(
#       media = mean(ecological_effect, na.rm = TRUE),
#       mediana = median(ecological_effect, na.rm = TRUE),
#       q25 = quantile(ecological_effect, 0.25, na.rm = TRUE),
#       q75 = quantile(ecological_effect, 0.75, na.rm = TRUE),
#       .groups = "drop"
#     ) |>
#     ggplot2::ggplot(ggplot2::aes(x = reduction_level, y = mediana)) +
#     ggplot2::geom_ribbon(ggplot2::aes(ymin = q25, ymax = q75), alpha = 0.3) +
#     ggplot2::geom_line() +
#     ggplot2::geom_point() +
#     ggplot2::scale_x_continuous(name = "% removed SIMPER contribution") +
#     ggplot2::scale_y_continuous(name = "Ecological effect size") +
#     ggplot2::theme_bw() +
#     ggplot2::theme(
#       panel.grid.minor.x = ggplot2::element_blank(),
#       panel.border = ggplot2::element_rect(linewidth = 0.4),
#       axis.ticks = ggplot2::element_line(linewidth = 0.2)
#     )

#   print(p1)
#   invisible(p1)
# }

# ==============

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
  show_hulls = FALSE,
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
      show_hulls = show_hulls,
      point_alpha = point_alpha,
      point_size = point_size,
      centroid_size = centroid_size,
      facet_ncol = facet_ncol
    )
  )
}

# =========================================================
# PCoA helpers
# =========================================================

#' @keywords internal
.extract_pcoa_data <- function(x) {
  pcoa <- attr(x, "pcoa")

  if (is.null(pcoa)) {
    stop(
      "This effect_size_data object does not contain PCoA data in `attr(x, 'pcoa')`.",
      call. = FALSE
    )
  }

  if (is.list(pcoa) && !is.data.frame(pcoa)) {
    if (!all(vapply(pcoa, is.data.frame, logical(1)))) {
      stop(
        "If `attr(x, 'pcoa')` is a list, all elements must be data.frames.",
        call. = FALSE
      )
    }
    pcoa <- dplyr::bind_rows(pcoa)
  }

  if (!is.data.frame(pcoa)) {
    stop(
      "`attr(x, 'pcoa')` must be a data.frame or a list of data.frames.",
      call. = FALSE
    )
  }

  # Flexible detection of axis names
  axis1_candidates <- c("Axis1", "PCoA1", "Dim1", "MDS1")
  axis2_candidates <- c("Axis2", "PCoA2", "Dim2", "MDS2")
  group_candidates <- c("group", "Group", "grp")

  axis1_col <- intersect(axis1_candidates, names(pcoa))[1]
  axis2_col <- intersect(axis2_candidates, names(pcoa))[1]
  group_col <- intersect(group_candidates, names(pcoa))[1]

  if (is.na(axis1_col) || is.na(axis2_col)) {
    stop(
      "PCoA data must contain axis columns such as `Axis1`/`Axis2` or `PCoA1`/`PCoA2`.",
      call. = FALSE
    )
  }

  if (is.na(group_col)) {
    stop(
      "PCoA data must contain a grouping column, preferably named `group`.",
      call. = FALSE
    )
  }

  names(pcoa)[names(pcoa) == axis1_col] <- "Axis1"
  names(pcoa)[names(pcoa) == axis2_col] <- "Axis2"
  names(pcoa)[names(pcoa) == group_col] <- "group"

  required_cols <- c("reduction_level", "Axis1", "Axis2", "group")
  missing_cols <- setdiff(required_cols, names(pcoa))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required PCoA columns: ",
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

  id_candidates <- c("reduction_level", "k", "m", "n")
  join_cols <- intersect(id_candidates, intersect(names(df), names(pcoa_df)))

  if (!"reduction_level" %in% join_cols) {
    stop(
      "`reduction_level` must exist in both the main object and PCoA data.",
      call. = FALSE
    )
  }

  extra_id_cols <- setdiff(join_cols, "reduction_level")

  if (length(extra_id_cols) == 0 && any(duplicated(df$reduction_level))) {
    stop(
      "Cannot select a representative simulation per reduction_level because no simulation identifiers ",
      "(such as `k`, `m`, or `n`) are shared between the main table and the PCoA data.",
      call. = FALSE
    )
  }

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
        dplyr::slice_min(
          order_by = .dist_to_target,
          n = 1,
          with_ties = FALSE
        ) |>
        dplyr::ungroup() |>
        dplyr::select(dplyr::all_of(join_cols))
    },
    first = {
      df |>
        dplyr::group_by(reduction_level) |>
        dplyr::slice_head(n = 1) |>
        dplyr::ungroup() |>
        dplyr::select(dplyr::all_of(join_cols))
    },
    min = {
      df |>
        dplyr::group_by(reduction_level) |>
        dplyr::slice_min(
          order_by = ecological_effect,
          n = 1,
          with_ties = FALSE
        ) |>
        dplyr::ungroup() |>
        dplyr::select(dplyr::all_of(join_cols))
    },
    max = {
      df |>
        dplyr::group_by(reduction_level) |>
        dplyr::slice_max(
          order_by = ecological_effect,
          n = 1,
          with_ties = FALSE
        ) |>
        dplyr::ungroup() |>
        dplyr::select(dplyr::all_of(join_cols))
    }
  )

  out <- dplyr::inner_join(pcoa_df, selected, by = join_cols)

  reduction_order <- selected |>
    dplyr::distinct(reduction_level) |>
    dplyr::arrange(reduction_level) |>
    dplyr::pull(reduction_level)

  out$panel_label <- factor(
    paste0("Reduction = ", out$reduction_level),
    levels = paste0("Reduction = ", reduction_order)
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

#' @keywords internal
.compute_pcoa_hulls <- function(df) {
  df |>
    dplyr::group_by(panel_label, reduction_level, group) |>
    dplyr::filter(dplyr::n() >= 3) |>
    dplyr::slice(chull(Axis1, Axis2)) |>
    dplyr::ungroup()
}

#' @keywords internal
.plot_effect_size_ordination <- function(
  x,
  reduction_levels = NULL,
  selection = c("median", "first", "min", "max"),
  show_centroids = TRUE,
  label_centroids = TRUE,
  show_hulls = FALSE,
  point_alpha = 0.7,
  point_size = 1.6,
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
  )

  if (show_hulls) {
    hulls <- .compute_pcoa_hulls(plot_df)

    if (nrow(hulls) > 0) {
      p <- p +
        ggplot2::geom_polygon(
          data = hulls,
          ggplot2::aes(fill = group, group = interaction(panel_label, group)),
          alpha = 0.12,
          colour = NA,
          inherit.aes = FALSE
        )
    }
  }

  p <- p +
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
      colour = "Group",
      fill = "Group"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(linewidth = 0.4),
      axis.ticks = ggplot2::element_line(linewidth = 0.2)
    )
}
