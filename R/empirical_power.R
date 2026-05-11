#' Estimate empirical power by rarefying pilot data
#'
#' @description
#' Estimates empirical detection probability from pilot community data by
#' repeatedly subsampling a nested sampling design and applying a PERMANOVA-based
#' test to each subsample. The function is intended as a field-data-only
#' rarefaction procedure: it evaluates sampling efforts that are smaller than or
#' equal to the effort already represented in the pilot data.
#'
#' @details
#' The function assumes a nested design with a fixed factor represented by
#' `sector_col` and a nested sampling unit represented by `site_col`. All levels
#' of the fixed factor are retained in every iteration. For each sampling effort,
#' the function selects `m` sites within each sector and `n` subsamples within
#' each selected site. The values of `m` and `n` are defined by the empirical
#' effort grid constructed from the observed pilot-data structure.
#'
#' For each `m` by `n` combination, the function performs `n_iter` independent
#' subsampling iterations. In each iteration, the selected subset is passed to
#' the internal PERMANOVA workflow through [run_pilot_rarefaction()]. Empirical
#' power is estimated as the proportion of successful subsamples with
#' `p_A < alpha`, where `p_A` corresponds to the fixed sector effect. The same
#' summary is also computed for the nested `site(sector)` component as
#' `power_BA`.
#'
#' This analysis should be interpreted as interpolation within the observed pilot
#' effort, not as extrapolation beyond the available data. It is useful for
#' comparing simulated power curves against curves obtained exclusively from the
#' original field data.
#'
#' @param comm_pilot A data frame, matrix, or object coercible to a data frame
#'   containing the community data. Rows must correspond to sampling units and
#'   columns to taxa, species, or response variables used to compute community
#'   dissimilarities.
#' @param factEnv A data frame containing the experimental or environmental
#'   factors associated with the rows of `comm_pilot`.
#' @param sector_col A character string giving the name of the column in
#'   `factEnv` that identifies the fixed factor of interest. This column is
#'   internally renamed to `sector`.
#' @param site_col A character string giving the name of the column in `factEnv`
#'   that identifies the nested sampling unit within each sector. This column is
#'   internally renamed to `site`.
#' @param method A character string specifying the dissimilarity method to be
#'   used by the PERMANOVA workflow. The default is `"bray"`.
#' @param transformation A character string specifying the transformation applied
#'   to the community data before calculating dissimilarities. Supported values
#'   depend on the internal PERMANOVA helper, but typically include `"none"`,
#'   `"square root"`, `"fourth root"`, `"Log (X+1)"`, and `"P/A"`.
#' @param dummy Logical. If `TRUE`, a dummy variable is added before calculating
#'   dissimilarities. This can help avoid undefined dissimilarities in rows with
#'   no observed taxa.
#' @param model A character string specifying the PERMANOVA model. The current
#'   workflow is intended for `"nested.symmetric"`.
#' @param n_iter Integer. Number of independent subsampling iterations to perform
#'   for each sampling-effort combination.
#' @param permutations Integer. Number of permutations used in each PERMANOVA
#'   test.
#' @param seed Integer. Random seed used for reproducible subsampling and
#'   permutation workflows.
#' @param min_m Integer. Minimum number of sites per sector to evaluate in the
#'   rarefaction grid. The default is `2`.
#' @param min_n Integer. Minimum number of subsamples per site to evaluate in the
#'   rarefaction grid. The default is `2`.
#' @param expected_n_sectors Optional integer. Expected number of sectors. If
#'   supplied, the function stops when the number of observed sectors differs
#'   from this value.
#' @param alpha Numeric. Significance threshold used to calculate empirical
#'   detection probability. The default is `0.05`.
#' @param parallel Logical. If `TRUE`, the rarefaction workflow is executed in
#'   parallel using the internal parallel implementation of
#'   [run_pilot_rarefaction()].
#' @param workers Integer. Number of parallel workers to use when
#'   `parallel = TRUE`. By default, it uses one fewer than the number of
#'   available cores.
#' @param progress Logical. If `TRUE`, a progress bar is displayed during the
#'   rarefaction workflow.
#' @param make_plot Logical. If `TRUE`, a `ggplot2` object showing the empirical
#'   rarefaction curve for the fixed sector effect is returned in the output.
#'
#' @return
#' An object of class `"pilot_rarefaction_workflow"`, implemented as a list with
#' the following components:
#'
#' \describe{
#'   \item{call}{The matched function call.}
#'   \item{parameters}{A list with the main analytical parameters used in the
#'   workflow, including the dissimilarity method, transformation, number of
#'   iterations, number of permutations, alpha level, and parallel settings.}
#'   \item{design}{A list describing the empirical design detected from
#'   `factEnv`, including the number of sectors, maximum number of sites per
#'   sector, maximum number of subsamples per site, balance diagnostics, and
#'   observed counts by sector and site.}
#'   \item{factEnv_nested}{A standardized factor data frame with columns
#'   `sector`, `site`, and `site_id`.}
#'   \item{effort_grid}{A data frame with all evaluated sampling-effort
#'   combinations, including `m`, `n`, total sites, total subsamples, effort
#'   labels, and the proportion of the full pilot effort.}
#'   \item{results}{A data frame with one row per subsampling iteration and
#'   effort combination. It includes pseudo-F statistics, p-values, error
#'   diagnostics, and sampling-effort metadata.}
#'   \item{power_summary}{A data frame summarising empirical detection
#'   probability and pseudo-F distributions for each sampling-effort
#'   combination.}
#'   \item{plot}{A `ggplot2` object with the empirical rarefaction curve when
#'   `make_plot = TRUE`; otherwise `NULL`.}
#' }
#'
#' @seealso
#' [run_pilot_rarefaction()]
#'
#' @examples
#' \donttest{
#' if (requireNamespace("SSP", quietly = TRUE)) {
#'   pilot_epibionts <- empirical_power(
#'     comm_pilot = SSP::epibionts[, -c(1, 2)],
#'     factEnv = SSP::epibionts[, c(1, 2)],
#'     sector_col = "sector",
#'     site_col = "site",
#'     method = "bray",
#'     transformation = "none",
#'     dummy = TRUE,
#'     model = "nested.symmetric",
#'     n_iter = 5,
#'     permutations = 9,
#'     seed = 123,
#'     min_m = 2,
#'     min_n = 2,
#'     expected_n_sectors = 3,
#'     alpha = 0.05,
#'     parallel = FALSE,
#'     progress = FALSE,
#'     make_plot = TRUE
#'   )
#'
#'   pilot_epibionts$power_summary
#'   pilot_epibionts$plot
#' }
#' }
#'
#' @export

empirical_power <- function(
  comm_pilot,
  factEnv,
  sector_col,
  site_col,
  method = "bray",
  transformation = "none",
  dummy = TRUE,
  model = "nested.symmetric",
  n_iter = 199,
  permutations = 199,
  seed = 123,
  min_m = 2,
  min_n = 2,
  expected_n_sectors = NULL,
  alpha = 0.05,
  parallel = FALSE,
  workers = max(1, parallelly::availableCores() - 1),
  progress = TRUE,
  make_plot = TRUE
) {
  # ------------------------------------------------------------
  # 1. Validating data
  # ------------------------------------------------------------

  comm_pilot <- as.data.frame(comm_pilot)
  factEnv <- as.data.frame(factEnv)

  if (nrow(comm_pilot) != nrow(factEnv)) {
    stop(
      "`comm_pilot` and `factEnv` must have the same number of rows."
    )
  }

  if (!sector_col %in% names(factEnv)) {
    stop("`sector_col` was not found in `factEnv`.")
  }

  if (!site_col %in% names(factEnv)) {
    stop("`site_col` was not found in `factEnv`.")
  }

  # ------------------------------------------------------------
  # 2. Environmental factors
  # ------------------------------------------------------------

  factEnv_nested <- factEnv |>
    dplyr::transmute(
      sector = as.factor(.data[[sector_col]]),
      site = as.factor(.data[[site_col]]),
      site_id = interaction(sector, site, drop = TRUE, sep = "_")
    )

  site_counts <- factEnv_nested |>
    dplyr::count(sector, site_id, name = "n_subsamples")

  sites_by_sector <- site_counts |>
    dplyr::count(sector, name = "n_sites")

  n_sectors <- dplyr::n_distinct(factEnv_nested$sector)
  m_max <- min(sites_by_sector$n_sites)
  n_max <- min(site_counts$n_subsamples)

  if (!is.null(expected_n_sectors)) {
    if (n_sectors != expected_n_sectors) {
      stop(
        "Expected ",
        expected_n_sectors,
        " sectors, but found ",
        n_sectors,
        "."
      )
    }
  }

  if (m_max < min_m) {
    stop(
      "The minimum number of sites per sector is ",
      m_max,
      ", which is lower than `min_m = ",
      min_m,
      "`."
    )
  }

  if (n_max < min_n) {
    stop(
      "The minimum number of subsamples per site is ",
      n_max,
      ", which is lower than `min_n = ",
      min_n,
      "`."
    )
  }

  balanced_sites <- length(unique(sites_by_sector$n_sites)) == 1L
  balanced_subsamples <- length(unique(site_counts$n_subsamples)) == 1L

  # ------------------------------------------------------------
  # 3. Create the effort grid
  # ------------------------------------------------------------

  effort_grid <- tidyr::expand_grid(
    m = seq(m_max, min_m, by = -1),
    n = seq(n_max, min_n, by = -1)
  ) |>
    dplyr::mutate(
      n_sectors = n_sectors,
      total_sites = n_sectors * m,
      total_subsamples = n_sectors * m * n,
      prop_full_pilot = total_subsamples / max(total_subsamples),
      effort_id = sprintf("m%02d_n%02d", m, n),
      effort_label = sprintf(
        "m = %d sites/sector; n = %d subsamples/site",
        m,
        n
      )
    ) |>
    dplyr::arrange(dplyr::desc(m), dplyr::desc(n))

  # ------------------------------------------------------------
  # 4. Execute the empirical rarefaction
  # ------------------------------------------------------------

  results <- run_pilot_rarefaction(
    comm_pilot = comm_pilot,
    factEnvP = factEnv_nested,
    pilot_effort_grid = effort_grid,
    n_iter = n_iter,
    permutations = permutations,
    seed = seed,
    method = method,
    transformation = transformation,
    dummy = dummy,
    model = model,
    parallel = parallel,
    workers = workers,
    progress = progress
  )

  # ------------------------------------------------------------
  # 5. Summarise the empirical power
  # ------------------------------------------------------------

  power_summary <- results |>
    dplyr::group_by(
      m,
      n,
      total_sites,
      total_subsamples,
      prop_full_pilot
    ) |>
    dplyr::summarise(
      n_success = sum(is.na(error)),
      n_failed = sum(!is.na(error)),
      power_A = mean(p_A < alpha, na.rm = TRUE),
      power_BA = mean(p_BA < alpha, na.rm = TRUE),
      pseudoF_A_median = stats::median(pseudoF_A, na.rm = TRUE),
      pseudoF_A_q025 = stats::quantile(pseudoF_A, 0.025, na.rm = TRUE),
      pseudoF_A_q975 = stats::quantile(pseudoF_A, 0.975, na.rm = TRUE),
      pseudoF_BA_median = stats::median(pseudoF_BA, na.rm = TRUE),
      pseudoF_BA_q025 = stats::quantile(pseudoF_BA, 0.025, na.rm = TRUE),
      pseudoF_BA_q975 = stats::quantile(pseudoF_BA, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(m), dplyr::desc(n))

  # ------------------------------------------------------------
  # 6. Plot
  # ------------------------------------------------------------

  p_power <- NULL

  if (make_plot) {
    p_power <- ggplot2::ggplot(
      power_summary,
      ggplot2::aes(
        x = n,
        y = power_A,
        group = factor(m),
        colour = factor(m)
      )
    ) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::scale_x_continuous(
        breaks = sort(unique(power_summary$n))
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, 1)
      ) +
      ggplot2::labs(
        x = "Subsamples per site",
        y = "Empirical detection probability",
        colour = "Sites per sector",
        title = "Pilot-data empirical rarefaction curve",
        subtitle = "Nested subsampling preserving all sector levels"
      ) +
      ggplot2::theme_bw()
  }

  # ------------------------------------------------------------
  # 7. Output
  # ------------------------------------------------------------

  out <- list(
    call = match.call(),
    parameters = list(
      method = method,
      transformation = transformation,
      dummy = dummy,
      model = model,
      n_iter = n_iter,
      permutations = permutations,
      seed = seed,
      min_m = min_m,
      min_n = min_n,
      alpha = alpha,
      parallel = parallel,
      workers = workers
    ),
    design = list(
      n_sectors = n_sectors,
      m_max = m_max,
      n_max = n_max,
      balanced_sites = balanced_sites,
      balanced_subsamples = balanced_subsamples,
      sites_by_sector = sites_by_sector,
      site_counts = site_counts
    ),
    factEnv_nested = factEnv_nested,
    effort_grid = effort_grid,
    results = results,
    power_summary = power_summary,
    plot = p_power
  )

  class(out) <- c("pilot_rarefaction_workflow", class(out))

  out
}


#' Print a pilot rarefaction workflow object
#'
#' @description
#' Prints a compact summary of an object returned by [empirical_power()].
#'
#' @param x An object of class `"pilot_rarefaction_workflow"`, usually returned
#'   by [empirical_power()].
#' @param ... Additional arguments passed to or from other methods. Currently
#'   ignored.
#'
#' @return
#' Invisibly returns `x`.
#'
#' @method print pilot_rarefaction_workflow
#' @export

print.pilot_rarefaction_workflow <- function(x, ...) {
  cat("Pilot-data empirical rarefaction workflow\n")
  cat("--------------------------------------------------\n")
  cat("Sectors:             ", x$design$n_sectors, "\n")
  cat("Max sites/sector:    ", x$design$m_max, "\n")
  cat("Max subsamples/site: ", x$design$n_max, "\n")
  cat("Effort combinations: ", nrow(x$effort_grid), "\n")
  cat("Iterations/effort:   ", x$parameters$n_iter, "\n")
  cat("Permutations/test:   ", x$parameters$permutations, "\n")
  cat("Alpha:               ", x$parameters$alpha, "\n")
  cat("Parallel:            ", x$parameters$parallel, "\n")
  cat("Workers:             ", x$parameters$workers, "\n")
  cat("Balanced sites:      ", x$design$balanced_sites, "\n")
  cat("Balanced subsamples: ", x$design$balanced_subsamples, "\n")

  invisible(x)
}
