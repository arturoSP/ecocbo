# =========================================================
# effect_size_data S3 class
# =========================================================

#' @keywords internal
.detect_effect_size_variant <- function(data, model = NULL) {
  if (!is.null(model) && identical(model, "nested.symmetric")) {
    return("nested")
  }

  if (all(c("ecological_effect_A", "ecological_effect_BA") %in% names(data))) {
    return("nested")
  }

  "single"
}

#' Create an effect_size_data object
#'
#' @param data A data frame with the effect size results.
#' @param centroid_distances Optional list or object storing centroid distance matrices.
#' @param pcoa Optional list or object storing PCoA results.
#' @param call Optional matched call.
#' @param model Optional model label. Expected values include
#'   `"single.factor"` and `"nested.symmetric"`.
#'
#' @return An object of class `effect_size_data`, or subclass
#'   `effect_size_data_nested` for nested outputs.
#' @export
new_effect_size_data <- function(
  data,
  centroid_distances = NULL,
  pcoa = NULL,
  call = NULL,
  model = NULL
) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)

  variant <- .detect_effect_size_variant(data, model = model)
  validate_effect_size_data(data, variant = variant)

  if (is.null(model)) {
    model <- if (variant == "nested") "nested.symmetric" else "single.factor"
  }

  cls <- if (variant == "nested") {
    c("effect_size_data_nested", "effect_size_data", class(data))
  } else {
    c("effect_size_data", class(data))
  }

  class(data) <- unique(cls)
  attr(data, "centroid_distances") <- centroid_distances
  attr(data, "pcoa") <- pcoa
  attr(data, "call") <- call
  attr(data, "model") <- model
  attr(data, "variant") <- variant

  data
}

#' Validate an effect_size_data object
#'
#' @param x A data frame or `effect_size_data` object.
#' @param variant Optional variant: `"single"` or `"nested"`.
#'
#' @return Invisibly returns `TRUE` if validation passes.
#' @keywords internal
validate_effect_size_data <- function(x, variant = NULL) {
  if (!is.data.frame(x)) {
    stop("`x` must be a data.frame or effect_size_data object.", call. = FALSE)
  }

  if (is.null(variant)) {
    variant <- .detect_effect_size_variant(x, model = attr(x, "model"))
  }

  required_common <- c("reduction_level")
  missing_common <- setdiff(required_common, names(x))

  if (length(missing_common) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_common, collapse = ", "),
      call. = FALSE
    )
  }

  if (variant == "single") {
    required_cols <- c("ecological_effect")
    missing_cols <- setdiff(required_cols, names(x))

    if (length(missing_cols) > 0) {
      stop(
        "Missing required columns for single-factor object: ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }

    if (!is.numeric(x$ecological_effect)) {
      stop("`ecological_effect` must be numeric.", call. = FALSE)
    }

    inferential_candidates <- c("omega2", "R2", "pseudoF")
  } else {
    required_cols <- c("ecological_effect_A", "ecological_effect_BA")
    missing_cols <- setdiff(required_cols, names(x))

    if (length(missing_cols) > 0) {
      stop(
        "Missing required columns for nested object: ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }

    if (
      !is.numeric(x$ecological_effect_A) || !is.numeric(x$ecological_effect_BA)
    ) {
      stop(
        "`ecological_effect_A` and `ecological_effect_BA` must be numeric.",
        call. = FALSE
      )
    }

    inferential_candidates <- c(
      "omega2_A",
      "omega2_BA",
      "R2_A",
      "R2_BA",
      "pseudoF_A",
      "pseudoF_BA"
    )
  }

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
  class(x) <- setdiff(
    class(x),
    c("effect_size_data_nested", "effect_size_data")
  )
  x
}

#' @export
as.data.frame.effect_size_data <- function(x, ...) {
  class(x) <- setdiff(class(x), "effect_size_data")
  x
}

# =========================================================
# print() methods
# =========================================================

#' @keywords internal
.print_effect_size_data_common <- function(
  x,
  digits = max(3L, getOption("digits") - 2L),
  preview_n = 6L,
  variant = c("single", "nested")
) {
  variant <- match.arg(variant)
  validate_effect_size_data(x, variant = variant)

  df <- as.data.frame(x)
  n_rows <- nrow(df)
  n_cols <- ncol(df)

  reduction_n <- length(unique(df$reduction_level))
  reduction_range <- tryCatch(
    range(df$reduction_level, na.rm = TRUE),
    error = function(e) NULL
  )

  has_pcoa <- !is.null(attr(x, "pcoa"))
  model <- attr(x, "model")

  if (variant == "single") {
    inferential_candidates <- c("omega2", "R2", "pseudoF")
    inferential_present <- intersect(inferential_candidates, names(df))
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

    cat("<effect_size_data>\n")
  } else {
    inferential_candidates <- c(
      "omega2_A",
      "omega2_BA",
      "R2_A",
      "R2_BA",
      "pseudoF_A",
      "pseudoF_BA"
    )
    inferential_present <- intersect(inferential_candidates, names(df))
    preview_cols <- intersect(
      c(
        "reduction_level",
        "k",
        "m",
        "n",
        "ecological_effect_A",
        "ecological_effect_BA",
        "omega2_A",
        "omega2_BA",
        "R2_A",
        "R2_BA",
        "pseudoF_A",
        "pseudoF_BA",
        "ok",
        "error_message"
      ),
      names(df)
    )

    cat("<effect_size_data_nested>\n")
  }

  cat("Model:", model, "\n")
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

  if (variant == "nested") {
    if ("n_groups_A" %in% names(df)) {
      cat(
        "Groups in A:",
        paste(sort(unique(stats::na.omit(df$n_groups_A))), collapse = ", "),
        "\n"
      )
    }
    if ("n_groups_BA" %in% names(df)) {
      cat(
        "Groups in B(A):",
        paste(sort(unique(stats::na.omit(df$n_groups_BA))), collapse = ", "),
        "\n"
      )
    }
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

  if ("ok" %in% names(df)) {
    cat(
      "Successful iterations:",
      sum(isTRUE(df$ok) | df$ok %in% TRUE, na.rm = TRUE),
      "/",
      n_rows,
      "\n"
    )
  }

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
  .print_effect_size_data_common(
    x = x,
    digits = digits,
    preview_n = preview_n,
    variant = "single"
  )
}

#' Print method for effect_size_data_nested objects
#'
#' @param x An object of class `effect_size_data_nested`.
#' @param digits Number of digits to print.
#' @param preview_n Number of rows to preview.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.effect_size_data_nested <- function(
  x,
  digits = max(3L, getOption("digits") - 2L),
  preview_n = 6L,
  ...
) {
  .print_effect_size_data_common(
    x = x,
    digits = digits,
    preview_n = preview_n,
    variant = "nested"
  )
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

#' @keywords internal
.summarise_effect_size_data_nested <- function(
  x,
  component = c("A", "BA"),
  inferential_var = c("omega2", "R2", "pseudoF"),
  summary_stat = c("median", "mean"),
  interval = c("iqr", "sd")
) {
  component <- match.arg(component)
  inferential_var <- match.arg(inferential_var)
  summary_stat <- match.arg(summary_stat)
  interval <- match.arg(interval)

  df <- as.data.frame(x)

  eco_col <- paste0("ecological_effect_", component)
  inf_col <- paste0(inferential_var, "_", component)

  if (!eco_col %in% names(df)) {
    stop(
      "Column `",
      eco_col,
      "` is not available in the object.",
      call. = FALSE
    )
  }
  if (!inf_col %in% names(df)) {
    stop(
      "Column `",
      inf_col,
      "` is not available in the object.",
      call. = FALSE
    )
  }

  dplyr::group_by(df, reduction_level) |>
    dplyr::summarise(
      eco_center = .es_center(.data[[eco_col]], summary_stat),
      eco_low = .es_low(.data[[eco_col]], summary_stat, interval),
      eco_high = .es_high(.data[[eco_col]], summary_stat, interval),
      inf_center = .es_center(.data[[inf_col]], summary_stat),
      inf_low = .es_low(.data[[inf_col]], summary_stat, interval),
      inf_high = .es_high(.data[[inf_col]], summary_stat, interval),
      .groups = "drop"
    )
}

#' @keywords internal
.plot_effect_size_nested_stacked <- function(
  x,
  component = c("A", "BA"),
  inferential_var = c("omega2", "R2", "pseudoF"),
  summary_stat = c("median", "mean"),
  interval = c("iqr", "sd"),
  ribbon_alpha = 0.25,
  line_size = 0.4,
  point_size = 1.5
) {
  component <- match.arg(component)
  inferential_var <- match.arg(inferential_var)
  summary_stat <- match.arg(summary_stat)
  interval <- match.arg(interval)

  df_sum <- .summarise_effect_size_data_nested(
    x = x,
    component = component,
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

  y_lab_eco <- if (component == "A") {
    "Ecological effect size (A)"
  } else {
    "Ecological effect size (B(A))"
  }

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
      y = y_lab_eco
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
      y = if (component == "A") {
        paste0(as.character(y_lab_inf), " (A)")
      } else {
        paste0(as.character(y_lab_inf), " (B(A))")
      }
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
      metric = y_lab_eco
    ),
    data.frame(
      reduction_level = df_sum$reduction_level,
      center = df_sum$inf_center,
      low = df_sum$inf_low,
      high = df_sum$inf_high,
      metric = paste0(inferential_var, "_", component)
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

#' @keywords internal
.plot_effect_size_nested_scatter <- function(
  x,
  component = c("A", "BA"),
  inferential_var = c("omega2", "R2", "pseudoF"),
  colour_by = c("reduction_level", "none"),
  scatter_data = c("summary", "raw"),
  summary_stat = c("median", "mean"),
  point_alpha = 0.7,
  point_size = 1.8,
  add_smooth = FALSE
) {
  component <- match.arg(component)
  inferential_var <- match.arg(inferential_var)
  colour_by <- match.arg(colour_by)
  scatter_data <- match.arg(scatter_data)
  summary_stat <- match.arg(summary_stat)

  df <- as.data.frame(x)

  eco_col <- paste0("ecological_effect_", component)
  inf_col <- paste0(inferential_var, "_", component)

  if (!eco_col %in% names(df)) {
    stop(
      "Column `",
      eco_col,
      "` is not available in the object.",
      call. = FALSE
    )
  }
  if (!inf_col %in% names(df)) {
    stop(
      "Column `",
      inf_col,
      "` is not available in the object.",
      call. = FALSE
    )
  }

  if (scatter_data == "summary") {
    df <- dplyr::group_by(df, reduction_level) |>
      dplyr::summarise(
        ecological_effect = .es_center(.data[[eco_col]], summary_stat),
        inferential = .es_center(.data[[inf_col]], summary_stat),
        .groups = "drop"
      )
  } else {
    df <- dplyr::transmute(
      df,
      reduction_level = reduction_level,
      ecological_effect = .data[[eco_col]],
      inferential = .data[[inf_col]]
    )
  }

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = ecological_effect,
      y = inferential,
      colour = if (colour_by == "reduction_level") reduction_level else NULL
    )
  ) +
    ggplot2::geom_point(alpha = point_alpha, size = point_size) +
    ggplot2::labs(
      x = if (component == "A") {
        "Ecological effect size (A)"
      } else {
        "Ecological effect size (B(A))"
      },
      y = switch(
        inferential_var,
        omega2 = if (component == "A") {
          expression(omega[A]^2)
        } else {
          expression(omega[B(A)]^2)
        },
        R2 = if (component == "A") {
          expression(R[A]^2)
        } else {
          expression(R[B(A)]^2)
        },
        pseudoF = if (component == "A") "pseudo-F (A)" else "pseudo-F (B(A))"
      ),
      colour = if (colour_by == "reduction_level") "Reduction level" else NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(linewidth = 0.4),
      axis.ticks = ggplot2::element_line(linewidth = 0.2)
    )

  if (colour_by == "none") {
    p <- p + ggplot2::guides(colour = "none")
  }

  if (add_smooth) {
    p <- p +
      ggplot2::geom_smooth(
        method = "lm",
        se = FALSE,
        linewidth = 0.4,
        inherit.aes = TRUE
      )
  }

  p
}

# =========================================================
# plot() method for nested objects
# =========================================================

#' Plot method for effect_size_data_nested objects
#'
#' @param x An object of class `effect_size_data_nested`.
#' @param type Type of plot. `"stacked"` returns the two-panel
#'   summary plot; `"scatter"` returns the inferential vs ecological plot;
#'   `"ordination"` returns a PCoA display.
#' @param component Which nested component to visualise. One of `"A"` or `"BA"`.
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
#' @param reduction_levels Reduction levels to keep in ordination plots.
#' @param selection How to select representative ordinations per reduction level.
#' @param show_centroids Logical; show centroids in ordination plot.
#' @param label_centroids Logical; label centroids in ordination plot.
#' @param ribbon_alpha Alpha transparency for ribbons in stacked plots.
#' @param point_alpha Alpha transparency for points in scatter plots.
#' @param line_size Line width in stacked plots.
#' @param point_size Point size.
#' @param centroid_size Point size for centroids in ordination plots.
#' @param facet_ncol Number of columns in ordination faceting.
#' @param add_smooth Logical; add linear smooth to scatter plot.
#' @param ... Unused, reserved for future extensions.
#'
#' @return A ggplot object or a patchwork object.
#' @export
plot.effect_size_data_nested <- function(
  x,
  type = c("stacked", "scatter", "ordination"),
  component = c("A", "BA"),
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
  validate_effect_size_data(x, variant = "nested")

  type <- match.arg(type)
  component <- match.arg(component)
  inferential_var <- match.arg(inferential_var)
  summary_stat <- match.arg(summary_stat)
  interval <- match.arg(interval)
  colour_by <- match.arg(colour_by)
  scatter_data <- match.arg(scatter_data)
  selection <- match.arg(selection)

  eco_col <- paste0("ecological_effect_", component)
  inf_col <- paste0(inferential_var, "_", component)

  if (type %in% c("stacked", "scatter") && !eco_col %in% names(x)) {
    stop(
      "Column `",
      eco_col,
      "` is not available in the object.",
      call. = FALSE
    )
  }

  if (type %in% c("stacked", "scatter") && !inf_col %in% names(x)) {
    stop(
      "Column `",
      inf_col,
      "` is not available in the object.",
      call. = FALSE
    )
  }

  switch(
    type,
    stacked = .plot_effect_size_nested_stacked(
      x = x,
      component = component,
      inferential_var = inferential_var,
      summary_stat = summary_stat,
      interval = interval,
      ribbon_alpha = ribbon_alpha,
      line_size = line_size,
      point_size = point_size
    ),
    scatter = .plot_effect_size_nested_scatter(
      x = x,
      component = component,
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
