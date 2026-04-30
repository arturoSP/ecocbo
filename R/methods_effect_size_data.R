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
  dataH0 = NULL,
  pcoaH0 = NULL,
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
  attr(data, "dataH0") <- dataH0
  attr(data, "pcoaH0") <- pcoaH0
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

#' Remove effect_size_data S3 classes before ordinary data-frame operations
#'
#' @keywords internal
.as_plain_effect_size_df <- function(x) {
  class(x) <- setdiff(
    class(x),
    c("effect_size_data_nested", "effect_size_data")
  )

  x
}

#' @export
as.data.frame.effect_size_data <- function(x, ...) {
  .as_plain_effect_size_df(x)
}

#' @export
as.data.frame.effect_size_data_nested <- function(x, ...) {
  .as_plain_effect_size_df(x)
}

#' @keywords internal
.extract_h0_data <- function(x) {
  dataH0 <- attr(x, "dataH0")

  if (is.null(dataH0)) {
    stop(
      "This object does not contain stored True H0 data in `attr(x, 'dataH0')`.",
      call. = FALSE
    )
  }

  dataH0 <- as.data.frame(dataH0, stringsAsFactors = FALSE)
  dataH0$reduction_level <- NA

  if (nrow(dataH0) == 0) {
    stop(
      "`attr(x, 'dataH0')` is empty.",
      call. = FALSE
    )
  }

  dataH0
}

#' @keywords internal
.es_ok_vec <- function(df) {
  if ("ok" %in% names(df)) {
    return(df$ok %in% TRUE)
  }

  rep(TRUE, nrow(df))
}

#' @keywords internal
.es_check_required_cols <- function(df, cols, object_name = "data") {
  missing_cols <- setdiff(cols, names(df))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in `",
      object_name,
      "`: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
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
        "pseudoF_BA"
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
    cat("Resampling (k):", length(unique(df$k)), "\n")
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
.es_base_theme <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(linewidth = 0.4),
      axis.ticks = ggplot2::element_line(linewidth = 0.2)
    )
}

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
      x = "Reduction level",
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
      x = "Reduction level",
      y = NULL
    ) +
    base_theme
}

#' @keywords internal
.plot_effect_size_scatter <- function(
  x,
  inferential_var = c("omega2", "R2", "pseudoF"),
  color_by = c("reduction_level", "none"),
  scatter_data = c("raw", "summary"),
  summary_stat = c("median", "mean"),
  point_alpha = 0.7,
  point_size = 1.8,
  add_smooth = FALSE
) {
  inferential_var <- match.arg(inferential_var)
  color_by <- match.arg(color_by)
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

  mapping <- if (color_by == "reduction_level") {
    ggplot2::aes(
      x = ecological_effect,
      y = inferential_value,
      color = reduction_level
    )
  } else {
    ggplot2::aes(
      x = ecological_effect,
      y = inferential_value
    )
  }

  p <- ggplot2::ggplot(df, mapping) +
    ggplot2::geom_point(alpha = point_alpha, size = point_size) +
    ggplot2::labs(
      x = "Ecological effect size",
      y = y_lab_inf,
      color = "Reduction level"
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
    effect_col = "ecological_effect",
    group_var = "group",
    color_var = group_var,
    label_var = group_var,
    legend_title = "Group",
    show_centroids = TRUE,
    label_centroids = TRUE,
    point_alpha = 0.7,
    point_size = 1.8,
    centroid_size = 2.8,
    facet_ncol = NULL
) {
  selection <- match.arg(selection)

  pcoa_df <- .extract_pcoa_data(x)

  required_pcoa_cols <- unique(c(group_var, color_var, label_var, "Axis1", "Axis2"))
  .es_check_required_cols(pcoa_df, required_pcoa_cols, object_name = "stored PCoA data")

  plot_df <- .select_representative_pcoa(
    x = x,
    pcoa_df = pcoa_df,
    effect_col = effect_col,
    reduction_levels = reduction_levels,
    selection = selection
  )

  plot_df$.color <- plot_df[[color_var]]

  centroids <- .compute_pcoa_centroids(
    plot_df,
    group_var = group_var,
    color_var = color_var,
    label_var = label_var
  )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = Axis1, y = Axis2, color = .color)
  ) +
    ggplot2::geom_point(alpha = point_alpha, size = point_size)

  if (show_centroids) {
    p <- p +
      ggplot2::geom_point(
        data = centroids,
        ggplot2::aes(x = Axis1, y = Axis2, color = .color),
        inherit.aes = FALSE,
        size = centroid_size,
        shape = 4
      )
  }

  if (label_centroids) {
    p <- p +
      ggplot2::geom_text(
        data = centroids,
        ggplot2::aes(x = Axis1, y = Axis2, label = .label, color = .color),
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
      color = legend_title
    ) +
    .es_base_theme()
}

#' @keywords internal
.plot_effect_size_nested_ordination <- function(
    x,
    component = c("A", "BA"),
    reduction_levels = NULL,
    selection = c("median", "first", "min", "max"),
    show_centroids = TRUE,
    label_centroids = TRUE,
    point_alpha = 0.7,
    point_size = 1.8,
    centroid_size = 2.8,
    facet_ncol = NULL
) {
  component <- match.arg(component)

  if (component == "A") {
    return(
      .plot_effect_size_ordination(
        x = x,
        reduction_levels = reduction_levels,
        selection = selection,
        effect_col = "ecological_effect_A",
        group_var = "sector",
        color_var = "sector",
        label_var = "sector",
        legend_title = "Sector",
        show_centroids = show_centroids,
        label_centroids = label_centroids,
        point_alpha = point_alpha,
        point_size = point_size,
        centroid_size = centroid_size,
        facet_ncol = facet_ncol
      )
    )
  }

  .plot_effect_size_ordination(
    x = x,
    reduction_levels = reduction_levels,
    selection = selection,
    effect_col = "ecological_effect_BA",
    group_var = "sector_site",
    color_var = "sector",
    label_var = "sector_site",
    legend_title = "Sector",
    show_centroids = show_centroids,
    label_centroids = label_centroids,
    point_alpha = point_alpha,
    point_size = point_size,
    centroid_size = centroid_size,
    facet_ncol = facet_ncol
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
      x = "Reduction level",
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
      )
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
      x = "Reduction level",
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
#' @param color_by Whether to color scatter points by
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
    type = c("stacked", "scatter", "ordination", "true_h0", "ES_H0"),
    inferential_var = c("omega2", "R2", "pseudoF"),
    summary_stat = c("median", "mean"),
    interval = c("iqr", "sd"),
    color_by = c("reduction_level", "none"),
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
    h0_alpha = 0.05,
    h0_x_var = c("pseudoF", "ecological_effect", "omega2", "R2"),
    h0_color_by = c("reduction_level", "step", "none"),
    facet_effort = FALSE,
    h0_quantile_type = 7,
    ...
) {
  validate_effect_size_data(x)

  type <- match.arg(type)
  if (type == "ES_H0") {
    type <- "true_h0"
  }

  inferential_var <- match.arg(inferential_var)
  summary_stat <- match.arg(summary_stat)
  interval <- match.arg(interval)
  color_by <- match.arg(color_by)
  scatter_data <- match.arg(scatter_data)
  selection <- match.arg(selection)
  h0_x_var <- match.arg(h0_x_var)
  h0_color_by <- match.arg(h0_color_by)

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
      color_by = color_by,
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
    ),
    true_h0 = .plot_effect_size_true_h0(
      x = x,
      h0_alpha = h0_alpha,
      h0_x_var = h0_x_var,
      h0_color_by = h0_color_by,
      facet_effort = facet_effort,
      facet_ncol = facet_ncol,
      quantile_type = h0_quantile_type,
      point_alpha = point_alpha,
      point_size = point_size
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
  color_by = c("reduction_level", "none"),
  scatter_data = c("summary", "raw"),
  summary_stat = c("median", "mean"),
  point_alpha = 0.7,
  point_size = 1.8,
  add_smooth = FALSE
) {
  component <- match.arg(component)
  inferential_var <- match.arg(inferential_var)
  color_by <- match.arg(color_by)
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
      color = if (color_by == "reduction_level") reduction_level else NULL
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
      color = if (color_by == "reduction_level") "Reduction level" else NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(linewidth = 0.4),
      axis.ticks = ggplot2::element_line(linewidth = 0.2)
    )

  if (color_by == "none") {
    p <- p + ggplot2::guides(color = "none")
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
#' @param color_by Whether to color scatter points by
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
    type = c("stacked", "scatter", "ordination", "true_h0", "ES_H0"),
    component = c("A", "BA"),
    inferential_var = c("omega2", "R2", "pseudoF"),
    summary_stat = c("median", "mean"),
    interval = c("iqr", "sd"),
    color_by = c("reduction_level", "none"),
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
    h0_alpha = 0.05,
    h0_x_var = c("pseudoF", "ecological_effect", "omega2", "R2"),
    h0_color_by = c("reduction_level", "step", "none"),
    facet_effort = FALSE,
    facet_subset = c("reduced", "all", "all_n", "all_m"),
    facet_start_at = 2,
    facet_m_value = NULL,
    facet_n_value = NULL,
    h0_quantile_type = 7,
    ...
) {
  validate_effect_size_data(x, variant = "nested")

  type <- match.arg(type)
  if (type == "ES_H0") {
    type <- "true_h0"
  }

  component <- match.arg(component)
  inferential_var <- match.arg(inferential_var)
  summary_stat <- match.arg(summary_stat)
  interval <- match.arg(interval)
  color_by <- match.arg(color_by)
  scatter_data <- match.arg(scatter_data)
  selection <- match.arg(selection)
  h0_x_var <- match.arg(h0_x_var)
  h0_color_by <- match.arg(h0_color_by)
  facet_subset <- match.arg(facet_subset)
  facet_subset <- match.arg(facet_subset)

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
      color_by = color_by,
      scatter_data = scatter_data,
      summary_stat = summary_stat,
      point_alpha = point_alpha,
      point_size = point_size,
      add_smooth = add_smooth
    ),
    ordination = .plot_effect_size_nested_ordination(
      x = x,
      component = component,
      reduction_levels = reduction_levels,
      selection = selection,
      show_centroids = show_centroids,
      label_centroids = label_centroids,
      point_alpha = point_alpha,
      point_size = point_size,
      centroid_size = centroid_size,
      facet_ncol = facet_ncol
    ),
    true_h0 = .plot_effect_size_nested_true_h0(
      x = x,
      component = component,
      h0_alpha = h0_alpha,
      h0_x_var = h0_x_var,
      h0_color_by = h0_color_by,
      facet_effort = facet_effort,
      facet_subset = facet_subset,
      facet_start_at = facet_start_at,
      facet_m_value = facet_m_value,
      facet_n_value = facet_n_value,
      facet_ncol = facet_ncol,
      quantile_type = h0_quantile_type,
      point_alpha = point_alpha,
      point_size = point_size
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
    effect_col = "ecological_effect",
    reduction_levels = NULL,
    selection = c("median", "first", "min", "max")
) {
  selection <- match.arg(selection)

  df <- as.data.frame(x)

  if (!effect_col %in% names(df)) {
    stop(
      "Column `", effect_col, "` is not available in the object.",
      call. = FALSE
    )
  }

  if (!is.null(reduction_levels)) {
    df <- df[df$reduction_level %in% reduction_levels, , drop = FALSE]
    pcoa_df <- pcoa_df[
      pcoa_df$reduction_level %in% reduction_levels,
      ,
      drop = FALSE
    ]
  }

  id_cols <- c("reduction_level", "step", "dat_sim", "k", "m", "n")
  id_cols <- intersect(id_cols, intersect(names(df), names(pcoa_df)))

  selected <- switch(
    selection,
    median = {
      targets <- df |>
        dplyr::group_by(reduction_level) |>
        dplyr::summarise(
          target_eco = stats::median(.data[[effect_col]], na.rm = TRUE),
          .groups = "drop"
        )

      df |>
        dplyr::inner_join(targets, by = "reduction_level") |>
        dplyr::mutate(
          .dist_to_target = abs(.data[[effect_col]] - target_eco)
        ) |>
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
        dplyr::slice_min(.data[[effect_col]], n = 1, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::select(dplyr::all_of(id_cols))
    },
    max = {
      df |>
        dplyr::group_by(reduction_level) |>
        dplyr::slice_max(.data[[effect_col]], n = 1, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::select(dplyr::all_of(id_cols))
    }
  )

  out <- dplyr::inner_join(pcoa_df, selected, by = id_cols)

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
.compute_pcoa_centroids <- function(
    df,
    group_var = "group",
    color_var = group_var,
    label_var = group_var
) {
  vars_to_keep <- unique(c(group_var, color_var, label_var))

  out <- df |>
    dplyr::group_by(
      panel_label,
      reduction_level,
      dplyr::across(dplyr::all_of(vars_to_keep))
    ) |>
    dplyr::summarise(
      Axis1 = mean(Axis1, na.rm = TRUE),
      Axis2 = mean(Axis2, na.rm = TRUE),
      .groups = "drop"
    )

  out$.group <- out[[group_var]]
  out$.color <- out[[color_var]]
  out$.label <- out[[label_var]]

  out
}

# Helper para resumen de H0
#' Summarise empirical rejection rates using the explicit True H0 scenario
#'
#' @keywords internal
.summarise_true_h0_rejection_single <- function(
    x,
    h0_alpha = 0.05,
    quantile_type = 7,
    include_true_h0 = TRUE
) {
  df <- as.data.frame(x)
  dataH0 <- .extract_h0_data(x)

  .es_check_required_cols(
    df,
    cols = c("m", "n", "pseudoF", "ecological_effect"),
    object_name = "x"
  )

  .es_check_required_cols(
    dataH0,
    cols = c("m", "n", "pseudoF"),
    object_name = "attr(x, 'dataH0')"
  )

  df$.ok <- .es_ok_vec(df)
  dataH0$.ok <- .es_ok_vec(dataH0)

  if (!"scenario" %in% names(df)) {
    df$scenario <- "effect path"
  }

  if (!"step" %in% names(df)) {
    df$step <- NA_integer_
  }

  if (!"reduction_level" %in% names(df)) {
    df$reduction_level <- NA_real_
  }

  crit_by_effort <- dataH0 |>
    dplyr::filter(
      .ok,
      is.finite(pseudoF)
    ) |>
    dplyr::group_by(m, n) |>
    dplyr::summarise(
      n_H0 = dplyr::n(),
      pseudoFcrit = unname(
        stats::quantile(
          pseudoF,
          probs = 1 - h0_alpha,
          na.rm = TRUE,
          type = quantile_type
        )
      ),
      .groups = "drop"
    )

  metric_cols <- intersect(
    c("ecological_effect", "omega2", "R2", "pseudoF"),
    names(df)
  )

  path_summary <- df |>
    dplyr::left_join(
      crit_by_effort,
      by = c("m", "n")
    ) |>
    dplyr::mutate(
      reject = dplyr::case_when(
        !.ok ~ NA,
        !is.finite(pseudoF) | is.na(pseudoFcrit) ~ NA,
        TRUE ~ pseudoF > pseudoFcrit
      )
    ) |>
    dplyr::filter(!is.na(reject)) |>
    dplyr::group_by(
      scenario,
      step,
      reduction_level,
      m,
      n
    ) |>
    dplyr::summarise(
      rejection_rate = mean(reject, na.rm = TRUE),
      n_valid = sum(!is.na(reject)),
      n_H0 = dplyr::first(n_H0),
      pseudoFcrit = dplyr::first(pseudoFcrit),
      dplyr::across(
        dplyr::all_of(metric_cols),
        ~ mean(.x, na.rm = TRUE),
        .names = "mean_{.col}"
      ),
      .groups = "drop"
    )

  if (!include_true_h0) {
    return(path_summary)
  }

  h0_metric_cols <- intersect(
    c("ecological_effect", "omega2", "R2", "pseudoF"),
    names(dataH0)
  )

  h0_summary <- dataH0 |>
    dplyr::left_join(
      crit_by_effort,
      by = c("m", "n")
    ) |>
    dplyr::mutate(
      scenario = "True H0",
      step = NA_integer_,
      reduction_level = NA_real_,
      reject = dplyr::case_when(
        !.ok ~ NA,
        !is.finite(pseudoF) | is.na(pseudoFcrit) ~ NA,
        TRUE ~ pseudoF > pseudoFcrit
      )
    ) |>
    dplyr::filter(!is.na(reject)) |>
    dplyr::group_by(
      scenario,
      step,
      reduction_level,
      m,
      n
    ) |>
    dplyr::summarise(
      rejection_rate = mean(reject, na.rm = TRUE),
      n_valid = sum(!is.na(reject)),
      n_H0 = dplyr::first(n_H0),
      pseudoFcrit = dplyr::first(pseudoFcrit),
      dplyr::across(
        dplyr::all_of(h0_metric_cols),
        ~ mean(.x, na.rm = TRUE),
        .names = "mean_{.col}"
      ),
      .groups = "drop"
    )

  out <- dplyr::bind_rows(path_summary, h0_summary)

  out$scenario <- factor(
    out$scenario,
    levels = c("attenuated", "True Ha", "True H0", "effect path")
  )

  out
}

#' @keywords internal
.plot_effect_size_true_h0 <- function(
    x,
    h0_alpha = 0.05,
    h0_x_var = c("pseudoF", "ecological_effect", "omega2", "R2"),
    h0_color_by = c("reduction_level", "step", "none"),
    facet_effort = FALSE,
    facet_ncol = NULL,
    quantile_type = 7,
    point_alpha = 0.7,
    point_size = 1.8
) {
  h0_x_var <- match.arg(h0_x_var)
  h0_color_by <- match.arg(h0_color_by)

  df_sum <- .summarise_true_h0_rejection_single(
    x = x,
    h0_alpha = h0_alpha,
    quantile_type = quantile_type
  )

  x_col <- paste0("mean_", h0_x_var)

  if (!x_col %in% names(df_sum)) {
    stop(
      "Column `",
      x_col,
      "` is not available after summarising the object.",
      call. = FALSE
    )
  }

  x_lab <- switch(
    h0_x_var,
    pseudoF = "Mean pseudo-F",
    ecological_effect = "Mean ecological effect size",
    omega2 = expression("Mean " * omega^2),
    R2 = expression("Mean " * R^2)
  )

  color_lab <- switch(
    h0_color_by,
    reduction_level = "Reduction level",
    step = "Step",
    none = NULL
  )

  if (facet_effort) {
    df_sum <- .order_effort_facets(df_sum)
  }

  if (h0_color_by == "none") {
    p <- ggplot2::ggplot(
      df_sum,
      ggplot2::aes(
        x = .data[[x_col]],
        y = rejection_rate,
        shape = scenario
      )
    )
  } else {
    p <- ggplot2::ggplot(
      df_sum,
      ggplot2::aes(
        x = .data[[x_col]],
        y = rejection_rate,
        color = .data[[h0_color_by]],
        shape = scenario
      )
    )
  }

  p <- p +
    ggplot2::geom_hline(
      yintercept = h0_alpha,
      linetype = "dashed",
      linewidth = 0.3
    ) +
    ggplot2::geom_point(
      alpha = point_alpha,
      size = point_size
    ) +
    ggplot2::scale_shape_manual(
      values = .es_scenario_shapes(),
      drop = FALSE
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = x_lab,
      y = "Empirical rejection rate",
      color = color_lab,
      shape = "Scenario"
    ) +
    .es_base_theme()

  if (h0_color_by == "none") {
    p <- p + ggplot2::guides(color = "none")
  }

  if (facet_effort) {
    p <- p +
      ggplot2::facet_wrap(
        ~effort_label,
        ncol = facet_ncol
      )
  }

  p
}

#' Summarise empirical rejection rates for nested effect-size objects
#'
#' @keywords internal
.summarise_true_h0_rejection_nested <- function(
    x,
    component = c("A", "BA"),
    h0_alpha = 0.05,
    quantile_type = 7,
    include_true_h0 = TRUE
) {
  component <- match.arg(component)

  df <- as.data.frame(x)
  dataH0 <- .extract_h0_data(x)

  pseudoF_col <- paste0("pseudoF_", component)
  eco_col <- paste0("ecological_effect_", component)
  omega_col <- paste0("omega2_", component)
  r2_col <- paste0("R2_", component)

  .es_check_required_cols(
    df,
    cols = c("m", "n", pseudoF_col, eco_col),
    object_name = "x"
  )

  .es_check_required_cols(
    dataH0,
    cols = c("m", "n", pseudoF_col),
    object_name = "attr(x, 'dataH0')"
  )

  df$.ok <- .es_ok_vec(df)
  dataH0$.ok <- .es_ok_vec(dataH0)

  if (!"scenario" %in% names(df)) {
    df$scenario <- "effect path"
  }

  if (!"step" %in% names(df)) {
    df$step <- NA_integer_
  }

  if (!"reduction_level" %in% names(df)) {
    df$reduction_level <- NA_real_
  }

  crit_by_effort <- dataH0 |>
    dplyr::filter(
      .ok,
      is.finite(.data[[pseudoF_col]])
    ) |>
    dplyr::group_by(m, n) |>
    dplyr::summarise(
      n_H0 = dplyr::n(),
      pseudoFcrit = unname(
        stats::quantile(
          .data[[pseudoF_col]],
          probs = 1 - h0_alpha,
          na.rm = TRUE,
          type = quantile_type
        )
      ),
      .groups = "drop"
    )

  metric_cols <- intersect(
    c(eco_col, omega_col, r2_col, pseudoF_col),
    names(df)
  )

  path_summary <- df |>
    dplyr::left_join(
      crit_by_effort,
      by = c("m", "n")
    ) |>
    dplyr::mutate(
      reject = dplyr::case_when(
        !.ok ~ NA,
        !is.finite(.data[[pseudoF_col]]) | is.na(pseudoFcrit) ~ NA,
        TRUE ~ .data[[pseudoF_col]] > pseudoFcrit
      )
    ) |>
    dplyr::filter(!is.na(reject)) |>
    dplyr::group_by(
      scenario,
      step,
      reduction_level,
      m,
      n
    ) |>
    dplyr::summarise(
      rejection_rate = mean(reject, na.rm = TRUE),
      n_valid = sum(!is.na(reject)),
      n_H0 = dplyr::first(n_H0),
      pseudoFcrit = dplyr::first(pseudoFcrit),
      dplyr::across(
        dplyr::all_of(metric_cols),
        ~ mean(.x, na.rm = TRUE),
        .names = "mean_{.col}"
      ),
      .groups = "drop"
    )

  if (!include_true_h0) {
    return(path_summary)
  }

  h0_metric_cols <- intersect(
    c(eco_col, omega_col, r2_col, pseudoF_col),
    names(dataH0)
  )

  h0_summary <- dataH0 |>
    dplyr::left_join(
      crit_by_effort,
      by = c("m", "n")
    ) |>
    dplyr::mutate(
      scenario = "True H0",
      step = NA_integer_,
      reduction_level = NA_real_,
      reject = dplyr::case_when(
        !.ok ~ NA,
        !is.finite(.data[[pseudoF_col]]) | is.na(pseudoFcrit) ~ NA,
        TRUE ~ .data[[pseudoF_col]] > pseudoFcrit
      )
    ) |>
    dplyr::filter(!is.na(reject)) |>
    dplyr::group_by(
      scenario,
      step,
      reduction_level,
      m,
      n
    ) |>
    dplyr::summarise(
      rejection_rate = mean(reject, na.rm = TRUE),
      n_valid = sum(!is.na(reject)),
      n_H0 = dplyr::first(n_H0),
      pseudoFcrit = dplyr::first(pseudoFcrit),
      dplyr::across(
        dplyr::all_of(h0_metric_cols),
        ~ mean(.x, na.rm = TRUE),
        .names = "mean_{.col}"
      ),
      .groups = "drop"
    )

  out <- dplyr::bind_rows(path_summary, h0_summary)

  out$scenario <- factor(
    out$scenario,
    levels = c("attenuated", "True Ha", "True H0", "effect path")
  )

  out
}

#' @keywords internal
.plot_effect_size_nested_true_h0 <- function(
    x,
    component = c("A", "BA"),
    h0_alpha = 0.05,
    h0_x_var = c("pseudoF", "ecological_effect", "omega2", "R2"),
    h0_color_by = c("reduction_level", "step", "none"),
    facet_effort = FALSE,
    facet_subset = c("reduced", "all", "all_n", "all_m"),
    facet_start_at = 2,
    facet_m_value = NULL,
    facet_n_value = NULL,
    facet_ncol = NULL,
    quantile_type = 7,
    point_alpha = 0.7,
    point_size = 1.8
) {
  component <- match.arg(component)
  h0_x_var <- match.arg(h0_x_var)
  h0_color_by <- match.arg(h0_color_by)
  facet_subset <- match.arg(facet_subset)
  facet_subset <- match.arg(facet_subset)

  df_sum <- .summarise_true_h0_rejection_nested(
    x = x,
    component = component,
    h0_alpha = h0_alpha,
    quantile_type = quantile_type
  )

  raw_x_col <- if (h0_x_var == "pseudoF") {
    paste0("pseudoF_", component)
  } else if (h0_x_var == "ecological_effect") {
    paste0("ecological_effect_", component)
  } else {
    paste0(h0_x_var, "_", component)
  }

  x_col <- paste0("mean_", raw_x_col)

  if (!x_col %in% names(df_sum)) {
    stop(
      "Column `",
      x_col,
      "` is not available after summarising the object.",
      call. = FALSE
    )
  }

  component_lab <- if (component == "A") "A" else "B(A)"

  x_lab <- switch(
    h0_x_var,
    pseudoF = paste0("Mean pseudo-F (", component_lab, ")"),
    ecological_effect = paste0("Mean ecological effect size (", component_lab, ")"),
    omega2 = if (component == "A") {
      expression("Mean " * omega[A]^2)
    } else {
      expression("Mean " * omega[B(A)]^2)
    },
    R2 = if (component == "A") {
      expression("Mean " * R[A]^2)
    } else {
      expression("Mean " * R[B(A)]^2)
    }
  )

  color_lab <- switch(
    h0_color_by,
    reduction_level = "Reduction level",
    step = "Step",
    none = NULL
  )

  if (facet_effort) {
    df_sum <- .order_effort_facets(
      df_sum,
      subset = facet_subset,
      start_at = facet_start_at,
      facet_m_value = facet_m_value,
      facet_n_value = facet_n_value
    )
  }

  if (h0_color_by == "none") {
    p <- ggplot2::ggplot(
      df_sum,
      ggplot2::aes(
        x = .data[[x_col]],
        y = rejection_rate,
        shape = scenario
      )
    )
  } else {
    p <- ggplot2::ggplot(
      df_sum,
      ggplot2::aes(
        x = .data[[x_col]],
        y = rejection_rate,
        color = .data[[h0_color_by]],
        shape = scenario
      )
    )
  }

  p <- p +
    ggplot2::geom_hline(
      yintercept = h0_alpha,
      linetype = "dashed",
      linewidth = 0.3
    ) +
    ggplot2::geom_point(
      alpha = point_alpha,
      size = point_size
    ) +
    ggplot2::scale_shape_manual(
      values = .es_scenario_shapes(),
      drop = FALSE
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = x_lab,
      y = paste0("Empirical rejection rate (", component_lab, ")"),
      color = color_lab,
      shape = "Scenario"
    ) +
    .es_base_theme()

  if (h0_color_by == "none") {
    p <- p + ggplot2::guides(color = "none")
  }

  if (facet_effort) {
    p <- p +
      ggplot2::facet_wrap(
        ~effort_label,
        ncol = facet_ncol
      )
  }

  p
}

#' Order sampling-effort facet labels numerically and optionally subset them
#'
#' subset modes:
#' - "all": keep all (m, n) combinations
#' - "reduced": keep diagonal + tail path
#' - "all_n": keep all n values for one selected m
#' - "all_m": keep all m values for one selected n
#'
#' @keywords internal
.order_effort_facets <- function(
    df,
    subset = c("all", "reduced", "all_n", "all_m"),
    start_at = 2,
    facet_m_value = NULL,
    facet_n_value = NULL
) {
  subset <- match.arg(subset)

  df$.m_num <- suppressWarnings(as.numeric(as.character(df$m)))
  df$.n_num <- suppressWarnings(as.numeric(as.character(df$n)))

  effort_df <- df |>
    dplyr::distinct(m, n, .m_num, .n_num) |>
    dplyr::filter(!is.na(.m_num), !is.na(.n_num)) |>
    dplyr::arrange(.m_num, .n_num)

  if (nrow(effort_df) == 0) {
    df$effort_label <- factor(character(0))
    df$.m_num <- NULL
    df$.n_num <- NULL
    return(df)
  }

  if (subset == "reduced") {
    keep_pairs <- .select_nested_effort_subset(df, start_at = start_at)

    df <- df |>
      dplyr::semi_join(keep_pairs, by = c("m", "n"))

  } else if (subset == "all_n") {
    if (is.null(facet_m_value)) {
      facet_m_value <- max(effort_df$.m_num, na.rm = TRUE)
    }

    df <- df |>
      dplyr::filter(.m_num == facet_m_value)

  } else if (subset == "all_m") {
    if (is.null(facet_n_value)) {
      facet_n_value <- max(effort_df$.n_num, na.rm = TRUE)
    }

    df <- df |>
      dplyr::filter(.n_num == facet_n_value)
  }

  df <- df |>
    dplyr::arrange(.m_num, .n_num, m, n)

  effort_levels <- df |>
    dplyr::distinct(m, n, .m_num, .n_num) |>
    dplyr::arrange(.m_num, .n_num, m, n) |>
    dplyr::transmute(label = paste0("m = ", m, ", n = ", n)) |>
    dplyr::pull(label)

  df$effort_label <- factor(
    paste0("m = ", df$m, ", n = ", df$n),
    levels = unique(effort_levels)
  )

  df$.m_num <- NULL
  df$.n_num <- NULL

  df
}

#' Shape values for effect-size scenarios
#'
#' @keywords internal
.es_scenario_shapes <- function() {
  c(
    "attenuated" = 16,
    "True Ha" = 17,
    "True H0" = 15,
    "effect path" = 16
  )
}
#' Select a reduced nested sampling-effort path for faceting
#'
#' Keeps diagonal pairs (2,2), (3,3), ..., up to min(max(m), max(n)),
#' and then extends along the remaining dimension using the maximum value
#' of the other one. For example, if max(m)=10 and max(n)=15, the selected
#' path is:
#' (2,2), (3,3), ..., (10,10), (10,11), ..., (10,15)
#'
#' @keywords internal
.select_nested_effort_subset <- function(df, start_at = 2) {
  effort_df <- df |>
    dplyr::distinct(m, n) |>
    dplyr::mutate(
      .m_num = suppressWarnings(as.numeric(as.character(m))),
      .n_num = suppressWarnings(as.numeric(as.character(n)))
    ) |>
    dplyr::filter(!is.na(.m_num), !is.na(.n_num)) |>
    dplyr::arrange(.m_num, .n_num)

  if (nrow(effort_df) == 0) {
    return(effort_df[, c("m", "n"), drop = FALSE])
  }

  max_m <- max(effort_df$.m_num, na.rm = TRUE)
  max_n <- max(effort_df$.n_num, na.rm = TRUE)
  diag_end <- min(max_m, max_n)

  start_at <- max(start_at, min(effort_df$.m_num, effort_df$.n_num, na.rm = TRUE))

  diag_pairs <- tibble::tibble(
    .m_num = seq.int(start_at, diag_end),
    .n_num = seq.int(start_at, diag_end)
  )

  tail_pairs <- if (max_n > diag_end) {
    tibble::tibble(
      .m_num = rep(max_m, max_n - diag_end),
      .n_num = seq.int(diag_end + 1, max_n)
    )
  } else if (max_m > diag_end) {
    tibble::tibble(
      .m_num = seq.int(diag_end + 1, max_m),
      .n_num = rep(max_n, max_m - diag_end)
    )
  } else {
    tibble::tibble(
      .m_num = numeric(0),
      .n_num = numeric(0)
    )
  }

  selected_num <- dplyr::bind_rows(diag_pairs, tail_pairs)

  effort_df |>
    dplyr::semi_join(selected_num, by = c(".m_num", ".n_num")) |>
    dplyr::distinct(m, n, .m_num, .n_num) |>
    dplyr::arrange(.m_num, .n_num) |>
    dplyr::select(m, n)
}
