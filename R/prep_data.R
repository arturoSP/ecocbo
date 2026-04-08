#' Prepare Data for Evaluation
#'
#' Formats and arranges the initial data so that it can be
#' readily used by the other functions in the package. The function first gets
#' the species names and the number of samples for each species from the input
#' data frame. Then, it permutes the sampling efforts and calculates the pseudo-F
#' statistic and the mean squares for each permutation. Finally, it returns a
#' data frame with the permutations, pseudo-F statistic, and mean squares.
#'
#' @param data Data frame where columns represent species names and rows correspond
#' to samples.
#'   - For `"single.factor"` analysis: The first column should indicate the replicate
#'   to which the sample belongs.
#'   - For `"nested.symmetric"` analysis: The first column should indicate the
#'   treatment, and the second column should indicate the replicate.
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
#' [sim_beta()]
#' [plot_power()]
#' [sim_cbo()]
#' [scompvar()]
#'
#' @aliases prepdata
#'
#' @export
#' @importFrom SSP assempar simdata
#'
#' @examples
#' \donttest{
#' simResults <- prep_data(data = epiDat, type = "counts", Sest.method = "average",
#'                         cases = 5, N = 100, M = 10,
#'                         n = 5, m = 5, k = 30,
#'                         transformation = "none", method = "bray",
#'                         dummy = FALSE, useParallel = FALSE,
#'                         model = "single.factor",
#'                         jitter.base = 0)
#' }
#' simResults
#'

prep_data <- function(
  data,
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
  # Check the inputs ----

  if (n > N) {
    stop("'n' must be equal or less than 'N' on simulated data")
  }
  if (ceiling(n) != floor(n)) {
    stop("n must be integer")
  }
  if (n <= 1) {
    stop("n must be larger than 1")
  }

  if (model != "single.factor") {
    if (m > M) {
      stop("'m' must be equal or less than 'M' on simulated data")
    }
    if (ceiling(m) != floor(m)) {
      stop("m must be integer")
    }
    if (m <= 1) {
      stop("m must be larger than 1")
    }
  }

  # The function to work with depends on the selected model
  Results <- if (model == "single.factor") {
    prep_data_single(
      data,
      type,
      Sest.method,
      cases,
      N,
      M,
      n,
      m,
      k,
      transformation,
      method,
      dummy,
      useParallel,
      jitter.base
    )
  } else if (model == "nested.symmetric") {
    prep_data_nestedsymmetric(
      data,
      type,
      Sest.method,
      cases,
      N,
      M,
      n,
      m,
      k,
      transformation,
      method,
      dummy,
      useParallel,
      model,
      jitter.base
    )
  }

  return(Results)
}


#-------------------------------------------
## S3Methods print()
#-------------------------------------------

#' Coerce ecocbo_data to data.frame
#'
#' @name as.data.frame.ecocbo_data
#' @method as.data.frame ecocbo_data
#'
#' @param x An object of class \code{ecocbo_data}.
#' @param row.names Passed to \code{as.data.frame()}.
#' @param optional Passed to \code{as.data.frame()}.
#' @param ... Additional arguments, ignored.
#'
#' @return A data.frame representation of the \code{Results} component.
#' @export
#' @keywords internal
as.data.frame.ecocbo_data <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  ...
) {
  if (is.null(x$Results)) {
    stop("`x` does not contain a `Results` component.")
  }

  as.data.frame(x$Results, row.names = row.names, optional = optional)
}

#' Print method for ecocbo_data
#'
#' @name print.ecocbo_data
#' @method print ecocbo_data
#'
#' @param x An object of class \code{ecocbo_data}.
#' @param ... Additional arguments, ignored.
#'
#' @return The input object, invisibly.
#' @export
#' @keywords internal
print.ecocbo_data <- function(x, ...) {
  dims <- if (!is.null(x$Results)) {
    dim(x$Results)
  } else {
    c(NA_integer_, NA_integer_)
  }

  cat("<ecocbo_data>\n")
  cat("  model :", x$model, "\n")
  cat("  a     :", x$a, "\n")
  cat("  Results:", dims[1], "rows x", dims[2], "columns\n")

  invisible(x)
}

#' Summary method for ecocbo_data
#'
#' @name summary.ecocbo_data
#' @method summary ecocbo_data
#'
#' @param object An object of class \code{ecocbo_data}.
#' @param ... Additional arguments, ignored.
#'
#' @return A summary list of the main object components.
#' @export
#' @keywords internal
summary.ecocbo_data <- function(object, ...) {
  out <- list(
    model = object$model,
    a = object$a,
    dimensions = if (!is.null(object$Results)) dim(object$Results) else NULL,
    colnames = if (!is.null(object$Results)) colnames(object$Results) else NULL
  )
  class(out) <- "summary.ecocbo_data"
  out
}

#' @export
print.summary.ecocbo_data <- function(x, ...) {
  cat("Summary of <ecocbo_data>\n")
  cat("  model      :", x$model, "\n")
  cat("  a          :", x$a, "\n")
  if (!is.null(x$dimensions)) {
    cat("  dimensions :", x$dimensions[1], "x", x$dimensions[2], "\n")
  }
  if (!is.null(x$colnames)) {
    cat("  variables   :", paste(x$colnames, collapse = ", "), "\n")
  }
  invisible(x)
}
