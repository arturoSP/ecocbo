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
#' @return \code{sim_es()} returns an object of class "effect_size_data".
#'
#' An object of class "effect_size_data" is a list containing:
#'   - \code{$resultOut}
#'   - \code{$pcoaOut}
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
#' ES1 <- sim_ES(data = epiDat, steps = 10, type = "counts", Sest.method = "average",
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
  n = 30,
  m = NULL,
  k = 50,
  transformation = "none",
  method = "bray",
  dummy = TRUE,
  useParallel = TRUE,
  model = match.arg(model, c("single.factor", "nested.symmetric")),
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
    sim_es_single(
      data,
      steps,
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
  } else if (model == "nested.symmetric") {
    sim_es_nested(
      data,
      steps,
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
