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
#' simResults <- sim_es(data = epiDat, type = "counts", Sest.method = "average",
#'                         cases = 5, N = 100, M = 10,
#'                         n = 5, m = 5, k = 30,
#'                         transformation = "none", method = "bray",
#'                         dummy = FALSE, useParallel = FALSE,
#'                         model = "single.factor",
#'                         jitter.base = 0)
#' }
#' simResults
#'

sim_ES <- function(
  data,
  steps = 14,
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
  # read data and store it in two objects, one for H0 and one for Ha ----
  datH0 <- data
  datH0[, 1] <- as.factor("zero")
  datHa <- data
  datHa[, 1] <- as.factor(data[, 1])
  a <- nlevels(datHa[, 1])
  M <- a

  ## Helper matrix to store labels ----
  NN <- casesHa * k * (n - 1)
  resultsHa <- matrix(NA_real_, nrow = NN, ncol = 6)
  colnames(resultsHa) <- c(
    "dat.sim",
    "k",
    "m",
    "n",
    "FobsHa",
    "SSa"
  )

  resultsHa[, 1] <- rep(seq(casesHa), times = 1, each = (k * (n - 1)))
  resultsHa[, 2] <- rep(1:k, times = (n - 1) * casesHa)
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

  ## Output matrix ----
  resultOut <- matrix(NA_real_, nrow = NN * (steps + 1), ncol = 4)
  colnames(resultOut) <- c("step", "m", "n", "SSa")

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

    # Loop to calculate pseudoF and SS ----

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
        variables = c("balanced_sampling", "dbmanova_oneway", "SS"),
        environment = asNamespace("ecocbo")
      )

      # Executing the loop in parallel
      result1 <- parabar::par_lapply(
        cl,
        x = 1:NN,
        fun = balanced_sampling,
        Y,
        mm,
        nn,
        YPU,
        HaSim,
        HaSim,
        resultsHa,
        transformation,
        method
      ) |>
        unlist() |>
        matrix(ncol = 7, byrow = TRUE)
      colnames(result1) <- c(
        "FobsH0",
        "FobsHa",
        "MSA",
        "MSR",
        "SSf",
        "SSr",
        "SSt"
      )

      # Assigning the results to the outcome matrix
      resultsHa[, c("FobsHa", "SSa")] <- result1[, c("FobsHa", "SSf")]

      parabar::stop_backend(cl)
      rm(result1)
    } else {
      pb <- txtProgressBar(max = NN, style = 3)

      for (i in seq_len(NN)) {
        # Performs the operation iteratively in a for loop
        result1 <- balanced_sampling(
          i,
          Y,
          mm,
          nn,
          YPU,
          HaSim,
          HaSim,
          resultsHa,
          transformation,
          method
        )
        # Assigning results to the results matrix
        resultsHa[i, c("FobsHa", "SSa")] <- result1[, c("FobsHa", "SSf")]

        # Updating the progress bar
        setTxtProgressBar(pb, i)
      }
      rm(result1)
      close(pb)
    }

    idx <- (st * NN + 1):((st + 1) * NN)

    resultOut[idx, 1] <- st
    resultOut[idx, 2:4] <- resultsHa[, c("m", "n", "SSa")]

    cat("Step ", st, ": Done! - ", steps - st, " remaining")
  }

  class(resultOut) <- c("effect_size_data", class(resultOut))

  return(resultOut)
}


#-------------------------------------------
## S3Methods plot()
#-------------------------------------------
#' Plot effect_size_data
#'
#' @method plot effect_size_data
#'
#' @param x An object of class \code{effect_size_data}.
#' @param ... Additional arguments, ignored.
#'
#' @return A ggplot representation of the \code{effect_size_data} object.
#' @export

plot.effect_size_data <- function(x, ...) {
  p1 <- x |>
    as.data.frame() |>
    dplyr::mutate(removed = step / max(step)) |>
    dplyr::group_by(removed) |>
    dplyr::summarise(
      media = mean(SSa, na.rm = TRUE),
      mediana = median(SSa, na.rm = TRUE),
      q25 = quantile(SSa, 0.25, na.rm = TRUE),
      q75 = quantile(SSa, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = removed, y = mediana)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q25, ymax = q75), alpha = 0.3) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(name = "% removed SIMPER contribution") +
    ggplot2::scale_y_continuous(name = "SSa") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(linewidth = 0.4),
      axis.ticks = ggplot2::element_line(linewidth = 0.2)
    )

  print(p1)
  invisible(p1)
}
