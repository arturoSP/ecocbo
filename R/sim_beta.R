#' Calculate Beta Error and Statistical Power from Simulated Samples
#'
#' Estimates the statistical power of a study by comparing variation under null
#' and alternative hypotheses. For instance, if the beta error is 0.25, there is
#' a 25% chance of failing to detect a real difference, and the power of the study
#' is \eqn{1 - \beta}, meaning 0.75 in this case.
#'
#' @param data An object of class `"ecocbo_data"` that results from applying
#' [prep_data()] to a community dataset.
#' @param alpha Numeric. Significance level for Type I error. Defaults to 0.05.
#'
#' @return A list of class "ecocbo_beta", containing:
#'   - `$Power`: a data frame with power and beta estimates across different
#' sampling efforts (`m` sites and `n` samples).
#'   - `$Results`: a data frame with pseudo-F estimates for `simH0` and `simHa`.
#'   - `$alpha`: significance level for Type I error.
#'
#' @details
#' The function displays a summary matrix with estimated power values for various
#' sampling efforts.
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
#' - Anderson, M. J. (2014). Permutational multivariate analysis of
#' variance (PERMANOVA). Wiley statsref: statistics reference online, 1-15.
#' - Guerra‐Castro, E. J., Cajas, J. C., Simões, N., Cruz‐Motta, J.
#'  J., & Mascaró, M. (2021). SSP: an R package to estimate sampling effort in
#'  studies of ecological communities. Ecography, 44(4), 561-573.
#'
#' @seealso
#' [plot_power()]
#' [scompvar()]
#' [sim_cbo()]
#' [prep_data()]
#' [SSP::assempar()]
#' [SSP::simdata()]
#'
#' @aliases simbeta
#'
#' @export
#' @importFrom stats reshape aggregate quantile
#' @importFrom tidyr separate
#' @importFrom dplyr rename
#'
#' @examples
#' sim_beta(data = simResults, alpha = 0.05)
#'

sim_beta <- function(data, alpha = 0.05){

  # Data validation ----
  if(!inherits(data, "ecocbo_data"))
    stop("data is not the right class(\"ecocbo_data\")")
  if(alpha >= 1)
    stop("alpha must be smaller than 1")

  # Preparing an index column combinating m & n, or using n ----
  m_n <- if(data$model == "single.factor"){
    data$Results[,"n"]
  } else {
    as.character(paste(data$Results[,"m"], data$Results[,"n"], sep="_"))
  }

  # Combining results in an intermediate variable ----
  resultsHH <- as.data.frame(cbind(m_n,
                                   pseudoFH0 = data$Results[,"pseudoFH0"],
                                   pseudoFHa = data$Results[,"pseudoFHa"]))
  resultsHH[,"pseudoFH0"] <- as.numeric(resultsHH[,"pseudoFH0"])
  resultsHH[,"pseudoFHa"] <- as.numeric(resultsHH[,"pseudoFHa"])

  # Calculating fCrit and total dimensions ----
  fCrit <- stats::aggregate(resultsHH[,2],
                            by = list(resultsHH[,1]),
                            stats::quantile, probs = (1 - alpha), type= 8, na.rm = T)
  totDim <- stats::aggregate(resultsHH[,2],
                             by = list(resultsHH[,1]),
                             length)

  # Creating power table ----
  powr <- data.frame(m_n = totDim[,1],
                     Power = totDim[,2],
                     Beta = NA,
                     fCrit = fCrit[,2],
                     index = seq(1, dim(fCrit)[1]))

  # Calculating power for each iteration ----
  for(i in powr[,5]){
    partialRes <- resultsHH[resultsHH$m_n == powr[i,1],]
    powr[i,2] <- dim(partialRes[partialRes[,3] >= powr[i,4],])[1] / totDim[i,2]
  }

  # Calculating beta ----
  powr[,3] <- 1 - powr[,2]

  # Adjust the results according to the model ----
  if(data$model == "single.factor"){
    rowidx <- order(powr[,1])

    powr <- dplyr::rename(powr, "n" = "m_n")
    powr <- powr[rowidx, -5]
  } else {
    powr <- powr |>
      tidyr::separate(m_n, into = c("m", "n"), sep = "_", convert = TRUE)

    powr$m <- as.numeric(powr$m)
    powr$n <- as.numeric(powr$n)

    rowidx <- order(powr[,1], powr[,2])

    powr <- powr[rowidx, -6]
  }

  # Creating the resulting object ----
  BetaResult <- list(Power = powr,
                     Results = as.data.frame(data$Results),
                     alpha = alpha,
                     model = data$model,
                     a = data$a)
  class(BetaResult) <- c("ecocbo_beta", class(BetaResult))

  return(BetaResult)
}


#-------------------------------------------
## S3Methods print()
#-------------------------------------------

#' S3Methods for Printing
#'
#' @name print.ecocbo_beta
#'
#' @method print ecocbo_beta
#'
#' @usage
#' \method{print}{ecocbo_beta}(x, ...)
#'
#' @description prints for \code{ecocbo::sim_beta()} objects.
#'
#' @param x Object from \code{ecocbo::sim_beta()} function.
#'
#' @param ... Additional arguments
#'
#' @return Prints the result of \code{ecocbo::sim_beta()} function, showing in an
#' ordered matrix the estimated power for the different experimental designs
#' that were considered.
#'
#' @export
#' @keywords internal

print.ecocbo_beta <- function(x, ...){
  x$Power[, "Power"] <- round(x$Power[, "Power"], 2)

  if (x$model == "single.factor") {
    # Caso: Experimento de un solo factor
    x1 <- x$Power[,c(1:2)]

    x1 <- as.data.frame(x1[,-1])
    rownames(x1) <- paste0("n = ", x$Power$n)
    colnames(x1) <- "Power"

    cat("Power at different sampling efforts (n):\n")
    print(x1)

  } else if (x$model == "nested.symmetric") {
    # Caso: Experimento de dos factores anidados
    x1 <- stats::reshape(
      x$Power[, c("m", "n", "Power")],
      direction = "wide",
      idvar = "m",
      timevar = "n",
      new.row.names = paste0("m = ", c(2:max(x$Power$m)))
    )

    x1 <- x1[, -1]  # Eliminar la columna de 'm'
    colnames(x1) <- paste0("n = ", c(3:max(x$Power$n)))

    cat("Power at different sampling efforts (m x n):\n")
    print(x1)

  } else {
    stop("Unknown model type. Accepted types are 'single.factor' or 'nested.symmetric'.")
  }
  invisible(x1)
}

