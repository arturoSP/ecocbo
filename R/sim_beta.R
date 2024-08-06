#' Calculate beta and power out of simulated samples
#'
#' \code{sim_beta()} can be used to assess the power of a study by comparing the
#' variation when one can assume whether an ecological community does not have
#' composition differences (H0 true) or it does (H0 false). For example, if the
#' beta error is 0.25, then there is a 25% chance of failing to detect a
#' difference even if the difference is real. The power of the study is
#' \eqn{1 - \beta}, so in this example, the power of the study is 0.75.
#'
#' @param data An object of class "ecocbo_data" that results from applying
#' [prep_data()] to a community data frame.
#' @param alpha Level of significance for Type I error. Defaults to 0.05.
#'
#' @return \code{sim_data()} returns an object of class "ecocbo_beta".
#'
#' The function \code{print()} is used to present a matrix that summarizes the
#' results by showing the estimate power according to different sampling efforts.
#'
#' An object of class "ecocbo_beta" is a list containing the following components:
#' \itemize{
#'   \item \code{$Power} a data frame containing the estimation of power and beta for
#' several combination of sampling efforts (\code{m} sites and \code{n} samples).
#'   \item \code{$Results} a data frame containing the estimates of pseudoF for \code{simH0}
#' and \code{simHa}.
#'   \item \code{$alpha} level of significance for Type I error.
#' }
#'
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @references Underwood, A. J. (1997). Experiments in ecology: their logical
#' design and interpretation using analysis of variance. Cambridge university
#' press.
#' @references Underwood, A. J., & Chapman, M. G. (2003). Power, precaution,
#' Type II error and sampling design in assessment of environmental impacts.
#' Journal of Experimental Marine Biology and Ecology, 296(1), 49-70.
#' @references Anderson, M. J. (2014). Permutational multivariate analysis of
#' variance (PERMANOVA). Wiley statsref: statistics reference online, 1-15.
#' @references  Guerra‐Castro, E. J., Cajas, J. C., Simões, N., Cruz‐Motta, J.
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
#'
#' @examples
#' sim_beta(data = simResults, alpha = 0.05)
#'

sim_beta <- function(data, alpha = 0.05){

  if(!inherits(data, "ecocbo_data"))
    stop("data is not the right class(\"ecocbo_data\")")
  if(alpha >= 1)
    stop("alpha must be smaller than 1")

  # Calculate power and beta ----
  resultsHa <- as.data.frame(data$Results)
  fCrit <- stats::aggregate(resultsHa[,5],
                     by = list(resultsHa[,3], resultsHa[,4]),
                     stats::quantile, probs = (1 - alpha), type = 8, na.rm = T)
  totDim <- stats::aggregate(resultsHa[,5],
                      by = list(resultsHa[,3], resultsHa[,4]),
                      length)
  powr <- matrix(nrow = dim(fCrit)[1], ncol = 6)
  colnames(powr) <- c("m", "n", "Power", "Beta",
                      "fCrit", "index")
  powr[,1] <- totDim[,1]
  powr[,2] <- totDim[,2]
  powr[,5] <- fCrit[,3]
  powr[,6] <- seq(1, nrow(powr))

  for(i in powr[,6]){
    partialRes <- resultsHa[resultsHa$m == powr[i,1] & resultsHa$n == powr[i,2],]
    powr[i,3] <- dim(partialRes[partialRes[,6] >= powr[i,5],])[1] / totDim[i,3]
  }

  powr[,4] <- 1 - powr[,3]
  rowidx <- order(powr[,1], powr[,2])
  powr <- as.data.frame(powr[rowidx, c(1:5)])

  BetaResult <- list(Power = powr, Results = resultsHa, alpha = alpha, model = data$model)
  class(BetaResult) <- "ecocbo_beta"

  return(BetaResult)
}

#-------------------------------------------
## S3Methods print()
#-------------------------------------------

#' S3Methods for Printing
#'
#' @name prints
#'
#' @aliases
#' print.ecocbo_beta
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
# Print ecocbo_beta
#' @export
print.ecocbo_beta <- function(x, ...){
  x$Power[,3] <- round(x$Power[,3], 2)
  x1 <- stats::reshape(x$Power[,c(1:3)],
                direction = "wide",
                idvar = "m", timevar = "n",
                new.row.names = paste0("m = ",c(2:max(x$Power$m))))
  x1 <- x1[,-1]
  colnames(x1) <- paste0("n = ", c(2:max(x$Power$n)))
  cat("Power at different sampling efforts (m x n):\n")
  print(x1)
}
