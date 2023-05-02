#' Simulated cost-benefit optimization
#'
#'\code{sim_cbo} can be used to apply a cost-benefit optimization model that
#' depends either on a desired level of precision or on a budgeted total cost,
#' as proposed by Underwood (1997).
#'
#' @param comp.var Data frame as obtained from \code{\link{scompvar}}.
#' @param multSE Optional. Required multivariate standard error for the
#' sampling experiment.
#' @param ct Optional. Total cost for the sampling experiment.
#' @param ck Cost per replicate.
#' @param cj Cost per unit.
#'
#' @return A data frame containing the optimized values for \code{b} number of
#' sites and \code{n} number of samples to consider.
#'
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @references Underwood, A. J., Underwood, A. J., & Wnderwood, A. J. (1997).
#' Experiments in ecology: their logical design and interpretation using
#' analysis of variance. Cambridge university press.
#' @references Underwood, A. J., & Chapman, M. G. (2003). Power, precaution,
#' Type II error and sampling design in assessment of environmental impacts.
#' Journal of Experimental Marine Biology and Ecology, 296(1), 49-70.
#'
#' @seealso
#' \code{\link{sim_beta}}
#' \code{\link{plot_power}}
#' \code{\link{scompvar}}
#'
#' @aliases simcbo
#'
#' @export
#'
#' @examples
#' compVar <- scompvar(data = epiBetaR)
#'
#' sim_cbo(comp.var = compVar, multSE = NULL, ct = 20000, ck = 100, cj = 2500)
#' sim_cbo(comp.var = compVar, multSE = 0.15, ct = NULL, ck = 100, cj = 2500)

sim_cbo <- function(comp.var, multSE = NULL, ct = NULL, ck, cj){
# Optimal cost-benefit model

# Helper functions ----
# Function to determine optimal b by setting costs.
cost_n <- function(n, ct, ck, cj){
  b <- data.frame(nOpt = n, bOpt = NA)

  # Using equation 9.19 (Underwood, 1997)
  b[,2] <- floor(ct / (n * ck + cj))
  b[,1] <- floor(b[,1])

  return(b)
}

# Function to determine optimal b by setting desired variability.
cost_v <- function(comp.var, multSE, n){
  b <- data.frame(nOpt = n, bOpt = NA)

  # Using equation 9.18 (Underwood, 1997)
  b[,2] <- floor((comp.var[,2] + b$nOpt * comp.var[,1]) /
                   (multSE * multSE * b$nOpt))
  b[,1] <- floor(b[,1])

  return(b)
}

# Main function ----
  ## Validating data ----
  if(is.null(multSE) & is.null(ct)){
    stop("It is necessary to provide either multSE or ct")
  }
  if(dim(comp.var)[1] != 1 | dim(comp.var)[2] != 2){
    stop("Variation components must be in a 1x2 matrix")
  }

  ## Calculate optimal n ----
  nOpt <- sqrt((cj * comp.var[,2]) / (ck * comp.var[,1]))

  ## Calculate optimal b ----
  if(is.null(multSE)) {
    b <- cost_n(nOpt, ct, ck, cj)
  } else {
    b <- cost_v(comp.var, multSE, nOpt)
  }
  return(b)
}

