#' Simulated cost-benefit optimization
#'
#' @param comp.var Data frame from scompvar().
#' @param multSE Optional. Required multivariate standard error for the sampling experiment.
#' @param ct Optional. Total cost for the sampling experiment.
#' @param ck Cost per replicate.
#' @param cj Cost per unit.
#'
#' @return A data frame containing the optimized values for b (number of sites) and n (number of samples) to consider.
#' @export
#'
#' @examples
#' compVar <- scompvar(data = epiBetaR)
#'
#' simcbo(comp.var = compVar, multSE = NULL, ct = 20000, ck = 100, cj = 2500)
#' simcbo(comp.var = compVar, multSE = 0.15, ct = NULL, ck = 100, cj = 2500)

simcbo <- function(comp.var, multSE = NULL, ct = NULL, ck, cj){
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
  b[,2] <- floor((comp.var[,2] + b$nOpt * comp.var[,1]) / (multSE * multSE * b$nOpt))
  b[,1] <- floor(b[,1])

  return(b)
}

# Main function ----
  ## Validating data ----
  if(is.null(multSE) & is.null(ct)){stop("It is necessary to provide either multSE or ct")}
  if(dim(comp.var)[1] != 1 | dim(comp.var)[2] != 2){stop("Variation components must be in a 1x2 matrix")}

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

