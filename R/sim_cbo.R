#' Simulated cost-benefit optimization
#'
#'\code{sim_cbo()} can be used to apply a cost-benefit optimization model that
#' depends either on a desired level of precision or on a budgeted total cost,
#' as proposed by Underwood (1997).
#'
#' @param comp.var Data frame as obtained from [scompvar()].
#' @param multSE Optional. Required multivariate standard error for the
#' sampling experiment.
#' @param budget Optional. Total cost for the sampling experiment.
#' @param ca Cost per treatment.
#' @param cm Cost per replicate.
#' @param cn Cost per unit.
#'
#' @return A data frame containing the optimized values for \code{m} number of
#' sites and \code{n} number of samples to consider.
#'
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @references Underwood, A. J. (1997). Experiments in ecology: their logical
#' design and interpretation using analysis of variance. Cambridge university
#' press.
#' @references Underwood, A. J., & Chapman, M. G. (2003). Power, precaution,
#' Type II error and sampling design in assessment of environmental impacts.
#' Journal of Experimental Marine Biology and Ecology, 296(1), 49-70.
#'
#' @seealso
#' [sim_beta()]
#' [plot_power()]
#' [scompvar()]
#'
#' @aliases simcbo
#'
#' @export
#'
#' @examples
#' compVar <- scompvar(data = simResults)
#'
#' sim_cbo(comp.var = compVar, multSE = NULL, ct = 20000, ck = 100, cj = 2500)
#' sim_cbo(comp.var = compVar, multSE = 0.15, ct = NULL, ck = 100, cj = 2500)

sim_cbo <- function(comp.var, multSE = NULL, budget = NULL, a = NULL,
                    ca = NULL, cm = NULL, cn){

# Optimal cost-benefit model

# Helper functions ----
# Function to determine optimal m by setting costs.
cost_n <- function(n, budget, cm, cn){
  m <- data.frame(nOpt = n, mOpt = NA)

  # Using equation 9.19 (Underwood, 1997)
  m[,2] <- floor(budget / (n * cn + cm))

  return(m)
  }

# Function to determine optimal m by setting desired variability.
cost_v <- function(n, comp.var, multSE){
  m <- data.frame(nOpt = n, mOpt = NA)

  # Using equation 9.18 (Underwood, 1997)
  m[,2] <- floor((comp.var[2,2] + n * comp.var[1,2]) /
                   (multSE * multSE * n))

  return(m)
  }

# Main function ----
  ## Validating data ----
  if(is.null(multSE) & is.null(budget)){
    stop("It is necessary to provide either multSE or ct")
  }

  if(dim(comp.var)[1] == 1 & is.null(ca) & is.null(multSE) |
     dim(comp.var)[1] == 1 & is.null(a) & is.null(multSE)) {
    stop("For single factor experiments, it is necessary to provide the
         number of treatments and its cost")
  }

  if(dim(comp.var)[1] == 2 & is.null(cm)){
    stop("Cost per unit is required for a multivariate model")
  }

  ## Calculate optimal n ----
  if(dim(comp.var)[1] == 1){
    if(is.null(multSE)){
      # when working for total cost, the optimal n will be what can be
      # done with the available economic resources
      nOpt <- floor((budget - (ca * a)) / cn)
    } else {
      # when working with SE, the optimal n comes from solving Var=MSR/n for n
      nOpt <- floor(comp.var[1,2] / (multSE * multSE))
    }
    m <- data.frame(nOpt)
  } else if(dim(comp.var)[1] == 2){
    nOpt <- floor(sqrt((cm * comp.var[2,2]) / (cn * comp.var[1,2])))

    ## Calculate optimal m ----
    if(is.null(multSE)) {
      m <- cost_n(nOpt, budget, cm, cn)
    } else {
      m <- cost_v(nOpt, comp.var, multSE)
    }
  }
  return(m)
}

