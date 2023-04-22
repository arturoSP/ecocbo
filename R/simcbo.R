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
#' library("SSP")
#' # Load data and adjust it.
#' epiH0 <- epiDat
#' epiH0[,"site"] <- as.factor("T0")
#' epiHa <- epiDat
#' epiHa[,"site"] <- as.factor(epiHa[,"site"])
#'
#' # Calculate simulation parameters.
#' parH0 <- SSP::assempar(data = epiH0, type = "counts", Sest.method = "average")
#' parHa <- SSP::assempar(data = epiHa, type = "counts", Sest.method = "average")
#'
#' # Simulation.
#' simH0 <- SSP::simdata(parH0, cases = 3, N = 1000, sites = 1)
#' simHa <- SSP::simdata(parHa, cases = 3, N = 100, sites = 10)
#'
#' betaResults <- beta(simH0, simHa, n = 10, m = 3, k = 50, alpha = 0.05)
#'
#' compVar <- scompvar(data = betaResults)
#'
#' simcbo(comp.var = compVar, multSE = NULL, ct = 600, ck = 11, cj = 108)
#' simcbo(comp.var = compVar, multSE = 0.15, ct = NULL, ck = 11, cj = 108)

simcbo <- function(comp.var, multSE, ct, ck, cj){
# función para calcular el model de optimización de costo-beneficio

# funciones de apoyo ----
# función para determinar el número de sitios a muestrear, conociendo
# costos de muestreo y número de muestras a tomar.
cost_n <- function(n, ct, ck, cj){
  b <- data.frame(nOpt = n, bOpt = NA)

  # aplica la ecuación 9.19 de Underwood, 1997
  b[,2] <- round(ct / (n * ck + cj))
  b[,1] <- round(b[,1])

  return(b)
}

# función para determinar el número de sitios a muestrear, conociendo
# la variabilidad desde sampsd()
cost_v <- function(comp.var, multSE, n){
  b <- data.frame(nOpt = n, bOpt = NA)

  # Aplica la ecuación 9.18 de Underwood, 1997
  b[,2] <- round((comp.var$compVarR + b$nOpt * comp.var$compVarA) / (multSE * multSE * b$nOpt))
  b[,1] <- round(b[,1])

  return(b)
}

# función principal ----
#simcbo <- function(comp.var, multSE, ct, ck, cj){
  # Cálculo de n óptimo ----
  nOpt <- sqrt((cj * comp.var$compVarR) / (ck * comp.var$compVarA))

  # Cálculo de b óptimo ----
  if(is.null(multSE)) {
    b <- cost_n(nOpt, ct, ck, cj)
  } else {
    b <- cost_v(comp.var, multSE, nOpt)
  }
  return(b)
}

