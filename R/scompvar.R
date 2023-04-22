#' Simulated components of variation
#'
#' @param data List that results from beta().
#' @param m Site label to be used as basis for the computation.
#' @param n Number of samples to be considered.
#'
#' @return A data frame containing the values for the variation component among sites and in the residuals.
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
#' scompvar(data = betaResults)
#' scompvar(data = betaResults, m = 2, n = 5)

scompvar <- function(data, m = NULL, n = NULL){
  # Función para determinar los componentes de variación previo al modelo de
  # optimización de costo-beneficio

  # lee la tabla de resultados desde la lista de Beta para poder usar los MeanSquares de Ha
  resultsBeta <- data$Results
  if(is.null(m)){m <- max(resultsBeta$m)}
  if(is.null(n)){n <- max(resultsBeta$n)}
  resultsBeta <- resultsBeta[resultsBeta$m == m & resultsBeta$n == n,]

  # Crea un data frame con tamaño: número juegos de datos simulados x 3
  compVar <- data.frame(compVarA = NA, compVarR = NA)

  # Calcula componentes de variación promedio (de la Tabla 9.3)
  # σe = RMS
  compVar[,2] <- mean(resultsBeta[,8], na.rm = T)
  # σB(A) = (AMS - σe) / n
  MSA <- mean(resultsBeta[,7], na.rm = T)
  compVar[,1] <- (MSA - compVar[,2]) / n

  return(compVar)
}


