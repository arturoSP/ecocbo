#' Calculate beta and power out of simulated samples
#'
#' @param simH0 Simulated community from SSP::simdata in which H0 is true.
#' @param simHa Simulated community from SSP::simdata in which H0 is false.
#' @param n Maximum number of samples to consider.
#' @param m Maximum number of sites.
#' @param k Number of resamples the process will take. Defaults to 50.
#' @param alpha Level of significance. Defaults to 5%.
#'
#' @return A list with two data frames: a data frame containing the values of beta for different sampling efforts, and a data frame containing the results of the sampling.
#' @export
#'
#' @examples
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
#' sim_beta(simH0, simHa, n = 10, m = 3, k = 50, alpha = 0.05)

sim_beta <- function(simH0, simHa, n, m, k= 50, alpha = 0.05){
  # Cálculo de potencia y simulación de valores pseudoF en múltiples iteraciones ----

  N <- max(simHa[[1]][,'N'])
  sites <- max(as.numeric(simHa[[1]][,'sites']))
  if (n > N){
    stop("'n' must be equal or less than 'N'")}
  if (m > sites){
    stop("'m' must be equal or less than 'sites'")}

  xH0 <- dim(simH0[[1]])[1]
  yH0 <- dim(simH0[[1]])[2]
  zH0 <- length(simH0)
  H0Sim <- array(unlist(simH0), dim = c(xH0, yH0, zH0))
  HaSim <- array(unlist(simHa), dim = c(xH0, yH0, zH0))

  # Parámetros de simulación ----
  casesHa <- dim(HaSim)[zH0]

  labHa <- HaSim[,c((yH0-1):yH0),1]
  colnames(labHa) <- c("N", "sites")
  labHa <- cbind(labHa, index = 1:xH0)
           #as.matrix(dplyr::bind_cols(labHa, row_number(labHa[,2])))

  popH0 <- max(labHa[,1])

  # Se asignan etiquetas de Ha a H0 (sitio y N)
  H0Sim[,c((yH0-1):yH0),] <- labHa[,c(1:2)]

  ## se crea matriz que almacena las etiquetas de sitio (factores) ----
  resultsHa <- matrix(nrow = casesHa * k * (m-1) * (n-1), ncol = 8)
  resultsHa[, 1] <- rep(seq(casesHa), times = 1, each = (k * (m-1) * (n-1)))
  resultsHa[, 2] <- rep(1:k, times = (n-1) * (m-1) * casesHa)
  resultsHa[, 3] <- rep(seq(2, m), times = (n-1), each = k)
  resultsHa[, 4] <- rep(seq(2, n), times = 1, each = k * (m-1))
  colnames(resultsHa) <- c("dat.sim", "k", "m", "n",
                           "pseudoFH0", "pseudoFHa",
                           "AMSHa", "RMSHa")

  # loop para calcular pseudoF ----
  # objetos para el loop
  Y <- cbind(1:(N * sites))
  YPU <- as.numeric(as.vector(gl(sites, N)))
  NN <- nrow(resultsHa)
  mm <- resultsHa[,3]
  nn <- rep(NA, NN)
  for (i in seq_len(NN)){
    nn[i] <- resultsHa[i, 3] * resultsHa[i, 4]
  }

  for (i in seq_len(NN)){
    sel <- sampling::balancedtwostage(Y, selection = 1, m = mm[i], n = nn[i], PU = YPU, FALSE)

    #replace incorrect probabilities (p < 0 and p > 1) produced by balancetwostage
    sel[sel[,1]<= -1, 1] <- 0
    sel[sel[,1]>= 2, 1] <- 1

    #getting data
    dat.H0 <- H0Sim
    dat.Ha <- HaSim
    rownames(sel) <- Y
    ones <- sel[, 1]
    y0 <- dat.H0[,,resultsHa[i,1]][ones == 1,]
    ya <- dat.Ha[,,resultsHa[i,1]][ones == 1,]
    perH0 <- as.data.frame(y0[ , 1:(yH0-2)])
    perHa <- as.data.frame(ya[ , 1:(yH0-2)])
    perH0Env <- y0[,yH0]

    result1 <- ecocbo::permanova_reduced(perH0, perH0Env)
    result2 <- ecocbo::permanova_reduced(perHa, perH0Env)
    resultsHa[i,5] <- result1$Fobs
    resultsHa[i,6] <- result2$Fobs
    resultsHa[i,7] <- result2$AMS
    resultsHa[i,8] <- result2$RMS
  }

  # Cálculo de potencia y beta ----
  resultsHa <- as.data.frame(resultsHa)
  fCrit <- stats::aggregate(resultsHa[,5],
                     by = list(resultsHa$m, resultsHa$n),
                     stats::quantile, probs = (1 - alpha), type = 8, na.rm = T)
  totDim <- stats::aggregate(resultsHa[,5],
                      by = list(resultsHa$m, resultsHa$n),
                      length)
  AMSHa <- stats::aggregate(resultsHa[,8],
                     by = list(resultsHa$m, resultsHa$n),
                     mean, na.rm = T)
  RMSHa <- stats::aggregate(resultsHa[,7],
                     by = list(resultsHa$m, resultsHa$n),
                     mean, na.rm = T)
  powr <- matrix(nrow = dim(fCrit)[1], ncol = 8)
  colnames(powr) <- c("m", "n", "Power", "Beta",
                      "fCrit", "AMSHa", "RMSHa", "index")
  powr[,1] <- totDim[,1]
  powr[,2] <- totDim[,2]
  powr[,5] <- fCrit[,3]
  powr[,6] <- AMSHa[,3]
  powr[,7] <- RMSHa[,3]
  powr[,8] <- seq(1, nrow(powr))

  for(i in powr[,8]){
    partialRes <- resultsHa[resultsHa$m == powr[i,1] & resultsHa$n == powr[i,2],]
    powr[i,3] <- dim(partialRes[partialRes[,6] >= powr[i,5],])[1] / totDim[i,3]
  }

  powr[,4] <- 1 - powr[,3]
  rowidx <- order(powr[,1], powr[,2])
  powr <- as.data.frame(powr[rowidx, c(1:7)])

  BetaResult <- list(Power = powr, Results = resultsHa)
  return(BetaResult)
}
