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
#' @importFrom stats reshape
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
#' sim_beta(simH0, simHa, n = 10, m = 3, k = 20, alpha = 0.05)

sim_beta <- function(simH0, simHa, n, m, k= 50, alpha = 0.05){
  # Calculate power and pseudoF in multiple iterations ----

  N <- max(simHa[[1]][,'N'])
  sites <- max(as.numeric(simHa[[1]][,'sites']))
  if (n > N){stop("'n' must be equal or less than 'N'")}
  if(ceiling(n) != floor(n)){stop("n must be integer")}
  if(n <= 1){stop("n must be larger than 1")}

  if (m > sites){stop("'m' must be equal or less than 'sites'")}
  if(ceiling(m) != floor(m)){stop("m must be integer")}
  if(m <= 1){stop("m must be larger than 1")}

  if(alpha >= 1){stop("alpha must be smaller than 1")}

  if(!is.list(simH0) | !is.list(simHa)){stop("simulation data must be lists")}


  xH0 <- dim(simH0[[1]])[1]
  yH0 <- dim(simH0[[1]])[2]
  zH0 <- length(simH0)
  H0Sim <- array(unlist(simH0), dim = c(xH0, yH0, zH0))
  HaSim <- array(unlist(simHa), dim = c(xH0, yH0, zH0))

  # Simulation parameters ----
  casesHa <- dim(HaSim)[3]

  labHa <- HaSim[,c((yH0-1):yH0),1]
  colnames(labHa) <- c("N", "sites")
  labHa <- cbind(labHa, index = 1:xH0)

  popH0 <- max(labHa[,1])

  # Labels from Ha to H0 (for sites and N)
  H0Sim[,c((yH0-1):yH0),] <- labHa[,c(1:2)]

  ## Helper matrix to store labels ----
  resultsHa <- matrix(nrow = casesHa * k * (m-1) * (n-1), ncol = 8)
  resultsHa[, 1] <- rep(seq(casesHa), times = 1, each = (k * (m-1) * (n-1)))
  resultsHa[, 2] <- rep(1:k, times = (n-1) * (m-1) * casesHa)
  resultsHa[, 3] <- rep(seq(2, m), times = (n-1), each = k)
  resultsHa[, 4] <- rep(seq(2, n), times = 1, each = k * (m-1))
  colnames(resultsHa) <- c("dat.sim", "k", "m", "n",
                           "pseudoFH0", "pseudoFHa",
                           "AMSHa", "RMSHa")

  # Loop to calculate pseudoF ----
  # Loop parameters
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

    # Replace incorrect probabilities (p < 0 and p > 1) produced by balancetwostage
    sel[sel[,1]<= -1, 1] <- 0
    sel[sel[,1]>= 2, 1] <- 1

    # Getting data
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

  # Calculate power and beta ----
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
#' @description Prints for \code{ecocbo::sim_beta} objects
#'
#' @param x Object from \code{ecocbo::sim_beta} package
#'
#' @param ... Additional arguments
#'
#' @return Prints \code{ecocbo::sim_beta} object
#'
# Print ecocbo_beta
#' @export
print.ecocbo_beta <- function(x, ...){
  x$Power[,3] <- round(x$Power[,3], 3)
  x1 <- reshape(x$Power[,c(1:3)],
                direction = "wide",
                idvar = "m", timevar = "n",
                new.row.names = paste0("m = ",c(2:max(x$Power$m))))
  x1 <- x1[,-1]
  colnames(x1) <- paste0("n = ", c(2:max(x$Power$n)))
  cat("Power at different sampling efforts (m x n):\n")
  print(x1)
}
