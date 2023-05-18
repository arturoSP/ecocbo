#' Calculate beta and power out of simulated samples
#'
#' \code{sim_beta} can be used to assess the power of a study by comparing the
#' variation when one can assume wether an ecological community does not have
#' composition differences (H0 true) or it does (H0 false). For example, if the
#' beta error is 0.25, then there is a 25% chance of failing to detect a
#' difference even if the difference is real. The power of the study is
#' \eqn{1 - \beta}, so in this example, the power of the study is 0.75.
#'
#' @param simH0 Simulated community from \code{SSP::simdata} in which H0 is
#' true.
#' @param simHa Simulated community from \code{SSP::simdata} in which H0 is
#' false.
#' @param n Maximum number of samples to consider.
#' @param m Maximum number of sites.
#' @param k Number of resamples the process will take. Defaults to 50.
#' @param alpha Level of significance for Type I error. Defaults to 0.05.
#' @param transformation Mathematical function to reduce the weight of very
#' dominant species: 'square root', 'fourth root', 'Log (X+1)', 'P/A', 'none'
#' @param method The appropriate distance/dissimilarity metric (e.g. Gower,
#' Bray–Curtis, Jaccard, etc). The function \code{\link[vegan]{vegdist}} is
#' called for that purpose.
#' @param dummy Logical. It is recommended to use TRUE in cases where there are
#' observations that are empty.
#' @param useParallel Logical. Perform the analysis in parallel? Defaults to TRUE.
#' @param nCores Number of cores to use if the analysis is performed in parallel.
#'
#' @return \code{sim_data} returns an object of class "ecocbo_beta".
#'
#' The function \code{print} is used to present a matrix that summarizes the
#' results by showing the estimate power according to different sampling efforts.
#'
#' An object of class "ecocbo_beta" is a list containing the following components:
#' \itemize{
#'   \item \code{Power} a data frame containing the estimation of power and beta for
#' several combination of sampling efforts (\code{m} sites and \code{n} samples).
#'   \item \code{Results} a data frame containing the estimates of pseudoF for \code{simH0}
#' and \code{simHa}.
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
#' \code{\link[SSP]{assempar}}
#' \code{\link[SSP]{simdata}}
#' \code{\link{plot_power}}
#' \code{\link{scompvar}}
#' \code{\link{sim_cbo}}
#'
#' @aliases simbeta
#'
#' @export
#' @importFrom stats reshape
#' @importFrom stats aggregate
#' @importFrom stats quantile
#' @importFrom sampling balancedtwostage
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#'
#' @examples
#' sim_beta(simH0Dat, simHaDat, n = 10, m = 3, k = 20, alpha = 0.05,
#' transformation = "square root", method = "bray", dummy = FALSE,
#' nCores = 2)
#'

sim_beta <- function(simH0, simHa, n, m, k= 50, alpha = 0.05,
                     transformation = "none", method = "bray", dummy = FALSE,
                     useParallel = TRUE, nCores = NULL){
  # Helper functions ----

  permanova_oneway <- function(x, factEnv, type = "P", method = "bray", transformation = "none"){

    # Helper functions ----
    ## Sum of Squares using Huygen theorem ----
    SS <- function (d) {
      ss <- numeric(2)
      ss[1] <- dim(as.matrix(d))[1]
      ss[2] <- sum(d^2)/ss[1]
      return(ss)
    }

    ## PERMANOVA function ----
    pseudoF_P <- function(x, factEnv, method = "bray", transformation = "none"){
      if (transformation == "square root") {
        x.t <- sqrt(x)
        rm(x)
        d <- vegan::vegdist(x.t, method = method)
      }
      if (transformation == "fourth root") {
        x.t <- sqrt(sqrt(x))
        rm(x)
        d <- vegan::vegdist(x.t, method = method)
      }
      if (transformation == "Log (X+1)") {
        x.t <- log(x + 1)
        rm(x)
        d <- vegan::vegdist(x.t, method = method)
      }
      if (transformation == "P/A") {
        x.t <- 1 * (x > 0)
        rm(x)
        d <- vegan::vegdist(x.t, method = method, binary = TRUE)
      }
      if (transformation == "none") {
        x.t <- x
        rm(x)
        d <- vegan::vegdist(x.t, method = method)
      }

      TSS <- SS(d)[2]

      # Size for labels
      nlev <- nlevels(as.factor(factEnv))

      # Calculate the SS for residuals
      lista <- split(as.data.frame(x.t), factEnv)
      dimxt <- dim(x.t)
      vdist <- lapply(lista, vegan::vegdist, method = method)
      SSi <- lapply(vdist, SS)
      SSi <- array(unlist(SSi), dim = c(1,2,nlev))

      # Calculate denominators
      denA <- nlev - 1
      denR <- dimxt[1] - nlev

      # Results
      RSS <- sum(SSi[,2,])
      ASS <- abs(TSS - RSS)
      AMS <- (ASS/denA)
      RMS <- (RSS/denR)
      Fobs <- AMS/RMS
      Fobs <- data.frame(ASS, RSS, TSS, denA, denR, AMS, RMS, Fobs)
      return(Fobs)
    }

    # Main function permanova ----
    Results <- pseudoF_P(x, factEnv, method, transformation)[1,c(6,7,8)]
    Fobs <- Results
    return(Fobs)
  }

  balanced_sampling <- function(i, Y, mm, nn, YPU, H0Sim, HaSim, resultsHa, transformation, method){
    # main function sampling ----
    # Get the samples index
    sel <- sampling::balancedtwostage(Y, selection = 1, m = mm[i],
                                      n = nn[i], PU = YPU, FALSE)
    ones <- which(sel[,1] %in% 1)
    y0 <- H0Sim[ones,,resultsHa[i,1]]
    ya <- HaSim[ones,,resultsHa[i,1]]
    yHa <- dim(y0)[2] - 2

    # Apply PERMANOVA to get F and mean squares
    result1 <- permanova_oneway(y0[, 1:yHa], y0[,yHa+2],
                                transformation = transformation, method = method)
    result2 <- permanova_oneway(ya[, 1:yHa], y0[,(yHa+2)],
                                transformation = transformation, method = method)
    result0 <- matrix(nrow = 1, ncol = 4)
    colnames(result0) <- c("Fobs", "Fobs", "AMS", "RMS")

    # Gather the results and return
    result0[,1] <- result1[,3]
    result0[,2] <- result2[,3]
    result0[,3] <- result2[,1]
    result0[,4] <- result2[,2]
    return(result0)
  }

  # Main function sim_beta ----
  ## Check the inputs ----

  if(useParallel){
    if(is.null(nCores)){
      nCores <- parallel::detectCores()-2
      doParallel::registerDoParallel(cores = nCores)
    }
    if(!is.null(nCores)){
      doParallel::registerDoParallel(cores = nCores)
    }
  }

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

  if(length(simH0) != length(simHa) |
     dim(simH0[[1]])[1] != dim(simHa[[1]])[1] |
     dim(simH0[[1]])[2] != dim(simHa[[1]])[2]){stop("dimensions for simH0 and simHa do not match")}

  if(nCores < 2){useParallel = FALSE}

  # Simulation parameters ---
  xH0 <- dim(simHa[[1]])[1]
  yH0 <- dim(simHa[[1]])[2]
  casesHa <- length(simHa)

  H0Sim <- simH0
  HaSim <- simHa

  if(dummy == TRUE){
    yH0 <- yH0 + 1
    for(i in seq_len(casesHa)){
      H0Sim[[i]] <- cbind(simH0[[i]], dummy = 1)
      H0Sim[[i]] <- H0Sim[[i]][,c(1:(yH0-3), (yH0), (yH0-2):(yH0-1))]
      HaSim[[i]] <- cbind(simHa[[i]], dummy = 1)
      HaSim[[i]] <- HaSim[[i]][,c(1:(yH0-3), (yH0), (yH0-2):(yH0-1))]
    }
  }

  H0Sim <- array(unlist(H0Sim), dim = c(xH0, yH0, casesHa))
  HaSim <- array(unlist(HaSim), dim = c(xH0, yH0, casesHa))

  labHa <- HaSim[,c((yH0-1):(yH0)),1]
  colnames(labHa) <- c("N", "sites")

  # Stamp Ha labels to H0 (for sites and N)
  H0Sim[,c((yH0-1):yH0),] <- labHa

  ## Helper matrix to store labels ----
  NN <- casesHa * k * (m-1) * (n-1)
  resultsHa <- matrix(nrow = NN, ncol = 8)
  resultsHa[, 1] <- rep(seq(casesHa), times = 1, each = (k * (m-1) * (n-1)))
  resultsHa[, 2] <- rep(1:k, times = (n-1) * (m-1) * casesHa)
  resultsHa[, 3] <- rep(seq(2, m), times = (n-1), each = k)
  resultsHa[, 4] <- rep(seq(2, n), times = 1, each = k * (m-1))
  colnames(resultsHa) <- c("dat.sim", "k", "m", "n",
                           "pseudoFH0", "pseudoFHa",
                           "AMSHa", "RMSHa")

  ## Loop to calculate pseudoF ----
  # Loop parameters
  Y <- cbind(1:(N * sites))
  YPU <- as.numeric(gl(sites, N))
  mm <- resultsHa[,3]
  nn <- resultsHa[,3] * resultsHa[,4]

  if(useParallel){
    result1 <- foreach::foreach(i=1:NN, .combine = rbind) %dopar% {
      balanced_sampling(i, Y, mm, nn, YPU,
                        H0Sim, HaSim, resultsHa,
                        transformation, method)
    }
  } else {
    result1 <- foreach::foreach(i=1:NN, .combine = rbind) %do% {
      balanced_sampling(i, Y, mm, nn, YPU,
                        H0Sim, HaSim, resultsHa,
                        transformation, method)
    }
  }

  resultsHa[,5] <- result1[,1]
  resultsHa[,6] <- result1[,2]
  resultsHa[,7] <- result1[,3]
  resultsHa[,8] <- result1[,4]

  resultsHa <- resultsHa[!is.na(resultsHa[,5] | !is.na(resultsHa[,6])),]

  ## Calculate power and beta ----
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
  powr <- as.data.frame(powr[rowidx, c(1:5)])

  BetaResult <- list(Power = powr, Results = resultsHa)
  class(BetaResult) <- "ecocbo_beta"
  return(BetaResult)
}

#
## S3Methods print() ----
#

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
#' @description Prints for \code{ecocbo::sim_beta} objects.
#'
#' @param x Object from \code{ecocbo::sim_beta} function.
#'
#' @param ... Additional arguments
#'
#' @return Prints \code{ecocbo::sim_beta} function.
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
