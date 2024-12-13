#' Prepare data for evaluation in single-factor experiments
#'
#' \code{prep_data()} formats and arranges the initial data so that it can be
#' readily used by the other functions in the package. The function first gets
#' the species names and the number of samples for each species from the input
#' data frame. Then, it permutes the sampling efforts and calculates the pseudo-F
#' statistic and the mean squares for each permutation. Finally, it returns a
#' data frame with the permutations, pseudo-F statistic, and mean squares.
#'
#' @param data Data frame with species names (columns) and samples (rows)
#' information. The first column should indicate the site to which the sample
#' belongs, regardless of whether a single site has been sampled.
#' @param type Nature of the data to be processed. It may be presence / absence
#' ("P/A"), counts of individuals ("counts"), or coverage ("cover")
#' @param Sest.method Method for estimating species richness. The function
#' specpool is used for this. Available methods are the incidence-based Chao
#' "chao", first order jackknife "jack1", second order jackknife "jack2" and
#' Bootstrap "boot". By default, the "average" of the four estimates is used.
#' @param cases Number of data sets to be simulated.
#' @param N Total number of samples to be simulated in each site.
#' @param sites Total number of sites to be simulated in each data set.
#' @param n Maximum number of samples to consider.
#' @param m Maximum number of sites.
#' @param k Number of resamples the process will take. Defaults to 50.
#' @param transformation Mathematical function to reduce the weight of very
#' dominant species: 'square root', 'fourth root', 'Log (X+1)', 'P/A', 'none'
#' @param method The appropriate distance/dissimilarity metric (e.g. Gower,
#' Bray–Curtis, Jaccard, etc). The function [vegan::vegdist()] is called for
#' that purpose.
#' @param dummy Logical. It is recommended to use TRUE in cases where there are
#' observations that are empty.
#' @param useParallel Logical. Perform the analysis in parallel? Defaults to FALSE.
#'
#' @return \code{prep_data()} returns an object of class "ecocbo_data".
#'
#' An object of class "ecocbo_data" is a list containing: \code{$Results}, a data
#' frame that lists the estimates of pseudoF for \code{simH0} and \code{simHa}
#' that can be used to compute the statistical power for different sampling
#' efforts, as well as the square means necessary for calculating the variation
#' components.
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
#' [prep_data()]
#' [sim_beta()]
#' [plot_power()]
#' [sim_cbo()]
#' [scompvar()]
#'
#' @export
#' @import parallel
#' @import doParallel
#' @import foreach
#' @importFrom SSP assempar simdata
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom doSNOW registerDoSNOW
#'
#' @keywords internal

prep_data_single <- function(data, type = "counts", Sest.method = "average",
                      cases = 5, N = 100, sites = 3,
                      n, m = NULL, k = 50,
                      transformation = "none", method = "bray",
                      dummy = FALSE, useParallel = TRUE){

  # read data and store it in two objects, one for H0 and one for Ha ----
  datH0 <- data
  datH0[,1] <- as.factor("T___0")
  datHa <- data
  datHa[,1] <- as.factor(data[,1])

  # calculate simulation parameters, then simulate communities ----
  parH0 <- SSP::assempar(data = datH0,
                         type = type, Sest.method = Sest.method)
  parHa <- SSP::assempar(data = datHa,
                         type = type, Sest.method = Sest.method)

  simH0 <- SSP::simdata(parH0, cases = cases,
                        N = (N * sites), sites = 1)
  simHa <- SSP::simdata(parHa, cases = cases,
                        N = N, sites = sites)

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

  rm(simH0, simHa)

  H0Sim <- array(unlist(H0Sim), dim = c(xH0, yH0, casesHa))
  HaSim <- array(unlist(HaSim), dim = c(xH0, yH0, casesHa))

  labHa <- HaSim[,c((yH0-1):(yH0)),1]
  colnames(labHa) <- c("N", "M")

  # Stamp Ha labels to H0 (for M and N)
  H0Sim[,c((yH0-1):yH0),] <- labHa

  ## Helper matrix to store labels ----
  # NN <- casesHa * k * (m-1) * (n-1)
  # resultsHa <- matrix(nrow = NN, ncol = 8)
  # resultsHa[, 1] <- rep(seq(casesHa), times = 1, each = (k * (m-1) * (n-1)))
  # resultsHa[, 2] <- rep(1:k, times = (n-1) * (m-1) * casesHa)
  # resultsHa[, 3] <- rep(seq(2, m), times = (n-1), each = k)
  # resultsHa[, 4] <- rep(seq(2, n), times = 1, each = k * (m-1))
  # colnames(resultsHa) <- c("dat.sim", "k", "m", "n",
  #                          "pseudoFH0", "pseudoFHa",
  #                          "MSA", "MSR")

  # Changing the routine so that the results matrix will only include one case of m
  # that corresponds to the one we have in the original experiment.
  # We are taking out MSA as it is not relevant for the package.
  NN <- casesHa * k * (n-1)
  resultsHa <- matrix(nrow = NN, ncol = 7)
  colnames(resultsHa) <- c("dat.sim", "k", "m", "n",
                           "pseudoFH0", "pseudoFHa",
                           "MSR")

  resultsHa[, 1] <- rep(seq(casesHa), times = 1, each = (k * (n-1)))
  resultsHa[, 2] <- rep(1:k, times = (n-1) * casesHa)
  resultsHa[, 3] <- sites
  resultsHa[, 4] <- rep(seq(2, n), times = 1, each = k)


  # Loop to calculate pseudoF ----
  # Loop parameters
  Y <- cbind(1:(N * sites))
  YPU <- as.numeric(gl(sites, N))
  mm <- resultsHa[,3]
  nn <- resultsHa[,3] * resultsHa[,4]

  pb <- txtProgressBar(max = NN, style = 3)

  if(useParallel){
    cl <- parallel::makeCluster(parallel::detectCores() - 2)
    registerDoSNOW(cl)
    # doParallel::registerDoParallel(cl)

    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)

    parallel::clusterExport(cl, list("balanced_sampling", "permanova_oneway"))
    result1 <- foreach::foreach(i=1:NN, .combine = rbind,
                                .options.snow = opts) %dopar% {
                                  balanced_sampling(i, Y, mm, nn, YPU,
                                                    H0Sim, HaSim, resultsHa,
                                                    transformation, method)
                                }
    resultsHa[,5] <- result1[,1]
    resultsHa[,6] <- result1[,2]
    resultsHa[,7] <- result1[,4]
  } else {
    for (i in seq_len(NN)){
      result1 <- balanced_sampling(i, Y, mm, nn, YPU,
                                   H0Sim, HaSim, resultsHa,
                                   transformation, method)
      resultsHa[i,5] <- result1[,1]
      resultsHa[i,6] <- result1[,2]
      resultsHa[i,7] <- result1[,4]

      setTxtProgressBar(pb, i)

    }
  }

  stopCluster(cl)
  close(pb)

  resultsHa <- resultsHa[!is.na(resultsHa[,5]) & !is.na(resultsHa[,6]),]
  resultsHa <- resultsHa[,-3]

  SimResults <- list(Results = resultsHa, model = "single.factor")
  class(SimResults) <- "ecocbo_data"

  return(SimResults)

}

#' Prepare data for evaluation in nested symmetric double-factor experiments
#'
#' \code{prep_data()} formats and arranges the initial data so that it can be
#' readily used by the other functions in the package. The function first gets
#' the species names and the number of samples for each species from the input
#' data frame. Then, it permutes the sampling efforts and calculates the pseudo-F
#' statistic and the mean squares for each permutation. Finally, it returns a
#' data frame with the permutations, pseudo-F statistic, and mean squares.
#'
#' @param data Data frame with species names (columns) and samples (rows)
#' information. The first column should indicate the site to which the sample
#' belongs, regardless of whether a single site has been sampled.
#' @param type Nature of the data to be processed. It may be presence / absence
#' ("P/A"), counts of individuals ("counts"), or coverage ("cover")
#' @param Sest.method Method for estimating species richness. The function
#' specpool is used for this. Available methods are the incidence-based Chao
#' "chao", first order jackknife "jack1", second order jackknife "jack2" and
#' Bootstrap "boot". By default, the "average" of the four estimates is used.
#' @param cases Number of data sets to be simulated.
#' @param N Total number of samples to be simulated in each site.
#' @param sites Total number of sites to be simulated in each data set.
#' @param n Maximum number of samples to consider.
#' @param m Maximum number of sites.
#' @param k Number of resamples the process will take. Defaults to 50.
#' @param transformation Mathematical function to reduce the weight of very
#' dominant species: 'square root', 'fourth root', 'Log (X+1)', 'P/A', 'none'
#' @param method The appropriate distance/dissimilarity metric (e.g. Gower,
#' Bray–Curtis, Jaccard, etc). The function [vegan::vegdist()] is called for
#' that purpose.
#' @param dummy Logical. It is recommended to use TRUE in cases where there are
#' observations that are empty.
#' @param useParallel Logical. Perform the analysis in parallel? Defaults to FALSE.
#'
#' @return \code{prep_data()} returns an object of class "ecocbo_data".
#'
#' An object of class "ecocbo_data" is a list containing: \code{$Results}, a data
#' frame that lists the estimates of pseudoF for \code{simH0} and \code{simHa}
#' that can be used to compute the statistical power for different sampling
#' efforts, as well as the square means necessary for calculating the variation
#' components.
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
#' [prep_data()]
#' [sim_beta()]
#' [plot_power()]
#' [sim_cbo()]
#' [scompvar()]
#'
#' @export
#' @import parallel
#' @import doParallel
#' @import foreach
#' @importFrom SSP assempar simdata
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom doSNOW registerDoSNOW
#'
#' @keywords internal
#'

prep_data_nestedsymmetric <- function(data, type = "counts",
                                      Sest.method = "average",
                                      cases = 5, N = 100, sites = 3,
                                      n, m, k = 50,
                                      transformation = "none", method = "bray",
                                      dummy = FALSE, useParallel = TRUE){
  # get values for size limits
  data[,1] <- as.factor(data[,1])
  factSect <- data[,1]
  nSect <- nlevels(factSect)         # number of sectors

  # adapt data set so it only has 1 label column
  datSim <- data
  datSim[,2] <- as.factor(paste0(data[,1], "___", data[,2]))
  datSim <- datSim[,-1]


  # generate data for H0 and Ha
  datH0 <- datSim
  datH0[,1] <- "T___0"

  datHa <- datSim

  model <- "nested.symmetric"

  # calculate simulation parameters
  parH0 <- assempar(datH0, type = type, Sest.method = Sest.method)
  parHa <- assempar(datHa, type = type, Sest.method = Sest.method)

  # calculate simulated communities
  simH0 <- simdata(parH0, cases = cases, N = (N * sites * nSect), sites = 1)
  simHa <- simdata(parHa, cases = cases, N = (N * nSect), sites = sites)

  # design and fill the results matrix ----
  NN <- cases * k * (m-1) * (n-1)
  resultsHa <- matrix(nrow = NN, ncol = 8)
  colnames(resultsHa) <- c("dat.sim", "k", "m", "n",
                           "pseudoFH0", "pseudoFHa",
                           "MSB(A)", "MSR")

  resultsHa[, 1] <- rep(seq(cases), times = 1, each = (k * (m-1) * (n-1)))
  resultsHa[, 2] <- rep(1:k, times = (n-1) * (m-1) * cases)
  resultsHa[, 3] <- rep(seq(2, m), times = (n-1), each = k)
  resultsHa[, 4] <- rep(seq(2, n), times = 1, each = k * (m-1))

  # design the arrays to store the lists ----
  # H0 comes from simdata as-is, it does not include a column for sectors given that
  # its data was simulated by merging sectors and sites together
  H0Sim <- array(unlist(simH0), dim = c((N * sites * nSect), length(simH0[[1]]), cases))
  colnames(H0Sim) <- colnames(simH0[[1]])
  H0Sim <- H0Sim[,1:(dim(H0Sim)[2]-2),]

  HaSim <- array(unlist(simHa), dim = c((N * sites * nSect), length(simHa[[1]]), cases))
  colnames(HaSim) <- colnames(simHa[[1]])
  HaSim <- HaSim[,1:(dim(HaSim)[2]-2),]

  # Factors matrix and data matrix for Ha
  factEnv <- as.data.frame(simHa[[1]][,(dim(simHa[[1]])[2]-1):dim(simHa[[1]])[2]])
  factEnv$sites <- as.factor(rep(seq_len(sites), each = N, times = nSect))
  factEnv$sector <- as.factor(rep(levels(factSect), each = N*sites))

  rm(simH0, simHa)

  # Set of parameters for using balancedtwostage ----
  # index marking the size of each resampled site
  Y <- cbind(1:(N * sites))
  # index for the sites repeated N times
  YPU <- as.numeric(gl(sites, N))
  # labels for sites (i.e. (m-1) sites, repeated k times)
  mm <- resultsHa[,3]
  # number of samples to consider (i.e. m * n)
  nn <- resultsHa[,3] * resultsHa[,4]

  # progress bar
  pb <- txtProgressBar(max = NN, style = 3)

  if(useParallel){
    cl <- parallel::makeCluster(parallel::detectCores() - 2)
    doSNOW::registerDoSNOW(cl)

    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)

    parallel::clusterExport(cl, list("balanced_sampling2", "permanova_twoway", "SS"))
    result1 <- foreach::foreach(i=1:NN, .combine = rbind,
                                .options.snow = opts) %dopar% {
                                  balanced_sampling2(i, Y, mm, nn, YPU,
                                                     H0Sim, HaSim, factEnv, resultsHa,
                                                     transformation, method, model,
                                                     nSect, sites, N)
                                }
    resultsHa[,5] <- result1[,1]
    resultsHa[,6] <- result1[,2]
    resultsHa[,7] <- result1[,4]
    resultsHa[,8] <- result1[,5]

  } else {
    for (i in seq_len(NN)){
      result1 <- balanced_sampling2(i, Y, mm, nn, YPU,
                                    H0Sim, HaSim, factEnv, resultsHa,
                                    transformation, method, model,
                                    nSect, sites, N)
      resultsHa[i,5] <- result1[,1]
      resultsHa[i,6] <- result1[,2]
      resultsHa[i,7] <- result1[,4]
      resultsHa[i,8] <- result1[,5]

      setTxtProgressBar(pb, i)
    }
  }

  stopCluster(cl)
  close(pb)

  resultsHa <- resultsHa[!is.na(resultsHa[,5]) & !is.na(resultsHa[,6]),]

  SimResults <- list(Results = resultsHa, model = model)
  class(SimResults) <- "ecocbo_data"

  return(SimResults)

}
