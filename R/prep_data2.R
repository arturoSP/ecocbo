#' Prepare data for evaluation
#'
#' \code{prep_data2()} formats and arranges the initial data so that it can be
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
#' An object of class "ecocbo_data" is a list containing: \code{$Results} a data
#' frame containing the estimates of pseudoF for \code{simH0} and \code{simHa}
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
#' [sim_beta()]
#' [plot_power()]
#' [sim_cbo()]
#' [scompvar()]
#'
#' @aliases prepdata2
#'
#' @export
#' @import parallel
#' @import doParallel
#' @import foreach
#' @importFrom doSNOW registerDoSNOW
#' @importFrom SSP assempar
#' @importFrom SSP simdata
#' @importFrom stringr str_replace_all
#' @importFrom dplyr select
#' @importFrom plyr rbind.fill
#'
#'
#' @examples
#' \donttest{
#' simResults <- prep_data2(data = epiDat, type = "counts", Sest.method = "average",
#'                         cases = 5, N = 100, sites = 10,
#'                         n = 5, m = 5, k = 30,
#'                         transformation = "none", method = "bray",
#'                         dummy = FALSE, useParallel = FALSE, nested = FALSE)
#' }
#' simResults
#'

prep_data2 <- function(data, type = "counts", Sest.method = "average",
                       cases = 5, N = 100, sites = 10,
                       n, m, k = 50,
                       transformation = "none", method = "bray",
                       dummy = FALSE, useParallel = TRUE,
                       model = "nested.symmetric"){
  # Check the inputs ----
  if (n > N){stop("'n' must be equal or less than 'N' on simulated data")}
  if(ceiling(n) != floor(n)){stop("n must be integer")}
  if(n <= 1){stop("n must be larger than 1")}

  if (m > sites){stop("'m' must be equal or less than 'sites' on simulated data")}
  if(ceiling(m) != floor(m)){stop("m must be integer")}
  if(m <= 1){stop("m must be larger than 1")}

  # split data by sectors ----
  data[,1] <- as.factor(data[,1])
  data[,2] <- as.factor(data[,2])
  factSect <- data[,1]
  dataSpl <- split(data, factSect)
  dataSpl <- lapply(dataSpl, dplyr::select, -1)

  nSect <- nlevels(factSect)
  nR <- ncol(data)
  nC <- max(summary(factSect))

  # simulate H0 and Ha ----
  dat2H0 <- lapply(dataSpl, function(df){
    df[,1] <- "T0"
    df[,1] <- as.factor(df[,1])
    return(df)
  })

  dat2Ha <- lapply(dataSpl, function(df){
    df[,1] <- as.factor(df[,1])
    return(df)
  })

  # calculate simulation parameters, then simulate communities H0 and Ha ----
  parH0 <- lapply(dat2H0, assempar, type = type, Sest.method = Sest.method)
  parHa <- lapply(dat2Ha, assempar, type = type, Sest.method = Sest.method)

  for(i in levels(factSect)){
    parH0[[i]]$par$Species <- stringr::str_replace_all(parH0[[i]]$par$Species,
                                                       "unseen.species ",
                                                       paste0("unseen.species.",i,"."))
    parHa[[i]]$par$Species <- stringr::str_replace_all(parHa[[i]]$par$Species,
                                                       "unseen.species ",
                                                       paste0("unseen.species.",i,"."))
  }

  simH0 <- lapply(parH0, simdata, cases = cases, N = (N*sites), sites = 1)
  simHa <- lapply(parHa, simdata, cases = cases, N = N, sites = sites)

  # create blank matrixes to store registries ----
  emptyH0 <- data.frame(matrix(nrow = 0, ncol = ncol(data)-2))
  colnames(emptyH0) <- colnames(data[,-c(1,2)])
  emptyHa <- emptyH0

  listH0 <- list()
  listHa <- list()

  for(i in seq_len(cases)){
    for(j in seq_len(nSect)){
      simH0[[j]][[i]][["sector"]] <- levels(factSect)[j]
      simHa[[j]][[i]][["sector"]] <- levels(factSect)[j]

      emptyH0 <- plyr::rbind.fill(emptyH0, simH0[[j]][[i]])
      emptyHa <- plyr::rbind.fill(emptyHa, simHa[[j]][[i]])
    }
    listH0[[i]] <- emptyH0
    listHa[[i]] <- emptyHa
    emptyH0 <- data.frame(matrix(nrow = 0, ncol = ncol(data)-2))
    colnames(emptyH0) <- colnames(data[,-c(1,2)])
    emptyHa <- emptyH0
  }

  # separate the lists to set data and labels apart ----
  listH0label <- lapply(listH0, dplyr::select, c("N", "sites", "sector"))[[1]]
  listHalabel <- lapply(listHa, dplyr::select, c("N", "sites", "sector"))[[1]]

  listH0data <- lapply(listH0, dplyr::select, -c("N", "sites", "sector"))
  listHadata <- lapply(listHa, dplyr::select, -c("N", "sites", "sector"))

  # change NA to 0 and remove columns that do not have data ----
  listH0data <- lapply(listH0data, function(x){
    x[is.na(x)] <- 0
    return(x)
  })
  listHadata <- lapply(listHadata, function(x){
    x[is.na(x)] <- 0
    return(x)
  })

  # design and fill the results matrix ----
  NN <- cases * k * (m-1) * (n-1)
  resultsHa <- matrix(nrow = NN, ncol = 8)
  resultsHa[, 1] <- rep(seq(cases), times = 1, each = (k * (m-1) * (n-1)))
  resultsHa[, 2] <- rep(1:k, times = (n-1) * (m-1) * cases)
  resultsHa[, 3] <- rep(seq(2, m), times = (n-1), each = k)
  resultsHa[, 4] <- rep(seq(2, n), times = 1, each = k * (m-1))
  colnames(resultsHa) <- c("dat.sim", "k", "m", "n",
                           "pseudoFH0", "pseudoFHa",
                           "MSAHa", "MSRHa")

  # Loop to calculate pseudoF ----
  listH0 <- lapply(listH0data, cbind, listH0label)
  listHa <- lapply(listHadata, cbind, listHalabel)

  H0Sim <- array(unlist(listH0), dim = c((N * sites * nSect), length(listH0[[1]]), cases))
  HaSim <- array(unlist(listHa), dim = c((N * sites * nSect), length(listHa[[1]]), cases))

  H0names <- colnames(listH0[[1]])
  colnames(H0Sim) <- H0names
  colnames(HaSim) <- H0names

  # lista de 1 a mil
  Y <- cbind(1:(N * sites))
  # lista de 10 sitios repetidos 100 veces
  YPU <- as.numeric(gl(sites, N))
  # etiquetas de los sitios ((m-1) sitios repetidos k veces)
  mm <- resultsHa[,3]
  # números de muestras a tomar (columna m * columna n)
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
                        H0Sim, HaSim, resultsHa,
                        transformation, method, model,
                        nSect, sites, N)
    }
    resultsHa[,5] <- result1[,1]
    resultsHa[,6] <- result1[,2]
    resultsHa[,7] <- result1[,3]
    resultsHa[,8] <- result1[,4]
  } else {
    for (i in seq_len(NN)){
      result1 <- balanced_sampling2(i, Y, mm, nn, YPU,
                                   H0Sim, HaSim, resultsHa,
                                   transformation, method, model,
                                   nSect, sites, N)
      resultsHa[i,5] <- result1[,1]
      resultsHa[i,6] <- result1[,2]
      resultsHa[i,7] <- result1[,3]
      resultsHa[i,8] <- result1[,4]

      setTxtProgressBar(pb, i)
    }
  }

  close(pb)

  # resultsHa <- resultsHa[!is.na(resultsHa[,5] | !is.na(resultsHa[,6])),]
  resultsHa <- resultsHa[!is.na(resultsHa[,5]) & !is.na(resultsHa[,6]),]

  SimResults <- list(Results = resultsHa)
  class(SimResults) <- "ecocbo_data"

  return(SimResults)

}

# # cargar datos de trabajo ----
# library(SSP)
# data("epibionts")
# data <- epibionts
# type <- "counts"
# Sest.method = "average"
# cases = 4
# N = 100
# sites = 10
# n = 5
# m = 6
# k = 30
# transformation = "square root"
# method = "bray"
# dummy = FALSE
# useParallel = TRUE
# model = "nested.symmetric"
