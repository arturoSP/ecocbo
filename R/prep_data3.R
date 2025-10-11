library(tidyr)
library(readr)
library(dplyr)
library(SSP)
library(ecocbo)
# data <- read_csv("pruebas/macrofauna.csv", col_select = -61)|>
#   relocate(locality:campaign, .after = labels)|>
#   filter(locality == "Sisal" & campaign == "3" & site == 1)|>
#   select(-c(labels, locality, site, campaign))|>
#   mutate(across(where(is.numeric), ~ round(.))) |> 
#   as.data.frame()
data <- SSP::epibionts
type = "counts"
Sest.method = "average"
cases = 2
N = 100
M = 15
n = 20
m = 10
k = 30
transformation = "square root"
method = "bray"
dummy = TRUE
useParallel = TRUE
model = "nested.symmetric"
jitter = FALSE

Res1 <- prep_data3(data, type = type, Sest.method = Sest.method,
                      cases = cases, N = N, M = M,
                      n = n, m = m, k = k,
                      transformation = transformation, method = method,
                      dummy = dummy, useParallel = useParallel,
                      model = model,
                      jitter = FALSE)

Res2 <- prep_data3(data, type = type, Sest.method = Sest.method,
                      cases = cases, N = N, M = M,
                      n = n, m = m, k = k,
                      transformation = transformation, method = method,
                      dummy = dummy, useParallel = useParallel,
                      model = model,
                      jitter = TRUE)

Res2$Results |> 
  as.data.frame() |> 
  mutate(cvarBA = (`MSB(A)` - MSR )/ n) |> 
  group_by(m) |> 
  # summarise(CV_BA = mean(cvarBA), CV_R = mean(MSR))
  summarise(CV_BA = max(cvarBA), CV_R = max(MSR))

beta1 <- sim_beta(Res1)
beta2 <- sim_beta(Res2)
scompvar(Res1, n = 2, m = 2)
scompvar(Res2, n = 2, m = 2)
plot_power(beta1, n = 3, m = 4, method = "both")
plot_power(beta2, n = 3, m = 4, method = "both")

prep_data3 <- function(data, type = "counts", Sest.method = "average",
                      cases = 5, N = 100, M = NULL,
                      n, m = NULL, k = 50,
                      transformation = "none", method = "bray",
                      dummy = FALSE, useParallel = TRUE,
                      model = "single.factor", jitter = FALSE){
  # Check the inputs ----

  if (n > N){stop("'n' must be equal or less than 'N' on simulated data")}
  if(ceiling(n) != floor(n)){stop("n must be integer")}
  if(n <= 1){stop("n must be larger than 1")}

  if(model != "single.factor"){
    if (m > M){stop("'m' must be equal or less than 'M' on simulated data")}
    if(ceiling(m) != floor(m)){stop("m must be integer")}
    if(m <= 1){stop("m must be larger than 1")}
  }

  # The function to work with depends on the selected model
  Results <- if(model == "single.factor"){
    prep_data_single3(data, type, Sest.method, cases, N,
                     M, n, m, k, transformation, method,
                     dummy, useParallel)
  } else if(model == "nested.symmetric"){
    prep_data_nestedsymmetric3(data, type, Sest.method, cases, N,
                              M, n, m, k, transformation, method,
                              dummy, useParallel, model, jitter)
  }

  return(Results)
}

prep_data_nestedsymmetric3 <- function(data, type = "counts",
                                      Sest.method = "average",
                                      cases = 5, N = 100, M = 3,
                                      n, m, k = 50,
                                      transformation = "none", method = "bray",
                                      dummy = FALSE, useParallel = TRUE, model, jitter){
  # Temporal change of value for useParallel, it'll be removed when I find the
  # fix for whatever the error is
  # useParallel <- FALSE

  # get values for size limits
  data[,1] <- as.factor(data[,1])
  factSect <- data[,1]
  nSect <- nlevels(factSect)         # number of treatments
  nivel <- levels(factSect)

  # Apply PERMANOVA to the original data, to get an estimated MSB(A)
  # to supplement the simulated MSB(A)
  # vegan::adonis2(data[,-c(1:2)]~ Locality / Site, data = data[,c(1:2)], strata = data$Locality, by = "terms")

  ## Simulated data for Ha ----
  # Create a list to store the results
  ListSimA <- vector(mode = "list", length = nSect)
  names(ListSimA) <- nivel
  ListSim0 <- vector(mode = "list", length = nSect)
  names(ListSim0) <- nivel
  # run the simulation for the different sectors we already have.
  for(i in nivel){
    # Prepare the data by dissecting by treatment
    dataTrimmedA <- data[data[,1] == i, -1]

    # Calculate simulation parameters
    dataParameterA <- SSP::assempar(dataTrimmedA, type = type, Sest.method = Sest.method)
    # Calculate simulated communities
    if(jitter){
      dataSimA <- simdata2(dataParameterA, cases = cases, N = N, sites = M)
    } else {
      dataSimA <- SSP::simdata(dataParameterA, cases = cases, N = N, sites = M)
    }
      

    dataSimA <- lapply(dataSimA, function(df){
      df <- mutate(df, sector = as.factor(i))
      new_names <- gsub("unseen\\.species\\s*(\\d+)",
                        paste0("unseen.species.", i, ".\\1"),
                        names(df))
      names(df) <- new_names
      df
    })

    if(dummy == TRUE){
      dataSimA <- lapply(dataSimA, mutate, dummy = 1)
    }

    # Store the simulations in the list
    ListSimA[[i]] <- dataSimA
  }

  for(i in nivel){
    # Prepare data by setting replicates to just one value
    dataTrimmed0 <- data[, -1]
    # dataTrimmed0[,1] <- "zero"

    # Calculate simulation parameters
    dataParameter0 <- assempar(dataTrimmed0, type = type, Sest.method = Sest.method)
    # Calculate simulated communities
    if(jitter) {
      dataSim0 <- simdata2(dataParameter0, cases = cases, N = N, sites = M)
    } else {
      dataSim0 <- simdata(dataParameter0, cases = cases, N = N, sites = M)
    }
       

    dataSim0 <- lapply(dataSim0, function(df){
      df <- mutate(df, sector = as.factor(i))
      new_names <- gsub("unseen\\.species\\s*(\\d+)",
                        paste0("unseen.species.", "\\1"),
                        names(df))
      names(df) <- new_names
      df
    })

    if(dummy == TRUE){
      dataSim0 <- lapply(dataSim0, mutate, dummy = 1)
    }

    # Store simulations in the list
    ListSim0[[i]] <- dataSim0
  }

  # Organize the data from different simulations in cases
  # Create a list to store results
  simHa <- vector(mode = "list", length = cases)
  simH0 <- vector(mode = "list", length = cases)
  # Using two nested loops we can arrange the data into the simulated cases
  for(j in seq_len(cases)){
    for(i in nivel){
      simHa[[j]] <- bind_rows(as.data.frame(simHa[[j]]),
                              as.data.frame(ListSimA[[i]][j]))
      simH0[[j]] <- bind_rows(as.data.frame(simH0[[j]]),
                              as.data.frame(ListSim0[[i]][j]))
    }
    # Fill NA spaces with 0 and then reorder the columns to have the labels at the end

    if(dummy == TRUE){
      simHa[[j]] <- mutate(simHa[[j]], across(everything(),
                                              ~replace_na(.x, 0))) |>
        relocate(c(starts_with("unseen"), dummy, sector, sites, N),
                 .after = last_col())
      simH0[[j]] <- mutate(simH0[[j]], across(everything(),
                                              ~replace_na(.x, 0))) |>
        relocate(c(starts_with("unseen"), dummy, sector, sites, N),
                 .after = last_col())
    } else {
      simHa[[j]] <- mutate(simHa[[j]], across(everything(),
                                              ~replace_na(.x, 0))) |>
        relocate(c(starts_with("unseen"), sector, sites, N),
                 .after = last_col())
      simH0[[j]] <- mutate(simH0[[j]], across(everything(),
                                              ~replace_na(.x, 0))) |>
        relocate(c(starts_with("unseen"), sector, sites, N),
                 .after = last_col())
    }
  }

  ## design and fill the results matrix ----
  NN <- cases * k * (m-1) * (n-1)
  resultsHa <- matrix(nrow = NN, ncol = 8)
  colnames(resultsHa) <- c("dat.sim", "k", "m", "n",
                           "pseudoFH0", "pseudoFHa",
                           "MSB(A)", "MSR")

  resultsHa[, 1] <- rep(seq(cases), times = 1, each = (k * (m-1) * (n-1)))
  resultsHa[, 2] <- rep(1:k, times = (n-1) * (m-1) * cases)
  resultsHa[, 3] <- rep(seq(2, m), times = (n-1) * cases, each = k)
  resultsHa[, 4] <- rep(seq(2, n), times = cases, each = k * (m-1))

  ## design the arrays to store the lists ----
  # H0 comes from simdata as-is, it does not include a column for sectors given that
  # its data was simulated by merging sectors and sites together
  H0Sim <- array(unlist(simH0), dim = c(dim(simH0[[1]])[1], length(simH0[[1]]), cases))
  colnames(H0Sim) <- colnames(simH0[[1]])
  H0Sim <- H0Sim[,1:(dim(H0Sim)[2]-3),]

  HaSim <- array(unlist(simHa), dim = c(dim(simHa[[1]])[1], length(simHa[[1]]), cases))
  colnames(HaSim) <- colnames(simHa[[1]])
  HaSim <- HaSim[,1:(dim(HaSim)[2]-3),]

  # Factors matrix and data matrix for Ha
  factEnv <- as.data.frame(simHa[[1]][,(dim(simHa[[1]])[2]-2):dim(simHa[[1]])[2]])
  names(factEnv) <- c("sector", "site", "N")

  rm(simH0, simHa)

  ## Set of parameters for using balancedtwostage ----
  # index marking the size of each resampled site
  Y <- cbind(1:(N * M))
  # index for the sites repeated N times
  YPU <- as.numeric(gl(M, N))
  # labels for sites (i.e. (m-1) sites, repeated k times)
  mm <- resultsHa[,3]
  # number of samples to consider (i.e. m * n)
  nn <- resultsHa[,3] * resultsHa[,4]

  # Se crean dos dataframes para el muestreo
  Y1 <- cbind(Y, YPU)
  mn <- cbind(mm, nn)

  if(useParallel){
    # Registering the cluster of workers with parabar
    parabar::configure_bar(type = "modern", format = "[:bar] :percent")
    cl <- parabar::start_backend(cores = parallelly::availableCores()/2,
                                 cluster_type="psock",
                                 backend_type="async")

    # Exporing functions needed for the parallel iterations
    parabar::export(cl,
                    variables = c("balanced_sampling23", "permanova_twoway3", "SS", 
                    "dbmanova_nested"),
                    # environment = asNamespace("ecocbo")
                    )

    # Executing the loop in parallel
    result1 <- parabar::par_lapply(cl, x = 1:NN, fun= balanced_sampling23,
                                   NN,
                                   Y1, mn, nSect, M, N, H0Sim, HaSim, resultsHa,
                                   factEnv, transformation, method, model) |>
      unlist() |>
      matrix(ncol=4, byrow=TRUE)
    colnames(result1) <- c("FobsH0", "FobsHa", "MSBA", "MSR")

    # Assigning the results to the outcome matrix
    resultsHa[,5:8] <- result1[,c(1,2,3,4)]
    parabar::stop_backend(cl)
    rm(result1)
  } else {
    # progress bar
    pb <- txtProgressBar(max = NN, style = 3)

    for (i in seq_len(NN)){
      result1 <- balanced_sampling23(i, NN, Y1, mn, nSect, M, N, H0Sim, HaSim,
                                    resultsHa, factEnv, transformation,
                                    method, model)
      resultsHa[i,5] <- result1[1]
      resultsHa[i,6] <- result1[2]
      resultsHa[i,7] <- result1[3]
      resultsHa[i,8] <- result1[4]
      setTxtProgressBar(pb, i)
    }
    rm(result1)
    close(pb)
  }

  resultsHa <- resultsHa[!is.na(resultsHa[,5]) & !is.na(resultsHa[,6]),]

  SimResults <- list(Results = resultsHa, model = model, a = nSect)
  class(SimResults) <- c("ecocbo_data", class(SimResults))

  return(SimResults)

}

permanova_twoway3 <- function(x, factEnv, method = "bray",
                             transformation = "none", model = "nested.symmetric"){

  pseudoF_2NestedSymmetric <- function(x, factEnv, method, transformation){
    # Apply transformation and calculate distance matrix
    if (transformation == "square root") {
      x.t <- sqrt(x)
      d <- vegan::vegdist(x.t, method = method)
    } else if (transformation == "fourth root") {
      x.t <- sqrt(sqrt(x))
      d <- vegan::vegdist(x.t, method = method)
    } else if (transformation == "Log (X+1)") {
      x.t <- log(x + 1)
      d <- vegan::vegdist(x.t, method = method)
    } else if (transformation == "P/A") {
      x.t <- 1 * (x > 0)
      d <- vegan::vegdist(x.t, method = method, binary = TRUE)
    } else {
      x.t <- x
      d <- vegan::vegdist(x.t, method = method)
    }
    # rm(x)

    # size for the labels we work with
    a = nlevels(as.factor(factEnv$sector))  # number of treatments (A)
    b = length(unique(as.factor(factEnv$site)))   # number of replicates (B)
    factEnv["secsit"] <- paste0(factEnv$sector, factEnv$site) # intersections AB
    nBA = nlevels(as.factor(factEnv$secsit))  # number of intersections AB
    nRep = dim(factEnv)[1] / nBA  # number of times we're repeating each intersection
    nNm = unique(factEnv$sector) # unique values for the sectors
    nScSt = unique(factEnv$secsit)  # unique values for the intersections site-sector

    # calculates SS for all
    # SST <- SS(d*100)[2]  # this *100 is added to get results that are comparable to those in PRIMER
    SST <- SS(d)[2]

    # calculates SS within replicates
    secsite_groups <- split(rownames(factEnv), factEnv$secsit)
    listR <- sapply(secsite_groups, function(rw) {
      SS(vegan::vegdist(x.t[rw,], method = method))
    }, simplify = "array")
    SSR <- sum(listR[2,])
    # rm(secsite_groups, listR)

    # calculates SS_B(A)

    # Bray-Curtis distance matrix is calculated for each sector, then the centroids
    # are estimated using `betadisper()`. Negative values are taken out, and then
    # the SS for the Euclidian distance matrix is calculated. The SS are summed
    # to calculate SS_B(A)
    sector_groups <- split(rownames(factEnv), factEnv$sector)
    # sector_groups <- lapply(sector_groups, as.numeric)

    listBA <- sapply(sector_groups, function(rw) {
      dBA <- vegan::vegdist(x.t[rw,], method = "bray")
      tCentroid <- vegan::betadisper(dBA,
                                     group = factEnv[rw, "secsit"],
                                     type = "centroid",
                                     bias.adjust = FALSE)
      Eig <- which(tCentroid$eig > 0)
      SS(vegan::vegdist(tCentroid$centroids[, Eig, drop = FALSE],
                        method = "euclidean"))
    }, simplify = "array")
    SSBA <- sum(listBA[2,]) * nRep  # * nRep is added so that when calculating
    # the variation between nested groups it is weighted by the size of the subgroup
    # rm(sector_groups, listBA)

    # calculates SSA
    SSA <- SST - SSBA - SSR

    # fill the permanova table
    # degrees of freedom
    DoFA <- a - 1
    DoFBA <- a * (b - 1)
    DoFR <- a * b * (nRep - 1)
    DoFT <- (a * b * nRep) - 1

    # mean squares

    MSA <- (SSA / DoFA)
    MSBA <- (SSBA / DoFBA)
    MSR <- SSR / DoFR

    # observed pseudoF
    FobsA <- MSA / MSBA
    FobsBA <- MSBA / MSR

    Fobs <- as.data.frame(matrix(nrow = 4, ncol = 4))
    colnames(Fobs) <- c("SS", "DoF", "MS", "F")
    rownames(Fobs) <- c("A", "B(A)", "R", "T")
    Fobs[1,1] <- SSA
    Fobs[1,2] <- DoFA
    Fobs[1,3] <- MSA
    Fobs[1,4] <- FobsA
    Fobs[2,1] <- SSBA
    Fobs[2,2] <- DoFBA
    Fobs[2,3] <- MSBA
    Fobs[2,4] <- FobsBA
    Fobs[3,1] <- SSR
    Fobs[3,2] <- DoFR
    Fobs[3,3] <- MSR
    Fobs[4,1] <- SST
    Fobs[4,2] <- DoFT

    # return(Fobs)
    Ffinal <- c(Fobs[1,4], Fobs[1,3], Fobs[2,3], Fobs[3,3])
    return(Ffinal)

  }

  ## Main function which returns pseudoF and MS ----
  if(dim(x)[1] != dim(factEnv)[1])
    stop("x and factEnv dimensions do not match")
  if(model != "nested.symmetric" & model != "nested.asymmetric" & model != "orthogonal")
    stop("Possible values for type are \"nested.symmetric\", \"nested.asymmetric\",
         \"orthogonal\" and \"single.factor\"")

  if(model == "nested.symmetric"){
    Results <- pseudoF_2NestedSymmetric(x, factEnv, method, transformation)
  } else if(model == "nested.asymmetric"){
    Results <- pseudoF_2NestedSymmetric(x, factEnv, method, transformation)
  } else {
    Results <- pseudoF_2Orthogonal(x, factEnv, method, transformation)
  }

  Fobs <- Results
  return(Fobs)
}

balanced_sampling23 <- function(i, NN, Y1, mn, nSect, M, N, H0Sim, HaSim,
                               resultsHa, factEnv, transformation, method,
                               model){
  # Determine index for sampling units
  indice <- sampling::balancedtwostage(as.matrix(Y1[,1]), selection = 1, m = mn[i,1],
                                       n = mn[i,2], PU = Y1[,2], comment = FALSE)[,1] |>
    unname() |>
    as.logical() |>
    which()

  ones_n <- rep(indice, nSect)
  ones_s <- rep(c(0:(nSect-1)) * M * N, each = length(indice))
  ones <- ones_n + ones_s
  rm(ones_n, ones_s, indice)

  # Extract samples from the datasets and evaluate with PERMANOVA
  y0 <- H0Sim[ones, , resultsHa[i,1]]
  rownames(y0) <- ones
  ya <- HaSim[ones, , resultsHa[i,1]]
  rownames(ya) <- ones
  factEnvX <- factEnv[ones,]

  result_0 <- dbmanova_nested(x = y0,
                               factEnv = factEnvX,
                               transformation=transformation,
                               method = method,
                               model = model)
  result_a <- dbmanova_nested(x = ya,
                               factEnv = factEnvX,
                               transformation=transformation,
                               method = method,
                               model = model)

  #if(length(result_0) != 4){result_0 <- c(NA, NA, NA, NA)}
  #if(length(result_a) != 4){result_a <- c(NA, NA, NA, NA)}
  # Assemble results with pseudoF for H0 and Ha and MSBA + MSR for Ha
  result1 <- c(result_0[1,5], result_a[1,5],
               result_a[2,4], result_a[3,4])


  return(result1)
}

SS <- function (d) {
  ss <- numeric(2)
  ss[1] <- dim(as.matrix(d))[1]
  ss[2] <- sum(d^2)/ss[1]
  return(ss)
}

label_permanova <- function(dataP, factEnvP, method, transformation,
                            dummy, model){

  if(dummy){
    dataP$dummy = 1
  }

  # Apply transformation and calculate distance matrix
  if (transformation == "square root") {
    x.t <- sqrt(dataP)
    d <- vegan::vegdist(x.t, method = method)
  } else if (transformation == "fourth root") {
    x.t <- sqrt(sqrt(dataP))
    d <- vegan::vegdist(x.t, method = method)
  } else if (transformation == "Log (X+1)") {
    x.t <- log(dataP + 1)
    d <- vegan::vegdist(x.t, method = method)
  } else if (transformation == "P/A") {
    x.t <- 1 * (dataP > 0)
    d <- vegan::vegdist(x.t, method = method, binary = TRUE)
  } else {
    x.t <- dataP
    d <- vegan::vegdist(x.t, method = method)
  }
  # rm(dataP)

  # Compute the number of permutations available to the experiment,
  # then compare it with the given k
  a = nlevels(as.factor(factEnvP[,1]))  # number of treatments (A)
  b = length(unique(as.factor(factEnvP[,2])))   # number of replicates (B)
  factEnvP["secsit"] <- paste0(factEnvP[,1], factEnvP[,2]) # intersections AB
  nBA = nlevels(as.factor(factEnvP$secsit))  # number of intersections AB
  nRep = dim(factEnvP)[1] / nBA  # number of times we're repeating each intersection
  nNm = unique(factEnvP[,1]) # unique values for the sectors
  nScSt = unique(factEnvP$secsit)  # unique values for the intersections site-sector

  # calculates SS for all
  SST <- SS(d)[2]

  # permList <- vector("list", 999 + 1)
  # degrees of freedom
  DoFA <- a - 1
  DoFBA <- a * (b - 1)
  DoFR <- a * b * (nRep - 1)
  DoFT <- (a * b * nRep) - 1

  # Compute ANOVA for the original data
  secsite_groups <- split(rownames(factEnvP), factEnvP$secsit)
  listR <- sapply(secsite_groups, function(rw) {
    SS(vegan::vegdist(x.t[rw,], method = method))
  }, simplify = "array")
  SSR <- sum(listR[2,])
  # calculates SS_B(A)
  sector_groups <- split(rownames(factEnvP), factEnvP[,1])

  listBA <- sapply(sector_groups, function(rw) {
    dBA <- vegan::vegdist(x.t[rw,], method = "bray")
    tCentroid <- vegan::betadisper(dBA,
                                   group = factEnvP[rw, "secsit"],
                                   type = "centroid",
                                   bias.adjust = FALSE)
    Eig <- which(tCentroid$eig > 0)
    SS(vegan::vegdist(tCentroid$centroids[, Eig, drop = FALSE],
                      method = "euclidean"))
  }, simplify = "array")
  SSBA <- sum(listBA[2,]) * nRep

  # calculates SSA
  SSA <- SST - SSBA - SSR

  # fill the permanova table
  # mean squares
  MSA <- SSA / DoFA
  MSBA <- SSBA / DoFBA
  MSBA <- (SSBA / DoFBA)
  MSR <- SSR / DoFR

  # observed pseudoF
  FobsA <- MSA / MSBA
  FobsBA <- MSBA / MSR

  Fobs <- as.data.frame(matrix(nrow = 4, ncol = 4))
  colnames(Fobs) <- c("SS", "DoF", "MS", "F")
  rownames(Fobs) <- c("A", "B(A)", "R", "T")
  Fobs[1,1] <- SSA
  Fobs[1,2] <- DoFA
  Fobs[1,3] <- MSA
  Fobs[1,4] <- FobsA
  Fobs[2,1] <- SSBA
  Fobs[2,2] <- DoFBA
  Fobs[2,3] <- MSBA
  Fobs[2,4] <- FobsBA
  Fobs[3,1] <- SSR
  Fobs[3,2] <- DoFR
  Fobs[3,3] <- MSR
  Fobs[4,1] <- SST
  Fobs[4,2] <- DoFT

  # If this were full PERMANOVA, it's necessary to change the return to permList
  return(Fobs)
}

