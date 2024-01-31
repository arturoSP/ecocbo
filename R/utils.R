#' Sum of squares using Huygen Theorem
#'
#' Calculates sum of squares using Huygen theorem as implemented by Anderson (2014).
#'
#' @param d distance matrix from which the sum of squares will be calculated.
#'
#' @return A numeric vector containing the dimension for the distance matrix, and
#' the value for the sum of squares for the matrix.
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @references Anderson, M. J. (2014). Permutational multivariate analysis of
#' variance (PERMANOVA). Wiley statsref: statistics reference online, 1-15.
#'
#' @export
#' @keywords internal
#'

SS <- function (d) {
  ss <- numeric(2)
  ss[1] <- dim(as.matrix(d))[1]
  ss[2] <- sum(d^2)/ss[1]
  return(ss)
}

#' PERMANOVA one-way
#'
#' Calculates observed F and mean squares for the residuals and among sites. This
#' function is a helper for [prep_data()].
#'
#' @param x ecological community data.
#' @param factEnv label for the community data.
#' @param type which algorithm to use for the calculation? At the moment, the only
#' option is "P".
#' @param method appropriate distance/dissimilarity metric (e.g. Gower,
#' Bray–Curtis, Jaccard, etc).
#' @param transformation Mathematical function to reduce the weight of very
#' dominant species.
#'
#' @return A data frame containing the resulting PERMANOVA table.
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @references Underwood, A. J. (1997). Experiments in ecology: their logical
#' design and interpretation using analysis of variance. Cambridge university
#' press.
#' @references Anderson, M. J. (2014). Permutational multivariate analysis of
#' variance (PERMANOVA). Wiley statsref: statistics reference online, 1-15.
#'
#' @importFrom vegan vegdist
#'
#' @seealso [vegan::vegdist()]
#'
#' @export
#' @keywords internal

## PERMANOVA ----
permanova_oneway <- function(x, factEnv, type = "P", method = "bray", transformation = "none"){
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

    SST <- SS(d)[2]

    # Size for labels
    nlev <- nlevels(as.factor(factEnv))

    # Calculate the SS for residuals
    lista <- split(as.data.frame(x.t), factEnv)
    dimxt <- dim(x.t)
    vdist <- lapply(lista, vegan::vegdist, method = method)
    SSi <- lapply(vdist, SS)
    SSi <- array(unlist(SSi), dim = c(1,2,nlev))

    # Calculate denominators
    denR <- dimxt[1] - nlev
    denA <- nlev - 1

    # Results
    SSR <- sum(SSi[,2,])
    SSA <- abs(SST - SSR)
    MSR <- (SSR/denR)
    MSA <- (SSA/denA)
    Fobs <- MSA/MSR
    Fobs <- data.frame(SSA, SSR, SST, denA, denR, MSA, MSR, Fobs)
    return(Fobs)
  }

  # Prueba Brown-Forsythe ---
  # está desactivada, pero guardada por si se quiere usar algún día
  # pseudoF_BF <- function(x, factEnv, method = "bray"){
  #   d <- vegan::vegdist(x, method = method)
  #   SST <- SS(d)[2]
  #
  #   # Size for labels
  #   lev <- table(factEnv)
  #   nlev <- length(lev)
  #
  #   # Empty helper vectors
  #   Var <- numeric(nlev)
  #   d.res <- Var
  #
  #   # Calculate SS and Var for residuals
  #   lista <- split(x, factEnv)
  #   vdist <- lapply(lista, vegan::vegdist)
  #   SSi <- lapply(vdist, SS)
  #   SSi <- array(unlist(SSi), dim = c(1,2,nlev))
  #
  #   Var <- SSi[,2,]/(SSi[,1,]-1)
  #   d.res <- (1-(lev/sum(lev)))*Var
  #   den <- sum(d.res)
  #
  #   # Resultados
  #   SSR <- sum(SSi[,2,])
  #   SSA <- abs(SST-SSR)
  #   Fobs<- SSA/den
  #
  #   Fobs <- data.frame(SSA, SSR, SST, den, Fobs)
  #   return(Fobs)
  # }

  # Main function which returns pseudoF and MS ----
  # Validating data
  if(dim(x)[1] != length(factEnv))
    stop("x and factEnv dimensions do not match")
  if(type != "P" & type != "BF")
    stop("Possible values for type are P and BF")

  Results <- pseudoF_P(x, factEnv, method, transformation)[1,c(6,7,8)]
  Fobs <- Results
  return(Fobs)
}

#' Balanced sampling
#'
#' Develops the experimental design based on the provided conditions
#'
#' @param i pointer to the index in the list of experimental designs to try.
#' @param Y index to the data.frame the function will work with.
#' @param mm number of site the function is working with in each iteration.
#' @param nn number of samples to consider in each iteration.
#' @param YPU label for the sites in each iteration, as used by
#' [sampling::balancedtwostage()]
#' @param H0Sim simulated community from \code{SSP::simdata()} in which H0 is
#' true.
#' @param HaSim simulated community from \code{SSP::simdata()} in which H0 is
#' false.
#' @param resultsHa helper matrix that stores labels and later the results.
#' @param method appropriate distance/dissimilarity metric (e.g. Gower,
#' Bray–Curtis, Jaccard, etc).
#' @param transformation Mathematical function to reduce the weight of very
#' dominant species.
#'
#' @return a data frame with values for observed F (for H0 and Ha), and the Ha mean
#' squares for residuals and variation among sites.
#'
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @references Underwood, A. J. (1997). Experiments in ecology: their logical
#' design and interpretation using analysis of variance. Cambridge university
#' press.
#'
#' @importFrom sampling balancedtwostage
#'
#' @seealso [sampling::balancedtwostage()]
#'
#' @export
#' @keywords internal

balanced_sampling <- function(i, Y, mm, nn, YPU, H0Sim, HaSim, resultsHa, transformation, method){
  # Get the samples index
  sel <- sampling::balancedtwostage(Y, selection = 1, m = mm[i],
                                    n = nn[i], PU = YPU, comment = FALSE)
  ones <- which(sel[,1] %in% 1)
  y0 <- H0Sim[ones,,resultsHa[i,1]]
  ya <- HaSim[ones,,resultsHa[i,1]]
  yHa <- dim(y0)[2] - 2

  # Apply PERMANOVA to get F and mean squares
  result1 <- permanova_oneway(x = y0[, 1:yHa], factEnv = y0[,yHa+2],
                              transformation = transformation, method = method)
  result2 <- permanova_oneway(x = ya[, 1:yHa], factEnv = y0[,(yHa+2)],
                              transformation = transformation, method = method)
  result0 <- matrix(nrow = 1, ncol = 4)
  colnames(result0) <- c("Fobs", "Fobs", "MSA", "MSR")

  # Gather the results and return
  result0[,1] <- result1[,3]
  result0[,2] <- result2[,3]
  result0[,3] <- result2[,1]
  result0[,4] <- result2[,2]
  return(result0)
}

#' PERMANOVA two-way
#'
#' Calculates observed F and mean squares for the residuals and among sites. This
#' function is a helper for [prep_data2()].
#'
#' @param x ecological community data.
#' @param factEnv label for the community data.
#' @param model which algorithm to use for the calculation? At the moment, the only
#' option is "nested.symmetric".
#' @param method appropriate distance/dissimilarity metric (e.g. Gower,
#' Bray–Curtis, Jaccard, etc).
#' @param transformation Mathematical function to reduce the weight of very
#' dominant species.
#'
#' @return A data frame containing the resulting PERMANOVA table.
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @references Underwood, A. J. (1997). Experiments in ecology: their logical
#' design and interpretation using analysis of variance. Cambridge university
#' press.
#' @references Anderson, M. J. (2014). Permutational multivariate analysis of
#' variance (PERMANOVA). Wiley statsref: statistics reference online, 1-15.
#'
#' @importFrom vegan vegdist
#' @importFrom vegan betadisper
#'
#' @seealso [vegan::vegdist()]
#'
#' @export
#' @keywords internal

## PERMANOVA Two factors ----
permanova_twoway <- function(x, factEnv, method = "bray", transformation = "none", model = "nested.symmetric"){

  pseudoF_2Orthogonal <- function(x, factEnv, method = "bray", transformation = "none"){
    # data transformation, if necessary, and calculate the distance matrix
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

    # size for the labels we work with
    a = nlevels(as.factor(factEnv$sector)) # number of sectors (A)
    b = nlevels(as.factor(factEnv$sites))  # number of sites (B)
    factEnv["secsit"] <- paste0(factEnv$sector, factEnv$site)  # labels for the intersection (AB)
    n = nlevels(as.factor(factEnv$secsit)) # number of intersections (AB)

    # sum of squares total, Huygen for all
    SST <- SS(d)[2]

    # calculates SS for A
    listA <- list()
    for(i in levels(as.factor(factEnv$sector))){
      RwNm <- rownames(factEnv[factEnv$sector == i,])
      listA[[i]] <- SS(vegdist(x.t[RwNm,], method = method))
    }

    listA <- array(unlist(listA), dim = c(1,2,a))
    SSdA <- sum(listA[,2,])
    SSA <- SST - SSdA

    # calculate SS for B
    listB <- list()
    for(i in levels(as.factor(factEnv$sites))){
      RwNm <- rownames(factEnv[factEnv$sites == i,])
      listB[[i]] <- SS(vegdist(x.t[RwNm,], method = method))
    }

    listB <- array(unlist(listB), dim = c(1,2,b))
    SSdB <- sum(listB[,2,])
    SSB <- SST - SSdB

    # calculate SS for residuals
    listR <- list()
    for(i in levels(as.factor(factEnv$secsit))){
      RwNm <- rownames(factEnv[factEnv$secsit == i,])
      listR[[i]] <- SS(vegdist(x.t[RwNm,], method = method))
    }

    listR <- array(unlist(listR), dim = c(1,2,n))
    SSR <- sum(listR[,2,])

    # calculate SS for interaction A-B
    SSAB <- SST - SSA - SSB - SSR

    # fill the permanova table
    # degrees of freedom
    DoFA <- a - 1
    DoFB <- b - 1
    DoFAB <- DoFA * DoFB
    DoFR <- a * b * (n - 1)
    DoFT <- (a * b * n) - 1

    # mean squares
    MSA <- SSA / DoFA
    MSB <- SSB / DoFB
    MSAB <- SSAB / DoFR
    MSR <- SSR / DoFT

    # observed pseudoF
    # these divisions depend on the type of factors we're dealing with.
    # it's best explained in Table 10.8 from (Underwood, 1997)
    FobsA <- MSA / MSAB
    FobsB <- MSB / MSAB
    FobsAB <- MSAB / MSR

    Fobs <- as.data.frame(matrix(nrow = 5, ncol = 4))
    colnames(Fobs) <- c("SS", "DoF", "MS", "F")
    rownames(Fobs) <- c("A", "B", "AB", "R", "T")
    Fobs[1,1] <- SSA
    Fobs[1,2] <- DoFA
    Fobs[1,3] <- MSA
    Fobs[1,4] <- FobsA
    Fobs[2,1] <- SSB
    Fobs[2,2] <- DoFB
    Fobs[2,3] <- MSB
    Fobs[2,4] <- FobsB
    Fobs[3,1] <- SSAB
    Fobs[3,2] <- DoFAB
    Fobs[3,3] <- MSAB
    Fobs[3,4] <- FobsAB
    Fobs[4,1] <- SSR
    Fobs[4,2] <- DoFR
    Fobs[4,3] <- MSR
    Fobs[5,1] <- SST
    Fobs[5,2] <- DoFT

    return(Fobs)
  }

  pseudoF_2NestedSymmetric <- function(x, factEnv, method, transformation){
    # data transformation, if necessary, and calculate the distance matrix
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

    # size for the labels we work with
    a = nlevels(as.factor(factEnv$sector))  # number of sectors (A)
    b = nlevels(as.factor(factEnv$sites))   # number of sites (B)
    factEnv["secsit"] <- paste0(factEnv$sector, factEnv$site) # intersections AB
    nBA = nlevels(as.factor(factEnv$secsit))  # number of intersections AB
    nRep = dim(factEnv)[1] / nBA  # number of times we're repeating each intersection
    nNm = unique(factEnv$sector)
    nScSt = unique(factEnv$secsit)

    # calculates SS for all
    # SST <- SS(d*100)[2]
    SST <- SS(d)[2]

    # calculates SS within replicates
    listR <- list()
    for(i in nScSt){
      RwNm <- rownames(factEnv[factEnv$secsit == i,])
      dR <- vegdist(x.t[RwNm,], method = method)
      # listR[[i]] <- SS(dR*100)
      listR[[i]] <- SS(dR)
    }

    listR <- array(unlist(listR), dim = c(1,2,nBA))
    SSR <- sum(listR[,2,])

    # calculates SS_B(A)

    # se calcula distancia bray-curtis para cada sector, se calculan centroides para cada sector
    # se eliminan componentes principales negativos, se calcula SS de la matriz de distancia
    # euclidiana de los centroides, se suman los SS para obtener SSBA
    listBA <- list()
    for(i in nNm){
      RwNm <- rownames(factEnv[factEnv$sector == i,])
      dBA <- vegan::vegdist(x.t[RwNm,], method = method)
      tCentroid <- vegan::betadisper(dBA, # * 100,
                                     group = factEnv[factEnv[,3] == i,4],
                                     type = "centroid",
                                     bias.adjust = F)
      Eig <- which(tCentroid$eig > 0)
      listBA[[i]] <- vegan::vegdist(tCentroid$centroids[,Eig],
                                    method = "euclidean")
    }
    listBA <- lapply(listBA, SS)
    listBA <- array(unlist(listBA), dim = c(1,2,a))
    SSBA <- sum(listBA[,2,]) * nRep # this *nRep product was added to match the results in PRIMER

    # calculates SSA
    SSA <- SST - SSBA - SSR

    # fill the permanova table
    # degrees of freedom
    DoFA <- a - 1
    DoFBA <- a * (b - 1)
    DoFR <- a * b * (nRep - 1)
    DoFT <- (a * b * nRep) - 1

    # mean squares
    MSA <- SSA / DoFA
    MSBA <- SSBA / DoFBA
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

    return(Fobs)

  }

  ## Main function which returns pseudoF and MS ----
  if(dim(x)[1] != dim(factEnv)[1])
    stop("x and factEnv dimensions do not match")
  if(model != "nested.symmetric" & model != "nested.asymmetric" & model != "orthogonal")
    stop("Possible values for type are \"nested.symmetric\", \"nested.asymmetric\",
         \"orthogonal\" and \"random.oneway\"")

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

#' Balanced sampling 2
#'
#' Develops the experimental design based on the provided conditions
#'
#' @param i pointer to the index in the list of experimental designs to try.
#' @param Y index to the data.frame the function will work with.
#' @param mm number of site the function is working with in each iteration.
#' @param nn number of samples to consider in each iteration.
#' @param YPU label for the sites in each iteration, as used by
#' [sampling::balancedtwostage()]
#' @param H0Sim simulated community from \code{SSP::simdata()} in which H0 is
#' true.
#' @param HaSim simulated community from \code{SSP::simdata()} in which H0 is
#' false.
#' @param resultsHa helper matrix that stores labels and later the results.
#' @param transformation Mathematical function to reduce the weight of very
#' dominant species.
#' @param method appropriate distance/dissimilarity metric (e.g. Gower,
#' Bray–Curtis, Jaccard, etc).
#' @param model which algorithm to use for the calculation? At the moment, the only
#' option is "nested.symmetric".
#' @param nSect Total number of sectors to be simulated in each data set.
#' @param sites Total number of sites to be simulated in each data set.
#' @param N Total number of samples to be simulated in each site.
#'
#' @return a data frame with values for observed F (for H0 and Ha), and the Ha mean
#' squares for residuals and variation among sites.
#'
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @references Underwood, A. J. (1997). Experiments in ecology: their logical
#' design and interpretation using analysis of variance. Cambridge university
#' press.
#'
#' @importFrom sampling balancedtwostage
#'
#' @seealso [sampling::balancedtwostage()]
#'
#' @export
#' @keywords internal
#'

balanced_sampling2 <- function(i, Y, mm, nn, YPU, H0Sim, HaSim, factEnv, resultsHa,
                               transformation, method, model, nSect, sites, N){
  # Get the samples index
  sel <- sampling::balancedtwostage(Y, selection = 1, m = mm[i],
                                    n = nn[i], PU = YPU, comment = FALSE)
  ones <- which(sel[,1] %in% 1)
  ones_n <- rep(ones, nSect)
  ones_s <- rep(c(0:(nSect-1)) * sites * N, each = length(ones))
  ones <- ones_n + ones_s

  y0 <- H0Sim[ones,,resultsHa[i,1]]
  rownames(y0) <- ones
  ya <- HaSim[ones,,resultsHa[i,1]]
  rownames(ya) <- ones
  factEnv <- factEnv[ones,]

  # Apply PERMANOVA to get F and mean squares
  result1 <- permanova_twoway(x = y0,
                              factEnv = factEnv,
                              transformation = transformation,
                              method = method,
                              model = model)
  result2 <- permanova_twoway(x = ya,
                              factEnv = factEnv,
                              transformation = transformation,
                              method = method,
                              model = model)
  result0 <- matrix(nrow = 1, ncol = 5)

  # Values of pseudoF for A are stored in the result, values for MSA and MSR come from the dataset with Ha
  colnames(result0) <- c("FobsH0", "FobsHa", "MSA", "MSBA", "MSR")

  # Gather the results and return
  result0[,1] <- result1["A", "F"]
  result0[,2] <- result2["A", "F"]
  result0[,3] <- result2["A", "MS"]
  result0[,4] <- result2["B(A)", "MS"]
  result0[,5] <- result2["R", "MS"]
  return(result0)

}

#' Balanced sampling 3
#'
#' Develops the experimental design based on the provided conditions
#'
#' @param i pointer to the index in the list of experimental designs to try.
#' @param Y index to the data.frame the function will work with.
#' @param mm number of site the function is working with in each iteration.
#' @param nn number of samples to consider in each iteration.
#' @param YPU label for the sites in each iteration, as used by
#' [sampling::balancedtwostage()]
#' @param H0Sim simulated community from \code{SSP::simdata()} in which H0 is
#' true.
#' @param HaSim simulated community from \code{SSP::simdata()} in which H0 is
#' false.
#' @param resultsHa helper matrix that stores labels and later the results.
#' @param transformation Mathematical function to reduce the weight of very
#' dominant species.
#' @param method appropriate distance/dissimilarity metric (e.g. Gower,
#' Bray–Curtis, Jaccard, etc).
#' @param model which algorithm to use for the calculation? At the moment, the only
#' option is "nested.symmetric".
#' @param nSect Total number of sectors to be simulated in each data set.
#' @param sites Total number of sites to be simulated in each data set.
#' @param N Total number of samples to be simulated in each site.
#'
#' @return a data frame with values for observed F (for H0 and Ha), and the Ha mean
#' squares for residuals and variation among sites.
#'
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @references Underwood, A. J. (1997). Experiments in ecology: their logical
#' design and interpretation using analysis of variance. Cambridge university
#' press.
#'
#' @importFrom sampling balancedtwostage
#'
#' @seealso [sampling::balancedtwostage()]
#'
#' @export
#' @keywords internal
#'

# balanced_sampling3 <- function(i, Y, mm, nn, YPU, H0Sim, HaSim, factEnv, resultsHa,
#                                transformation, method, model, nSect, sites, N){
#   # Get the samples index
#   sel <- sampling::balancedtwostage(Y, selection = 1, m = mm[i],
#                                     n = nn[i], PU = YPU, comment = FALSE)
#   ones <- which(sel[,1] %in% 1)
#   ones_n <- rep(ones, nSect)
#   ones_s <- rep(c(0:(nSect-1)) * sites * N, each = length(ones))
#   ones <- ones_n + ones_s
#
#   y0 <- H0Sim[[resultsHa[i,1]]][ones,]
#   rownames(y0) <- ones
#   ya <- HaSim[[resultsHa[i,1]]][ones,]
#   rownames(ya) <- ones
#   factEnv <- factEnv[ones,]
#
#   # Apply PERMANOVA to get F and mean squares
#   result1 <- permanova_twoway(x = y0,
#                               factEnv = factEnv,
#                               transformation = transformation,
#                               method = method,
#                               model = model)
#   result2 <- permanova_twoway(x = ya,
#                               factEnv = factEnv,
#                               transformation = transformation,
#                               method = method,
#                               model = model)
#   result0 <- matrix(nrow = 1, ncol = 5)
#   # Values of pseudoF for A are stored in the result, values for MSA and MSR come from the
#   # dataset with Ha
#   colnames(result0) <- c("FobsH0", "FobsHa", "MSA", "MSBA", "MSR")
#
#   # Gather the results and return
#   result0[,1] <- result1["A", "F"]
#   result0[,2] <- result2["A", "F"]
#   result0[,3] <- result2["A", "MS"]
#   result0[,4] <- result2["B(A)", "MS"]
#   result0[,5] <- result2["R", "MS"]
#   return(result0)
#
# }



