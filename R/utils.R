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
#' @param transformation Mathematical function to reduce the weight of
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
    # Validate transformation method
    valid_transformations <- c("none", "square root", "fourth root", "Log (X+1)", "P/A")
    if (!transformation %in% valid_transformations) {
      stop("Invalid transformation method. Choose from: none, square root, fourth root, Log (X+1), P/A.")
    }

    # Validate distance method
    valid_methods <- c("euclidean", "bray", "jaccard", "gower")
    if (!method %in% valid_methods) {
      stop("Invalid distance method. Choose from: euclidean, bray, jaccard, gower.")
    }

    # Ensure dimensions match
    if (dim(x)[1] != length(factEnv)) {
      stop("The dimensions of the data matrix and the factor labels do not match.")
    }

    # Apply transformation and calculate distance matrix
    if (transformation == "square root") {
      x.t <- sqrt(x)
    } else if (transformation == "fourth root") {
      x.t <- sqrt(sqrt(x))
    } else if (transformation == "Log (X+1)") {
      x.t <- log(x + 1)
    } else if (transformation == "P/A") {
      x.t <- 1 * (x > 0)
    } else {
      x.t <- x
    }
    d <- vegan::vegdist(x.t, method = method)

    # Calculate sum of squares total
    SST <- SS(d)[2]

    # Size for labels
    nlev <- nlevels(as.factor(factEnv))

    # Split data by factors and calculate residuals
    lista <- split(as.data.frame(x.t), factEnv)
    dimxt <- dim(x.t)
    vdist <- lapply(lista, vegan::vegdist, method = method)
    SSi <- lapply(vdist, SS)
    SSi <- array(unlist(SSi), dim = c(1,2,nlev))

    # Calculate mean squares and pseudoF
    denR <- dimxt[1] - nlev
    denA <- nlev - 1

    SSR <- sum(SSi[,2,])
    SSA <- abs(SST - SSR)
    MSR <- (SSR/denR)
    MSA <- (SSA/denA)

    # Validate MS values
    if (MSR <= 0 || MSA <= 0) {
      stop("Mean squares (MSA or MSR) contain invalid values (e.g., non-positive).")
    }

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
  sel[sel[,1] <= -1, 1] <- 0
  sel[sel[,1] >= 2, 1] <- 1

  ones <- which(sel[,1] %in% 1)
  y0 <- H0Sim[ones,,resultsHa[i,1]]
  ya <- HaSim[ones,,resultsHa[i,1]]

  # Validate dimensions of H0 and Ha matrices
  if (!all(dim(y0) == dim(ya))) {
    stop("The dimensions of H0 and Ha data matrices do not match.")
  }

  yHa <- dim(y0)[2] - 2

  # Apply PERMANOVA to get pseudoF and mean squares
  result1 <- permanova_oneway(x = y0[, 1:yHa], factEnv = y0[,yHa+2],
                              transformation = transformation, method = method)
  result2 <- permanova_oneway(x = ya[, 1:yHa], factEnv = y0[,(yHa+2)],
                              transformation = transformation, method = method)

  # Create result matrix
  result0 <- matrix(nrow = 1, ncol = 4)
  colnames(result0) <- c("FobsH0", "FobsHa", "MSA", "MSR")

  # Gather the results and return
  result0[,1] <- result1$Fobs
  result0[,2] <- result2$Fobs
  result0[,3] <- result2$MSA
  result0[,4] <- result2$MSR
  return(result0)
}

#' PERMANOVA two-way
#'
#' Calculates observed F and mean squares for the residuals and among sites. This
#' function is a helper for [prep_data_nestedsymmetric()].
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
#' @importFrom vegan vegdist betadisper
#'
#' @seealso [vegan::vegdist()]
#'
#' @export
#' @keywords internal

## PERMANOVA Two factors ----
permanova_twoway <- function(x, factEnv, method = "bray",
                             transformation = "none", model = "nested.symmetric"){

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
    b = length(unique(as.factor(factEnv$sites)))  # number of sites (B)
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

    # fill the PERMANOVA table
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
    # Validate transformation method
    valid_transformations <- c("none", "square root", "fourth root", "Log (X+1)", "P/A")
    if (!transformation %in% valid_transformations) {
      stop("Invalid transformation method. Choose from: none, square root, fourth root, Log (X+1), P/A.")
    }

    # Validate distance method
    valid_methods <- c("euclidean", "bray", "jaccard", "gower")
    if (!method %in% valid_methods) {
      stop("Invalid distance method. Choose from: euclidean, bray, jaccard, gower.")
    }

    # Ensure dimensions match
    if (dim(x)[1] != dim(factEnv)[1]) {
      stop("The dimensions of the data matrix and the factor labels do not match.")
    }
    rm(valid_transformations, valid_methods)

    # Apply transformation and calculate distance matrix
    if (transformation == "square root") {
      x.t <- sqrt(x)
    } else if (transformation == "fourth root") {
      x.t <- sqrt(sqrt(x))
    } else if (transformation == "Log (X+1)") {
      x.t <- log(x + 1)
    } else if (transformation == "P/A") {
      x.t <- 1 * (x > 0)
    } else {
      x.t <- x
    }
    rm(x)
    d <- vegan::vegdist(x.t, method = method)

    # size for the labels we work with
    a = nlevels(as.factor(factEnv$sector))  # number of sectors (A)
    b = length(unique(as.factor(factEnv$sites)))   # number of sites (B)
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
    rm(secsite_groups, listR)

    # calculates SS_B(A)

    # Bray-Curtis distance matrix is calculated for each sector, then the centroids
    # are estimated using `betadisper()`. Negative values are taken out, and then
    # the SS for the Euclidian distance matrix is calculated. The SS are summed
    # to calculate SS_B(A)
    sector_groups <- split(rownames(factEnv), factEnv$sector)
    listBA <- sapply(sector_groups, function(rw) {
      dBA <- vegan::vegdist(x.t[rw,], method = "bray")
      tCentroid <- vegan::betadisper(dBA,
                                     group = factEnv[rw, 4],
                                     type = "centroid",
                                     bias.adjust = FALSE)
      Eig <- which(tCentroid$eig > 0)
      SS(vegan::vegdist(tCentroid$centroids[, Eig, drop = FALSE], method = "euclidean"))
    }, simplify = "array")
    SSBA <- sum(listBA[2,]) * nRep  # * nRep is added so that when calculating
    # the variation between nested groups it is weighted by the size of the subgroup
    rm(sector_groups, listBA)

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
  sel[sel[,1] <= -1, 1] <- 0
  sel[sel[,1] >= 2, 1] <- 1

  ones <- which(sel[,1] %in% 1)
  ones_n <- rep(ones, nSect)
  ones_s <- rep(c(0:(nSect-1)) * sites * N, each = length(ones))
  ones <- ones_n + ones_s

  y0 <- H0Sim[ones,,resultsHa[i,1]]
  rownames(y0) <- ones
  ya <- HaSim[ones,,resultsHa[i,1]]
  rownames(ya) <- ones
  factEnv_s <- factEnv[ones,]

  # Apply PERMANOVA to get F and mean squares
  result1 <- permanova_twoway(x = y0,
                              factEnv = factEnv_s,
                              transformation = transformation,
                              method = method,
                              model = model)
  result2 <- permanova_twoway(x = ya,
                              factEnv = factEnv_s,
                              transformation = transformation,
                              method = method,
                              model = model)
  # result0 <- matrix(nrow = 1, ncol = 5)
  result0 <- matrix(nrow = 1, ncol = 4)

  # Values of pseudoF for A are stored in the result, values for MSA and MSR
  # come from the dataset with Ha
  # colnames(result0) <- c("FobsH0", "FobsHa", "MSA", "MSBA", "MSR")
  colnames(result0) <- c("FobsH0", "FobsHa", "MSBA", "MSR")

  # Gather the results and return
  # result0[,1] <- result1["A", "F"]
  # result0[,2] <- result2["A", "F"]
  # result0[,3] <- result2["A", "MS"]
  # result0[,4] <- result2["B(A)", "MS"]
  # result0[,5] <- result2["R", "MS"]
  # return(result0)
  result0[,1] <- result1[1]
  result0[,2] <- result2[1]
  result0[,3] <- result2[3]
  result0[,4] <- result2[4]

  return(result0)

}

#' Balanced sampling 2_3
#'
#' Develops the experimental design based on the provided conditions
#'
#' @param i pointer to the index in the list of experimental designs to try.
#' @param NN Total number of iterations that the experiment will consider.
#' @param Y1 A data frame with two columns, one indicates the auxiliary variables
#' on which the sample must be balanced and the other contains the vector of integers
#' that defines the primary sampling units. This is used by
#' \code{sampling::balancedtwostage()}
#' @param mn A data frame with two columns, one indicates the number of primary
#' sampling units to be selected and the other the number of second-stage sampling
#' units to be selected in the iteration. This is used by
#' \code{sampling::balancedtwostage()}
#' @param nSect Total number of sectors to be simulated in each data set.
#' @param M Total number of replicates to be simulated in each data set.
#' @param N Total number of samples to be simulated in each site.
#' @param H0Sim simulated community from \code{SSP::simdata()} in which H0 is true.
#' @param HaSim simulated community from \code{SSP::simdata()} in which H0 is false.
#' @param resultsHa helper matrix that stores labels and later the results.
#' @param factEnv a data frame for indicating the treatment, replicate and sampling
#' unit lables in each experiment.
#' @param transformation Mathematical function to reduce the weight of very
#' dominant species.
#' @param method appropriate distance/dissimilarity metric (e.g. Gower,
#' Bray–Curtis, Jaccard, etc).
#' @param model which algorithm to use for the calculation? At the moment, the only
#' option is "nested.symmetric".
#'
#' @return a data frame with values for observed F (for H0 and Ha), and the Ha mean
#' squares for residuals (MS_R) and variation among sites (MS_B(A)).
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

balanced_sampling2_3 <- function(i, NN, Y1, mn, nSect, M, N, H0Sim, HaSim, resultsHa,
                                 factEnv, transformation, method, model){
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

  result_0 <- permanova_twoway(x = y0,
                               factEnv = factEnvX,
                               transformation=transformation,
                               method = method,
                               model = model)
  result_a <- permanova_twoway(x = ya,
                               factEnv = factEnvX,
                               transformation=transformation,
                               method = method,
                               model = model)

  if(length(result_0) != 4){result_0 <- c(NA, NA, NA, NA)}
  if(length(result_a) != 4){result_a <- c(NA, NA, NA, NA)}
  # Assemble results
  result1 <- c(result_0[1], result_a[1],
               result_a[3], result_a[4])


  return(result1)
}


