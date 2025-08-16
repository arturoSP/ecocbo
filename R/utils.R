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
#'
#' @keywords internal
#' @noRd
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
#'
#' @keywords internal
#' @noRd
#'

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
#'
#' @keywords internal
#' @noRd
#'

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
#' @param model Character. Select the model to use.
#' @param sigma2_est An estimation on the value for sigma_e^2.
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
#'
#' @keywords internal
#' @noRd
#'

## PERMANOVA Two factors ----
permanova_twoway <- function(x, factEnv, method = "bray",
                             transformation = "none", model = "nested.symmetric",
                             sigma2_est = NULL){

  pseudoF_2Orthogonal <- function(x, factEnv, method = "bray", transformation = "none",
                                  sigma2_est){
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

  pseudoF_2NestedSymmetric <- function(x, factEnv, method, transformation,
                                       sigma2_est){
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
    MSA <- SSA / DoFA
    # MSBA <- SSBA / DoFBA
    MSBAsim <- (sigma2_est / DoFBA) * rchisq(1, df = DoFBA) # sigma2_est <- scompvar(permanova(datospiloto))$BA
    MSBA <- (SSBA / DoFBA) + MSBAsim
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
    Results <- pseudoF_2NestedSymmetric(x, factEnv, method, transformation,
                                        sigma2_est)
  } else if(model == "nested.asymmetric"){
    Results <- pseudoF_2NestedSymmetric(x, factEnv, method, transformation,
                                        sigma2_est)
  } else {
    Results <- pseudoF_2Orthogonal(x, factEnv, method, transformation,
                                   sigma2_est)
  }

  Fobs <- Results
  return(Fobs)
}

#' Balanced sampling 2
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
#'
#' @keywords internal
#' @noRd
#'

balanced_sampling2 <- function(i, NN, Y1, mn, nSect, M, N, H0Sim, HaSim,
                               resultsHa, factEnv, transformation, method,
                               model, sigma2_est){
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
                               model = model,
                               sigma2_est = sigma2_est)
  result_a <- permanova_twoway(x = ya,
                               factEnv = factEnvX,
                               transformation=transformation,
                               method = method,
                               model = model,
                               sigma2_est = sigma2_est)

  if(length(result_0) != 4){result_0 <- c(NA, NA, NA, NA)}
  if(length(result_a) != 4){result_a <- c(NA, NA, NA, NA)}
  # Assemble results
  result1 <- c(result_0[1], result_a[1],
               result_a[3], result_a[4])


  return(result1)
}

#' mimimun permutations
#'
#' Determine if sampling effort allows for at least 100 permutations
#'
#' @param model which algorithm to use for the calculation? At the moment, the only
#' option is "nested.symmetric".
#' @param a Integer. Levels for the treatment factor.
#' @param m Integer. Levels for site within treatment. Only used in Nested Symmetric
#' experiments.
#' @param n Integer. Replicates in the experiment (either per treatment or site).
#'
#' @return Logical. TRUE if at least 100 permutations are guaranteed.
#'
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @keywords internal
#' @noRd
#'

minimum_cbo <- function(model, a, n, m = NULL){

  thr <- c(log(100))

  if(model == "single.factor"){
    permA <- lgamma(a * n+1) - lgamma(a + 1) - a * lgamma(n + 1)

    return(permA > thr)
  } else {
    if(is.null(m)) stop("m is required for the nested model")

    permA <- lgamma(a * m + 1) - lgamma(a + 1) - a * lgamma(m + 1)
    permBA <- a * (lgamma(m * n + 1) - lgamma(m + 1) - m * lgamma(n + 1))

    return((permA > thr) & (permBA > thr))
  }
}

#' PERMANOVA by exchanging labels
#'
#' Classic Permutational Multivariate Analysis of Variance
#'
#' @param data Data frame where columns represent species names and rows correspond
#' to samples.
#' @param factEnv Data frame containing the environmental factors for the data.
#' @param method Character. Dissimilarity metric used [vegan::vegdist()]. Common
#' options include: "Gower", "Bray–Curtis", "Jaccard", etc.
#' @param transformation Character. Transformation applied to reduce the weight
#' of dominant species: "square root", "fourth root", "Log (X+1)", "P/A", "none".
#' @param dummy Character. Transformation applied to reduce the weight
#' of dominant species: "square root", "fourth root", "Log (X+1)", "P/A", "none".
#' @param model Character. Select the model to use. Options, so far, are
#' @param k Integer. Number of resampling iterations. Defaults to 50.
#'
#' @return PERMANOVA table
#'
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @keywords internal
#' @noRd
#'

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

  # Only necessary for full PERMANOVA
  # permutaciones_rep <- replicate(999,
  #                                  expr = factEnvP[sample(nrow(factEnvP)), ],
  #                                  simplify= FALSE)

  # calculates SS for all
  SST <- SS(d)[2]

  # permList <- vector("list", 999 + 1)
  # degrees of freedom
  DoFA <- a - 1
  DoFBA <- a * (b - 1)
  DoFR <- a * b * (nRep - 1)
  DoFT <- (a * b * nRep) - 1

  # Only necessary for full PERMANOVA

  # for(i in c(1:999)){
  #   # dataframe en esta iteración
  #   currentPerm <- permutaciones_rep[[i]]
  #
  #   # calculates SS within replicates
  #   # secsite_groups <- split(rownames(permutaciones_rep[[i]]), factEnvP$secsit)
  #   secsite_groups <- split(rownames(x.t), currentPerm$secsit)
  #   listR <- sapply(secsite_groups, function(rw) {
  #     SS(vegan::vegdist(x.t[rw,], method = method))
  #   }, simplify = "array")
  #   SSR <- sum(listR[2,])
  #   # calculates SS_B(A)
  #   # sector_groups <- split(rownames(permutaciones_rep[[i]]), factEnvP[,1])
  #   sector_groups <- split(rownames(x.t), currentPerm$Locality)
  #
  #   listBA <- sapply(sector_groups, function(rw) {
  #     dBA <- vegan::vegdist(x.t[rw,], method = "bray")
  #     tCentroid <- vegan::betadisper(dBA,
  #                                    group = currentPerm[rw, "secsit"],
  #                                    type = "centroid",
  #                                    bias.adjust = FALSE)
  #     Eig <- which(tCentroid$eig > 0)
  #     SS(vegan::vegdist(tCentroid$centroids[, Eig, drop = FALSE],
  #                       method = "euclidean"))
  #   }, simplify = "array")
  #   SSBA <- sum(listBA[2,]) * nRep  # * nRep is added so that when calculating
  #   # the variation between nested groups it is weighted by the size of the subgroup
  #   # rm(sector_groups, listBA)
  #
  #   # calculates SSA
  #   SSA <- SST - SSBA - SSR
  #
  #   # fill the permanova table
  #   # mean squares
  #   MSA <- SSA / DoFA
  #   MSBA <- SSBA / DoFBA
  #   # MSBAsim <- (sigma2_est / DoFBA) * rchisq(1, df = DoFBA) # sigma2_est <- scompvar(permanova(datospiloto))$BA
  #   MSBA <- (SSBA / DoFBA) #+ MSBAsim
  #   MSR <- SSR / DoFR
  #
  #   # observed pseudoF
  #   FobsA <- MSA / MSBA
  #   FobsBA <- MSBA / MSR
  #
  #   Fobs <- as.data.frame(matrix(nrow = 4, ncol = 4))
  #   colnames(Fobs) <- c("SS", "DoF", "MS", "F")
  #   rownames(Fobs) <- c("A", "B(A)", "R", "T")
  #   Fobs[1,1] <- SSA
  #   Fobs[1,2] <- DoFA
  #   Fobs[1,3] <- MSA
  #   Fobs[1,4] <- FobsA
  #   Fobs[2,1] <- SSBA
  #   Fobs[2,2] <- DoFBA
  #   Fobs[2,3] <- MSBA
  #   Fobs[2,4] <- FobsBA
  #   Fobs[3,1] <- SSR
  #   Fobs[3,2] <- DoFR
  #   Fobs[3,3] <- MSR
  #   Fobs[4,1] <- SST
  #   Fobs[4,2] <- DoFT
  #
  #   permList[[i]] <- Fobs
  # }

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

  # Only necessary for full PERMANOVA
  # permList[[1000]] <- Fobs

  # If this were full PERMANOVA, it's necessary to change the return to permList
  return(Fobs)
}
