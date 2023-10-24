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

permanova_oneway <- function(x, factEnv, type = "P", method = "bray", transformation = "none"){
  ## PERMANOVA ----
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
    denA <- nlev - 1
    denR <- dimxt[1] - nlev

    # Results
    SSR <- sum(SSi[,2,])
    SSA <- abs(SST - SSR)
    MSA <- (SSA/denA)
    MSR <- (SSR/denR)
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

permanova_twoway <- function(x, factEnv, method = "bray", transformation = "none", type = "P"){
  # PERMANOVA two way nested ----
  pseudoF_P2 <- function(x, factEnv, method = "bray", transformation = "none"){
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
    a = nlevels(as.factor(factEnv$sector))
    b = nlevels(as.factor(factEnv$sites))
    # n = max(factEnv$N)
    factEnv["secsit"] <- paste0(factEnv$sector, factEnv$site)
    n = nlevels(as.factor(factEnv$secsit))
    N = n * a * b

    # sum of squares total, Huygen for all
    SST <- SS(d)
    SST <- SST[2]
    #SST <- (SST[2] * SST[1]) / N

    # calculate SS for A
    # A1 <- rownames(factEnv[factEnv$sector == "E",])
    # A2 <- rownames(factEnv[factEnv$sector == "I",])
    #listA <- lapply(list(A1 = vegdist(x.t[A1,]), A2 = vegdist(x.t[A2,])), SS)
    # SSdA <- sum(listA[,2,]) / (b * n)
    # SSdA <- sum(listA[,2,]) / b

    rwNmA <- as.data.frame(matrix(nrow = (nrow(factEnv)/a), ncol = a))
    colnames(rwNmA) <- levels(as.factor(factEnv$sector))

    listA <- list()
    for(i in levels(as.factor(factEnv$sector))){
      RwNm <- rownames(factEnv[factEnv$sector == i,])
      listA[[i]] <- SS(vegdist(x.t[RwNm,]))
    }

    listA <- array(unlist(listA), dim = c(1,2,a))
    SSdA <- sum(listA[,2,])
    SSA <- SST - SSdA

    # calculate SS for B
    # B1 <- rownames(factEnv[factEnv$sites == 1,])
    # B2 <- rownames(factEnv[factEnv$sites == 2,])
    # listB <- lapply(list(B1 = vegdist(x.t[B1,]), B2 = vegdist(x.t[B2,])), SS)
    # listB <- array(unlist(listB), dim = c(1,2,2))
    # SSdB <- sum(listB[,2,]) / (a * n)
    # SSdB <- sum(listB[,2,]) / a

    rwNmB <- as.data.frame(matrix(nrow = (nrow(factEnv)/b), ncol = b))
    colnames(rwNmB) <- levels(as.factor(factEnv$sites))

    listB <- list()
    for(i in levels(as.factor(factEnv$sites))){
      RwNm <- rownames(factEnv[factEnv$sites == i,])
      listB[[i]] <- SS(vegdist(x.t[RwNm,]))
    }

    listB <- array(unlist(listB), dim = c(1,2,b))
    SSdB <- sum(listB[,2,])
    SSB <- SST - SSdB

    # calculate SS for residuals
    # RE1 <- rownames(factEnv[factEnv$secsit == "E1",])
    # RE2 <- rownames(factEnv[factEnv$secsit == "E2",])
    # RI1 <- rownames(factEnv[factEnv$secsit == "I1",])
    # RI2 <- rownames(factEnv[factEnv$secsit == "I2",])
    # listR <- lapply(list(vegdist(x.t[RE1,], method = method), vegdist(x.t[RE2,], method = method),
    #                      vegdist(x.t[RI1,], method = method), vegdist(x.t[RI2,], method = method)),
    #                 SS)
    # listR <- array(unlist(listR), dim = c(1,2,(a * b)))
    # SSR <- sum(listR[,2,]) / n

    r <- length(unique(factEnv$secsit))
    rwNmR <- as.data.frame(matrix(nrow = (nrow(factEnv) / r), ncol = r))
    colnames(rwNmR) <- levels(as.factor(factEnv$secsit))

    listR <- list()
    for(i in levels(as.factor(factEnv$secsit))){
      RwNm <- rownames(factEnv[factEnv$secsit == i,])
      listR[[i]] <- SS(vegdist(x.t[RwNm,]))
    }

    listR <- array(unlist(listR), dim = c(1,2,r))
    SSR <- sum(listR[,2,])

    # calculate SS for interaction A-B
    SSAB <- SST - SSA - SSB - SSR

    # SSA + SSB + SSAB + SSR

    # fill the permanova table
    DoFA <- a - 1
    DoFB <- b - 1
    DoFAB <- DoFA * DoFB
    DoFR <- a * b * (n - 1)
    DoFT <- (a * b * n) - 1

    MSA <- SSA / DoFA
    MSB <- SSB / DoFB
    MSAB <- SSAB / DoFR
    MSR <- SSR / DoFT

    FobsA <- MSAB / MSA
    FobsB <- MSAB / MSB
    # FobsA <- MSA / MSAB
    # FobsB <- MSB / MSAB

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
    Fobs[4,1] <- SSR
    Fobs[4,2] <- DoFR
    Fobs[4,3] <- MSR
    Fobs[5,1] <- SST
    Fobs[5,2] <- DoFT

    return(Fobs)
  }

  ## Main function which returns pseudoF and MS ----
  if(dim(x)[1] != dim(factEnv)[1])
    stop("x and factEnv dimensions do not match")
  if(type != "P" & type != "BF")
    stop("Possible values for type are P and BF")

  Results <- pseudoF_P2(x, factEnv, method, transformation)
  Fobs <- Results
  return(Fobs)
}



