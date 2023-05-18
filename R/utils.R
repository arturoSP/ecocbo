permanova_oneway <- function(x, factEnv, type = "P", method = "bray", transformation = "none"){

  # Helper functions ----
  ## Sum of Squares using Huygen theorem ----
  SS <- function (d) {
    ss <- numeric(2)
    ss[1] <- dim(as.matrix(d))[1]
    ss[2] <- sum(d^2)/ss[1]
    return(ss)
  }

  ## PERMANOVA functions ----
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

  # Prueba Brown-Forsythe ---
  # está desactivada, pero guardada por si se quiere usar algún día
  # pseudoF_BF <- function(x, factor, method = "bray"){
  #   d <- vegan::vegdist(x, method = method)
  #   TSS <- SS(d)[2]
  #
  #   # Size for labels
  #   lev <- table(factor)
  #   nlev <- length(lev)
  #
  #   # Empty helper vectors
  #   Var <- numeric(nlev)
  #   d.res <- Var
  #
  #   # Calculate SS and Var for residuals
  #   lista <- split(x, factor)
  #   vdist <- lapply(lista, vegan::vegdist)
  #   SSi <- lapply(vdist, SS)
  #   SSi <- array(unlist(SSi), dim = c(1,2,nlev))
  #
  #   Var <- SSi[,2,]/(SSi[,1,]-1)
  #   d.res <- (1-(lev/sum(lev)))*Var
  #   den <- sum(d.res)
  #
  #   # Resultados
  #   RSS <- sum(SSi[,2,])
  #   ASS<- abs(TSS-RSS)
  #   Fobs<- ASS/den
  #
  #   Fobs <- data.frame(ASS, RSS, TSS, den, Fobs)
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

balanced_sampling <- function(i, Y, mm, nn, YPU, H0Sim, HaSim, resultsHa, transformation, method){
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
