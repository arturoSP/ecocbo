# Helper functions ----
## Sum of Squares using Huygen theorem ----
SS <- function (d) {
  ss <- numeric(2)
  ss[1] <- dim(as.matrix(d))[1]
  ss[2] <- sum(d^2)/ss[1]
  return(ss)
}

## PERMANOVA function ----
pseudoF_P <- function(x, factor, method = "bray"){
  d <- vegan::vegdist(x, method = method)
  TSS <- SS(d)[2]

  # size for labels
  nlev <- nlevels(as.factor(as.matrix(factor)))

  # Calculate the SS for residuals
  lista <- split(x, factor)
  vdist <- lapply(lista, vegdist)
  SSi <- lapply(vdist, SS)
  SSi <- array(unlist(SSi), dim = c(1,2,nlev))

  # Calculate denominators
  denA <- nlev - 1
  denR <- dim(x)[1] - nlev

  # Results
  RSS <- sum(SSi[,2,])
  ASS <- TSS - RSS
  AMS <- (ASS/denA)
  RMS <- (RSS/denR)
  Fobs <- AMS/RMS
  Fobs <- round(data.frame(ASS, RSS, TSS, denA, denR, AMS, RMS, Fobs), 4)
  return(Fobs)
}

# Main function which returns pseudoF and MS ----
permanova_reduced <- function(x, factor, type = "P", method = "bray"){
  Results <- if(type == "P"){
    pseudoF_P(x, factor)
  } else {
    pseudoF_BF(x, factor)
  }
  Fobs <- Results[1,c(6,7,8)]
  return(Fobs)
}
