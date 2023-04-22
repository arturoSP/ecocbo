#' Calculate pseudoF by using PERMANOVA
#'
#' @param x Ecological community to be used.
#' @param factor List of labels for the sites.
#' @param type PERMANOVA algorithm to be used. Options are "P" and "BF". "P" will also throw out the values for MeanSquares, as they are used for other functions in the package.
#' @param method Dissimilaty index to be used to calculate dissimilarity matrix.
#'
#' @return A data frame with information regarding pseudoF and, if required, mean square values.
#' @export
#'
#' @examples
#' perHa1 <- epiDat[,-1]
#' envHa1 <- epiDat[,1]
#'
#' permanova_reduced(x = perHa1, factor = envHa1, type = "P", method = "bray")
#' permanova_reduced(x = perHa1, factor = envHa1, type = "BF", method = "bray")

permanova_reduced <- function(x, factor, type = "P", method = "bray"){
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
  vdist <- lapply(lista, vegan::vegdist)
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

pseudoF_BF <- function(x, factor, method = "bray"){
  d <- vegan::vegdist(x, method = method)
  TSS <- SS(d)[2]

  # Tamaños de etiqueta
  lev <- table(factor)
  nlev <- length(lev)

  # Vectores de apoyo vacíos
  Var <- numeric(nlev)
  d.res <- Var

  # Cálculo de suma de cuadrados y varianzas para residuales
  lista <- split(x, factor)
  vdist <- lapply(lista, vegan::vegdist)
  SSi <- lapply(vdist, SS)
  SSi <- array(unlist(SSi), dim = c(1,2,nlev))

  Var <- SSi[,2,]/(SSi[,1,]-1)
  d.res <- (1-(lev/sum(lev)))*Var
  den <- sum(d.res)

  # Resultados
  RSS <- sum(SSi[,2,])
  ASS<-TSS-RSS
  Fobs<- ASS/den

  Fobs <- round(data.frame(ASS, RSS, TSS, den, Fobs), 4)
  return(Fobs)
}

# Main function which returns pseudoF and MS ----
#permanova_reduced <- function(x, factor, type = "P", method = "bray"){
  Results <- if(type == "P"){
    pseudoF_P(x, factor)[1,c(6,7,8)]
  } else {
    pseudoF_BF(x, factor)[1,5]
  }
  Fobs <- Results
  return(Fobs)
}
