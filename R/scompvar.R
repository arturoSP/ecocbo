#' Simulated components of variation
#'
#' @param data List that results from beta().
#' @param m Site label to be used as basis for the computation.
#' @param n Number of samples to be considered.
#'
#' @return A data frame containing the values for the variation component among sites and in the residuals.
#' @export
#'
#' @examples
#' scompvar(data = epiBetaR)
#' scompvar(data = epiBetaR, n = 5, m = 2)

scompvar <- function(data, n = NULL, m = NULL){
  # Determine variation components  ----
  # as these values are necessary for the cost-benefit optimization model.

  if(class(data) != "ecocbo_beta")
    stop("data is not the right class(\"ecocbo_beta\")")

  # read the results matrix to use Ha mean squares
  resultsBeta <- data$Results

  ## Validating data ----
  if(is.null(n)){n <- max(resultsBeta$n)}
  if(ceiling(n) != floor(n)){stop("n must be integer")}
  if(n <= 1){stop("n must be larger than 1")}
  if(n > max(resultsBeta$n)){stop("n is larger than the simulated n value")}

  if(is.null(m)){m <- max(resultsBeta$m)}
  if(ceiling(m) != floor(m)){stop("m must be integer")}
  if(m <= 1){stop("m must be larger than 1")}
  if(m > max(resultsBeta$m)){stop("m is larger than the simulated m value")}

  resultsBeta <- resultsBeta[resultsBeta$m == m & resultsBeta$n == n,]

  # Creates an empty matrix to store the results
  compVar <- data.frame(compVarA = NA, compVarR = NA)

  # Computes the mean for variation components (as per Table 9.3)
  # σe = RMS
  compVar[,2] <- mean(resultsBeta[,8], na.rm = T)
  # σB(A) = (AMS - σe) / n
  MSA <- mean(resultsBeta[,7], na.rm = T)
  compVar[,1] <- (MSA - compVar[,2]) / n

  return(compVar)
}
