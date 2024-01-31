#' Simulated components of variation
#'
#' \code{scompvar} can be used to calculate the average component of variation
#' among units and the average component of variation within samples in terms
#' of sampling effort.
#'
#' @param data Object of class "ecocbo_data" that results from [prep_data()].
#' @param m Site label to be used as basis for the computation. Defaults to NULL.
#' @param n Number of samples to be considered. Defaults to NULL.
#'
#' @note
#' If \code{m} or \code{n} are left as NULL, the function will calculate
#' the components of variation using the largest available values as set in
#' the experimental design in [sim_beta()].
#'
#' @return A data frame containing the values for the variation component
#' among sites \code{compVarA} and in the residuals \code{compVarR}.
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
#' [prep_data()]
#'
#' @aliases compvar
#'
#' @export
#'
#' @examples
#' scompvar(data = simResults)
#' scompvar(data = simResults, n = 5, m = 2)

scompvar <- function(data, n = NULL, m = NULL){
  # Determine variation components  ----
  # as these values are necessary for the cost-benefit optimization model.

  if(!inherits(data, "ecocbo_data"))
    stop("Data is not the right class(\"ecocbo_data\")")

  # read the results matrix to use Ha mean squares
  resultsBeta <- as.data.frame(data$Results)

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

  if(data$model == "single.factor"){
    # Creates an empty matrix to store the results
    compVar <- data.frame(Source = c("A", "Residual"),
                          Est.var.comp = NA)

    # Computes the mean for variation components as per Table 8.5 (Quinn & Keough, 2002) ----
    # σe = MSR
    compVar[2,2] <- mean(resultsBeta[,8], na.rm = T)
    # σA = (MSA - MSR) / n
    compVar[1,2] <- (mean(resultsBeta[,7], na.rm = T) - compVar[2,2]) / n
  } else {
    # Creates an empty matrix to store the results
    compVar <- data.frame(Source = c("A", "B(A)", "Residual"),
                          Est.var.comp = NA)

    # Computes the mean for variation components as per Table 9.5 (Quinn & Keough, 2002) ----
    # σe = MSR
    compVar[3,2] <- mean(resultsBeta[,9], na.rm = T)
    # σB(A) = (MSBA - MSR) / n
    compVar[2,2] <- (mean(resultsBeta[,8], na.rm = T) - compVar[3,2]) / n
    # σA = (MSA - MSBA) / (n * m)
    compVar[1,2] <- (mean(resultsBeta[,7], na.rm = T) - compVar[2,2]) / (n * m)
  }

  return(compVar)
}
