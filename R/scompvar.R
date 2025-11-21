#' Simulated Components of Variation
#'
#' Computes the average components of variation among sampling units and within
#' samples in relation to sampling effort.
#'
#' @param data Object of class `"ecocbo_data"` obtained from [prep_data()].
#' @param m Optional. Integer. Number of replicates to consider.
#' @param n Optional. Integer. Number of samples to consider.
#'
#' @details
#' If \code{m} or \code{n} are set to `NULL`, the function automatically uses the
#' largest available values from the experimental design set in [sim_beta()].
#'
#' @return A data frame containing the values for the variation component
#' among sites \code{compVarA} and in the residuals \code{compVarR}.
#'
#' @author Edlin Guerra-Castro (\email{edlinguerra@@gmail.com}), Arturo Sanchez-Porras
#'
#' @references
#' - Underwood, A. J. (1997). Experiments in ecology: their logical
#' design and interpretation using analysis of variance. Cambridge university
#' press.
#' - Underwood, A. J., & Chapman, M. G. (2003). Power, precaution,
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
#' @importFrom dplyr group_by summarise
#'
#' @examples
#' scompvar(data = simResults)
#' scompvar(data = simResults, n = 5, m = 2)
#'

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

  if(data$model != "single.factor"){
    if(is.null(m)){m <- max(resultsBeta$m)}
    if(ceiling(m) != floor(m)){stop("m must be integer")}
    if(m <= 1){stop("m must be larger than 1")}
    if(m > max(resultsBeta$m)){stop("m is larger than the simulated m value")}
  }

  if(data$model == "single.factor"){

    resultsBeta <- resultsBeta[resultsBeta$n == n,]

    # Creates an empty matrix to store the results
    compVar <- data.frame(Source = c("Residual"),
                          Est.var.comp = NA)

    # Computes the mean for variation components as per Table 8.5 (Quinn & Keough, 2002) ----
    # σe = MSR
    # This overall average is being ditched for a sectorized average
    # compVar[1,2] <- mean(resultsBeta[,7], na.rm = T)
    compVar[1,2] <- resultsBeta |>
      dplyr::group_by(dat.sim) |>
      dplyr::summarise(MSR = mean(MSR)) |>
      dplyr::summarise(max(MSR))

    # We skip this calculation given that it isn't relevant for subsequent steps
    # σA = (MSA - MSR) / n
    # compVar[1,2] <- (mean(resultsBeta[,7], na.rm = T) - compVar[2,2]) / n
  } else {

    resultsBeta <- resultsBeta[resultsBeta$m == m & resultsBeta$n == n,]

    # Creates an empty matrix to store the results
    compVar <- data.frame(Source = c("B(A)", "Residual"),
                          Est.var.comp = NA)

    # Computes the mean for variation components as per Table 9.5 (Quinn & Keough, 2002) ----
    # σe = MSR
    # This overall average is being ditched for a sectorized average
    # compVar[2,2] <- mean(resultsBeta[,8], na.rm = T)
    compVar[2,2] <- resultsBeta |>
      dplyr::group_by(dat.sim) |>
      dplyr::summarise(MSR = mean(MSR)) |>
      dplyr::summarise(max(MSR))

    # σB(A) = (MSBA - MSR) / n
    # This overall average is being ditched for a sectorized average
    # compVar[1,2] <- (mean(resultsBeta[,7], na.rm = T) - compVar[2,2]) / n
    compVar[1,2] <- resultsBeta |>
      dplyr::group_by(dat.sim) |>
      dplyr::summarise(MSBA = mean(`MSB(A)`)) |>
      dplyr::summarise(max(MSBA))

    compVar[1,2] <- (compVar[1,2] - compVar[2,2]) / n

    # We skip this calculation given that it isn't relevant for subsequent steps
    # σA = (MSA - MSBA) / (n * m)
    # compVar[1,2] <- (mean(resultsBeta[,7], na.rm = T) - compVar[2,2]) / (n * m)
  }

  return(compVar)
}
