#' Cost-Benefit Optimization for Sampling Effort
#'
#' Given a table of statistical power estimates produced by \code{\link{sim_beta}},
#' \code{sim_cbo} finds the sampling design (number of replicates/site and sites)
#' that minimizes total cost while achieving a user‐specified power threshold.
#'
#' @param data Object of class \code{"ecocbo_beta"}, as returned by
#' \code{\link{sim_beta}}.
#' @param cm Numeric. Fixed cost per replicate.
#' @param cn Numeric. Cost per sampling unit.
#' @param perm Integer. Minimum number of permutations needed to reject the null
#' hypothesis. Defaults to 100, as it would allow for rejecting with alpha = 0.05,
#' the user can change this value to make the testing more strict (e.g. 200 for
#' testing alpha = 0.01 or 5000 for testing alpha = 0.001).
#'
#' @return A data frame with one row per candidate design. In the single factor
#' case, the results include the available \code{n} values, their statistical
#' power and cost. For the nested symmetric experiments, the results include all
#' the available values for \code{m}, the optimal \code{n}, according to the
#' power, and the associated cost. The results also mark a suggested sampling
#' effort, based on the cost and power range as selected by the user.
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
#' [scompvar()]
#' [Underwood_cbo()]
#'
#' @aliases simcbo
#'
#' @export
#' @importFrom dplyr filter group_by arrange slice mutate
#'
#' @examples
#' # Optimization of single factor experiment
#' sim_cbo(data = epiBetaR, cn = 80)
#'
#' # Optimization of a nested factor experiment
#' sim_cbo(data = betaNested, cn = 80, cm = 180)
#'

sim_cbo <- function(data, cn, cm = NULL, perm = 100){
  # Obtaining parameters from the ecocbo_beta object
  powr <- subset(data$Power, select = -c(Beta, fCrit))
  objective <- 1 - data$alpha
  model <- data$model
  a <- data$a

  if(model == "nested.symmetric"){
    if(is.null(cm)){
      stop("Cost for sites is missing.")
    }

    # Calculates total cost
    powr$Cost <- powr$m * cm + powr$m * powr$n * cn

    # Find the combinations that achieve the criteria
    powr$OptPower <- powr$Power >= objective
    powr$OptPerm <- minimum_cbo(model, a = a,
                                m = powr$m, n = powr$n,
                                perm)

    # Filter out the cases when n = 2
    powr[powr$n == 2, 5] <- FALSE

    # Find the cases when power and permutations are valid
    ideal <- powr |>
      dplyr::filter(OptPower, OptPerm)

    if(nrow(ideal) == 0){
      warning("No power values above the specified precision.")
      ideal <- powr |>
        group_by(m) |>
        arrange(desc(Power)) |>
        slice(1) |>
        mutate(OptPower = TRUE) |>
        filter(OptPerm)

      powr <- merge(powr[,-c(5)], ideal[,c(1,2,5)],
            all.x=TRUE) |>
        mutate(OptPower = ifelse(is.na(OptPower), FALSE, TRUE))
    }

    idCost <- which.min(ideal$Cost)
    idBest <- ideal[idCost,]

    powr$OptCost <- FALSE
    powr[powr$m == idBest$m & powr$n == idBest$n, "OptCost"] <- TRUE

  } else {
    # Calculating the cost
    powr$Cost <- powr$n * cn

    # Find the iterations that meet the desired power range
    powr$OptPower <- powr$Power >= objective
    powr$OptPerm <- minimum_cbo(model, a = a,
                                n = powr$n, perm = perm)

    # Filter out n = 2
    powr[powr$n == 2, 5] <- FALSE

    # Find the cases when power and permutations are valid
    ideal <- powr |>
      dplyr::filter(OptPower, OptPerm)

    if(nrow(ideal) == 0){
      warning("No power values above the specified precision.")
    }

    idCost <- which.min(ideal$Cost)
    idBest <- ideal[idCost,]

    powr$OptCost <- FALSE
    powr[powr$n == idBest$n, "OptCost"] <- TRUE

  }

  class(powr) <- c("cbo_result", class(powr))
  return(powr)
}


#-------------------------------------------
## S3Methods print()
#-------------------------------------------

#' S3Methods for Printing
#'
#' @name print.cbo_result
#'
#' @method print cbo_result
#'
#' @usage
#' \method{print}{cbo_result}(x, ...)
#'
#' @description prints for \code{ecocbo::sim_cbo()} objects.
#'
#' @param x Object from \code{ecocbo::sim_cbo()} function.
#'
#' @param ... Additional arguments
#'
#' @return Prints a summary for the results of \code{ecocbo::sim_cbo()} function,
#' showing in an ordered matrix the suggested experimental design, according to
#' cost and estimated power.
#'
#' @importFrom dplyr mutate group_by ungroup arrange slice filter select desc
#'
#' @export
#' @keywords internal

print.cbo_result <- function(x, ...){
  # Create the subset of suggested options
  if(length(x) == 6){
    x1 <- x |>
      mutate(Suggested = ifelse(OptCost, "***", "")) |>
      select(n, Power, Cost, Suggested)
  } else {
    x1 <- x |>
      filter(OptPower) |>
      group_by(m) |>
      arrange(n, Cost) |>
      slice(1) |>
      ungroup() |>
      mutate(Suggested = ifelse(OptCost, "***", "")) |>
      select(m, n, Power, Cost, Suggested)
  }

  # Print
  cat("Sampling designs that meet the required power:\n")
  print.data.frame(x1, row.names=FALSE, right=TRUE)
  cat("\nThe listed cost is per treatment.")
  invisible(x1)

}
