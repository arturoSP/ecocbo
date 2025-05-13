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
#'
#' @examples
#' # Optimization of single factor experiment
#' sim_cbo(data = epiBetaR, cn = 80)
#'
#' # Optimization of a nested factor experiment
#' sim_cbo(data = betaNested, cn = 80, cm = 180)
#'

sim_cbo <- function(data, cn, cm = NULL){
  # Obtaining parameters from the ecocbo_beta object
  # data <- beta2
  # beta1 <- sim_beta(data = simResults, alpha = 0.05)
  powr <- data$Power
  objective <- 1 - data$alpha
  model <- data$model

  if(model == "nested.symmetric"){
    # Calculates total cost
    powr$Cost <- powr$m * cm + powr$m * powr$n * cn

    # Find the combinations that achieve the criteria
    powr$OptPower <- powr$Power >= objective

    costOK <- powr$Cost[powr$OptPower]

    if(length(costOK) == 0){
      warning("No power values above the specified precision.")
      for(i in unique(powr[,1])){
        powr[powr[,1] ==i & powr$Power == max(powr[powr[,1] == i,3]), 7] <- TRUE
        ideal <- powr[powr[,7] == TRUE, -c(4,5)]
      }
    } else {
      ideal <- powr[powr[,3] >= objective, -c(4,5)]
    }

    # Obtaining the unique identifiers for n or m
    emes <- unique(ideal[,1])

    # Selects the minimum sampling effort, considering the power requirements
    # then it calculates the cost for such an experiment.
    ideally <- sapply(emes, function(rw){
      idealTrimm <- ideal[ideal[,1] == rw,]
      # Selecting the minimun n for each m
      enes <- min(idealTrimm[,2])
      idealTrimm <- idealTrimm[idealTrimm[,2] == enes,]
      # Calculating the associated cost
      idealTrimm[,"Cost"] <- idealTrimm[,1] * cm +
        idealTrimm[,2] * cn
      return(idealTrimm)
    }, simplify=F)

    # Puts together all the iterations
    powerCost <- do.call(rbind, ideally)
    rownames(powerCost) <- NULL

    # Find the minimum number of samples within the desired range
    powerCost$OptCost <- powerCost$Cost == min(powerCost$Cost[powerCost$OptPower == TRUE])

    # Completes the dataframe in case any m didn't reach the range
    powerCost <- merge(x = data.frame(m = unique(powr[,1])),
                       y = powerCost,
                       by = "m",
                       all.x = TRUE,
                       sort = TRUE)

    # Fills the empty data with the highest attainable power, then calculates
    # the associated cost
    for(i in powerCost[,1]){
      if(is.na(powerCost[powerCost[,1] == i,2])){
        # Finds the maximum n
        powerCost[powerCost[,1] == i,2] <- max(powr[powr[,1] == i,2])
        # Finds the maximum power
        powerCost[powerCost[,1] == i,3] <- powr$Power[powr$m == i &
                                                        powr$n == powerCost[powerCost[,1] == i,2]]
        # Calculates the associated cost
        powerCost[powerCost[,1] == i,4] <- powerCost[powerCost[,1] == i,1] * cm +
          powerCost[powerCost[,1] == i,2] * cn

        # Fills logical values
        powerCost[powerCost[,1] == i, 5] <- FALSE
        powerCost[powerCost[,1] == i, 6] <- FALSE
      }
    }

    # Tags the optimal case (minimum cost within the identified sampling effort)
    powerCost$Suggested <- powerCost$OptPower & powerCost$OptCost
    powerCost$Suggested <- ifelse(powerCost$Suggested == TRUE, "***", "")
    powerCost <- powerCost[,-c(5,6)]

  } else {
    # Calculating the cost
    powerCost <- powr[,-c(3,4)]
    powerCost$Cost <- powerCost[,1] * cn
    # Find the iterations that meet the desired power range
    powerCost$OptPower <- powerCost$Power >= objective

    # Safe search of the minimum cost among cases where objective is met
    costOK <- powerCost$Cost[powerCost$OptPower]

    if(length(costOK) == 0){
      warning("No power values above the specified precision.")
      powerCost[length(powr), 4] <- TRUE
    }
    # Find the minimum number of samples within the desired range
    powerCost$OptCost <- powerCost$Cost == min(powerCost$Cost[powerCost$OptPower == TRUE])
    # Tags the optimal case (minimum cost within the identified sampling effort)
    powerCost$Suggested <- powerCost$OptPower * powerCost$OptCost
    powerCost$Suggested <- ifelse(powerCost$Suggested == TRUE, "***", "")

    powerCost <- powerCost[,-c(4,5)]
  }

  return(powerCost)
}
