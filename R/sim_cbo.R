#' Simulated Cost-Benefit Optimization
#'
#' Applies a cost-benefit optimization based on a desired level of statistical
#' power and the sampling cost.
#'
#' @param data Object of class `"ecocbo_beta"` obtained from [sim_beta()].
#' @param cm Numeric. Cost per replicate.
#' @param cn Numeric. Cost per sampling unit.
#'
#' @return A data frame containing the optimized values for \code{m} number of
#' sites to sample and \code{n} number of samples per site.
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
#' compVar <- scompvar(data = simResults)
#'
#' # Optimization of single factor experiment
#' sim_cbo(data = epiBetaR, cn = 80)
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
    powerCost$Optimal <- powerCost$OptPower & powerCost$OptCost
    powerCost$Optimal <- ifelse(powerCost$Optimal == TRUE, "***", "")
    powerCost <- powerCost[,-c(5,6)]

    # # Subsetting the power dataframe to get data that is within the range:
    # # (1 - alpha) <--> 1
    # ideal <- powr[powr[,3] >= objective, -c(4,5)]
    #

    #
    #
    #

    #


    #
    # # Finds the minimum cost and marks it with an asterisc
    # powerCost$Optimal <- powerCost$Cost == min(powerCost$Cost)
    # powerCost$Optimal <- ifelse(powerCost$Optimal == TRUE, "***", "")

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
    powerCost$Optimal <- powerCost$OptPower * powerCost$OptCost
    powerCost$Optimal <- ifelse(powerCost$Optimal == TRUE, "***", "")

    powerCost <- powerCost[,-c(4,5)]
  }

  return(powerCost)
}
