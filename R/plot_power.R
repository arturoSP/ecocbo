#' Power curves for different sampling efforts
#'
#' \code{plot_power()} can be used to visualize the power of a study as a
#' function of the sampling effort. The power curve plot shows that the
#' power of the study increases as the sample size increases, and the density
#' plot shows the overlapping areas where \eqn{\alpha} and \eqn{\beta} are
#' significant.
#'
#' @param data Object of class "ecocbo_beta" that results from [sim_beta()].
#' @param m Defaults to NULL, and then the function computes the number of
#' sites 'm' that result in a sampling effort that is close to (1 - alpha) in
#' power. If provided, said number of site will be used.
#' @param n Defaults to NULL, and then the function computes the number of
#' samples 'n', within the selected 'm', that result in a sampling effort close
#' to (1 - alpha) in power. If provided, said number of samples will be used.
#' @param method The desired plot. Options are "power", "density" or "both".
#' "power" plots the power curve, "density" plots the density distribution of
#' pseudoF, and "both" draws both plots one next to the other.
#'
#' @return  If the method is "power", then the power curves for the different values
#' of 'm'. The selected, or computed, 'n' is marked in red. If the method is "density", then a
#' density plot for the observed pseudoF values and a line marking the value of
#' pseudoF that marks the significance level indicated in [sim_beta()].
#' If the method is "both", then a composite with power curves and a
#' density plot side by side.
#'
#' The value of the selected 'm', 'n' and the corresponding component of variation
#' are presented in all methods.
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
#' [scompvar()]
#' [sim_cbo()]
#' [prep_data()]
#'
#' @aliases plotpower
#'
#' @export
#' @importFrom graphics hist
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_label theme_bw
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous scale_color_manual
#' @importFrom ggplot2 scale_linetype_manual element_blank element_rect element_line
#' @importFrom ggplot2 theme guides geom_col geom_area geom_vline coord_cartesian
#' @importFrom ggplot2 annotate
#' @importFrom stats density
#' @importFrom rlang .data
#'
#' @examples
#' epiBetaR <- sim_beta(simResults, alpha = 0.05)
#'
#' plot_power(data = epiBetaR, n = NULL, m = 3, method = "power")
#' plot_power(data = epiBetaR, n = NULL, m = 3, method = "density")
#' plot_power(data = epiBetaR, n = 4, m = 3, method = "both")

plot_power <- function(data, n = NULL, m = NULL, method = "power"){
# Función para graficar curvas de frecuencia de F para H0 y Ha ----
  ## Reading data ----
  if(!inherits(data, "ecocbo_beta")){
    stop("data is not the right class(\"ecocbo_beta\")")
  }

  powr <- data[["Power"]]
  results <- data[["Results"]]
  alpha <- data[["alpha"]]
  model <- data[["model"]]

  ## Validating data  ----
    if(model == "single.factor"){     ### Single-factor-model validation ----
    # Cancelling m value the user could have provided
    m <- NULL

    if(is.null(n)){
      # Find optimal value for n
      n <- powr[which.min(abs(powr$Power - (1 - alpha))),1]
    } else {
      # Validating n provided by user
      if(ceiling(n) != floor(n)){stop("n must be integer")}
      if(n <= 1){stop("n must be larger than 1")}
      if(n > max(powr$n)){stop("n is larger than the simulated n value")}
    }

  } else {                                 ### Double-factor-model validation ----
    if(is.null(m)){
      # Find optimal value for m
      m <- powr[which.min(abs(powr$Power - (1 - alpha))),1]

    } else {
      # Validating m provided by user
      if(ceiling(m) != floor(m)){stop("m must be integer")}
      if(m <= 1){stop("m must be larger than 1")}
      if(m > max(powr$m)){stop("m is langer than the simulated m value")}
    }

    if(is.null(n)){
      # Find optimal value for n
      powm <- powr[powr$m == m,]
      n <- powm[which.min(abs(powm$Power - (1 - alpha))),2]
      remove(powm)
    } else {
      # Validating n provided by user
      if(ceiling(n) != floor(n)){stop("n must be integer")}
      if(n <= 1){stop("n must be larger than 1")}
      if(n > max(powr$n)){stop("n is larger than the simulated n value")}
    }
  }

  # Validating selected method
  if(method != "both" & method != "power" & method != "density"){
    stop("Available methods are \"power\", \"density\" and \"both\"")
  }

  ## Plot according to the parameters ----
  dataRes <- list(Results = results, model = data$model)
  class(dataRes) <- "ecocbo_data"
  cVar <- round(scompvar(dataRes, n, m)[,2],2)

  if(method == "both") {
    p1 <- power_curve(powr, m, n, cVar, model)
    p2 <- density_plot(results, powr, m, n, method, cVar = NULL, model)
    plotF <- ggpubr::ggarrange(p1, p2)
  } else if(method == "power") {
    plotF <- power_curve(powr, m, n, cVar, model)
  } else if(method == "density") {
    plotF <- density_plot(results, powr, m, n, method, cVar, model)
  }
  return(plotF)
}
