#' Plot Statistical Power and Pseudo-F Distributions
#'
#' Visualizes the statistical power of a study as a function of the sampling effort.
#' The power curve plot illustrates how power increases with sample size, while
#' the density plot highlights overlapping areas where \eqn{\alpha} and
#' \eqn{\beta} are significant.
#'
#' @param data Object of class `"ecocbo_beta"` obtained from [sim_beta()].
#' @param cbo Optional. Object of class `"cbo_result"` obtained from [sim_cbo()].
#' If this is included, `plot_power()` uses the optimal values that have been already
#' calculated.
#' @param m Optional. Integer. Number of replicates `m` to use for power computation.
#' Defaults to `NULL`, in which case the function selects the number of sites that
#' result in a sampling effort that is close to \eqn{1 - \alpha}.
#' @param n Optional. Integer. Number of samples `n` within the selected `m`.
#' Defaults to `NULL`, and the function selects the number of samples yielding a
#' power close to \eqn{1 - \alpha}.
#' @param method Character. Type of plot to generate:
#'   - "power": Plots the power curve.
#'   - "density": Plots the density distribution of pseudo-F values.
#'   - "both": Displays both plots side by side.
#'   - "surface": Displays a 3d surface plot of the power curves for nested
#'   factors experiments.
#' @param completePlot Logical. Is the plot to be drawn complete? If TRUE the
#' plot will be trimmed to present a better distribution of the density plot.
#'
#' @return A plot displaying:
#'   - If `method = "power"`, power curves for different values of `m`, with the
#'   selected `n` highlighted in red.
#'   - If `method = "density"`: a density plot of observed pseudo-F values with
#'   a vertical line indicating significance from [sim_beta()].
#'   - If `method = "both"`: a composite figure with both the power curve and the
#'   density plot.
#'   - If `method = "surface"`: a surface plot for the statistical power in different
#'   sampling designs.
#'
#' The selected values of `m`, `n`, and the corresponding component of variation
#' are displayed in all cases.
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
#' @importFrom plotly plot_ly add_surface layout
#' @importFrom stats density
#' @importFrom rlang .data
#'
#' @examples
#' # Power curve visualization
#' plot_power(data = epiBetaR, method = "power")
#'
#' # Density plot of pseudo-F values
#' plot_power(data = betaNested, method = "density")
#'
#' # Composite plot with both power curve and density plot
#' plot_power(data = betaNested, method = "both")
#'

plot_power <- function(data, cbo = NULL, n = NULL, m = NULL,
                       method = "power", completePlot = TRUE){
# Función para graficar curvas de frecuencia de F para H0 y Ha ----
  ## Reading data ----
  if(!inherits(data, "ecocbo_beta")){
    stop("data is not the right class (\"ecocbo_beta\")")
  }

  if(!is.null(cbo)){
    if(!inherits(cbo, "cbo_result")){
      stop("cbo data is not the right class (\"cbo_result\")")
    }
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
      if(is.null(cbo)){
        # Find optimal value for n
        n <- powr[which.min(abs(powr$Power - (1 - alpha))),1]
      } else {
        # Use the optimal value
        n <- cbo[cbo$OptCost == TRUE, 1]
      }

    } else {
      # Validating n provided by user
      if(ceiling(n) != floor(n)){stop("n must be integer")}
      if(n <= 1){stop("n must be larger than 1")}
      if(n > max(powr$n)){stop("n is larger than the simulated n value")}
    }
    if(method == "surface"){
      stop("Surface is only available for nested factors experiments.")
    }

  } else {                                 ### Double-factor-model validation ----
    if(is.null(m)){
      if(is.null(cbo)){
        # Find optimal value for m
        m <- powr[which.min(abs(powr$Power - (1 - alpha))),1]
      } else {
        # Use the optimal m value
        m <- cbo[cbo$OptCost == TRUE, 1]
      }
    } else {
      # Validating m provided by user
      if(ceiling(m) != floor(m)){stop("m must be integer")}
      if(m <= 1){stop("m must be larger than 1")}
      if(m > max(powr$m)){stop("m is langer than the simulated m value")}
    }

    if(is.null(n)){
      if(is.null(cbo)){
        # Find optimal value for n
        powm <- powr[powr$m == m,]
        n <- powm[which.min(abs(powm$Power - (1 - alpha))),2]
        remove(powm)
      } else {
        # Use the optimal n value
        n <- cbo[cbo$OptCost == TRUE, 2]
      }

    } else {
      # Validating n provided by user
      if(ceiling(n) != floor(n)){stop("n must be integer")}
      if(n <= 1){stop("n must be larger than 1")}
      if(n > max(powr$n)){stop("n is larger than the simulated n value")}
    }
  }

  # Validating selected method
  if(method != "both" & method != "power" & method != "density" & method != "surface"){
    stop("Available methods are \"power\", \"density\", \"both\" and \"surface\"")
  }

  ## Plot according to the parameters ----
  dataRes <- list(Results = results, model = data$model)
  class(dataRes) <- "ecocbo_data"
  cVar <- round(scompvar(dataRes, n, m)[,2],2)

  if(method == "both") {
    p1 <- power_curve(powr, m, n, cVar, model)
    p2 <- density_plot(results, powr, m, n, method, cVar = NULL, model, completePlot)
    plotF <- ggpubr::ggarrange(p1, p2)
  } else if(method == "power") {
    plotF <- power_curve(powr, m, n, cVar, model)
  } else if(method == "density") {
    plotF <- density_plot(results, powr, m, n, method, cVar, model, completePlot)
  } else if(method == "surface") {
    plotF <- surface_plot(powr, model)
  }
  return(plotF)
}
