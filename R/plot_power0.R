#' Power curves for different sampling efforts
#'
#' \code{plot_power()} can be used to visualize the power of a study as a
#' function of the sampling effort. The power curve plot shows that the
#' power of the study increases as the sample size increases, and the density
#' plot shows the overlapping areas where \eqn{\alpha} and \eqn{\beta} are
#' significant.
#'
#' @param powr Part of the object of class "ecocbo_beta" that results from
#' [sim_beta()].
#' @param m Calculated in [plot_power()]. When using the `single.factor` model,
#' m is `NULL`.
#' @param n Calculated in [plot_power()].
#' @param cVar Calculated variation components.
#' @param model Model used for calculating power. Options, so far, are
#' 'single.factor' and 'nested.symmetric'.
#'
#' @return  Power curves for the different values of 'm'. The selected, or computed,
#' 'n' is marked in white with a bold outline.
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
#' [plot_power()]
#'
#' @importFrom graphics hist
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_label theme_bw
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous scale_color_manual
#' @importFrom ggplot2 scale_linetype_manual element_blank element_rect element_line
#' @importFrom ggplot2 theme guides geom_col geom_area geom_vline coord_cartesian
#' @importFrom ggplot2 annotate
#' @importFrom stats density
#' @importFrom rlang .data
#'
#' @keywords internal

## Power curve ----
power_curve <- function(powr, m = NULL, n, cVar, model){

  if(max(powr$n) <= 5){
    labx <- max(powr$n) - 1
  } else {
    labx <- max(powr$n) - 3
  }

  if(model == "single.factor"){
    nn <- n
    dummy <- data.frame(n = 1, Power = 0)
    powrPl <- data.frame(powr[,c(1:2)])
    powrPl <- rbind(dummy, powrPl)

    p1 <- ggplot(data = powrPl, aes(x = n, y = .data$Power))+
      geom_line()+
      geom_point(size = 1)+
      geom_point(data = powrPl[powrPl$n == nn,],
                 shape = 21, size = 2, stroke = 1.5, fill = "white")+
      annotate("label", x = labx, y = 0.1,
               label = paste0("n = ", nn))
  } else {
    mTot <- c(2:max(powr$m))
    mm <- m
    nn <- n
    dummy <- data.frame(m = mTot, n = 1, Power = 0)
    powrPl <- data.frame(powr[,c(1:3)])
    powrPl <- rbind(dummy, powrPl)
    powrPl$m <- factor(powrPl$m, ordered = TRUE)

    p1 <- ggplot(data = powrPl[powrPl$m == mm,], aes(x = n, y = .data$Power))+
      geom_line() +
      geom_point(size = 1, show.legend = FALSE)+
      geom_point(data = powrPl[powrPl$m == mm & powrPl$n == nn,],
                 shape = 21, size = 2, stroke = 1.5, fill = "white")+
      annotate("label", x = labx, y = 0.1,
               label = paste0("m = ", mm,
                              "\nn = ", nn,
                              "\nCV_B(A) = ", as.numeric(cVar[1])))
  }

  brks <- if(max(powr$n) <= 15){
    c(1:max(powr$n))
  } else {
    seq(0,max(powr$n),5)
  }

  p1 <- p1 +
    scale_linetype_manual(values = c(2,1))+
    scale_x_continuous(breaks = brks)+
    scale_y_continuous(name = "Power", limits = c(0, 1), breaks = seq(0, 1, 0.1),
                       minor_breaks=seq(0,1,0.05))+
    theme_bw()+
    theme(panel.grid.minor.x = element_blank(),
          panel.border = element_rect(linewidth = 0.4),
          axis.ticks= element_line(linewidth = 0.2),
          )+
    guides(alpha = "none", size = "none")

  return(p1)
}

#' Power curves for different sampling efforts
#'
#' \code{plot_power()} can be used to visualize the power of a study as a
#' function of the sampling effort. The power curve plot shows that the
#' power of the study increases as the sample size increases, and the density
#' plot shows the overlapping areas where \eqn{\alpha} and \eqn{\beta} are
#' significant.
#'
#' @param results Part of the object of class "ecocbo_beta" that results from
#' [sim_beta()].
#' @param powr Part of the object of class "ecocbo_beta" that results from
#' [sim_beta()].
#' @param m Calculated in [plot_power()]. When using the `single.factor` model,
#' m is `NULL`.
#' @param n Calculated in [plot_power()].
#' @param method Which plot is to be drawn? It is used to omit the text label when
#' the user selects `both` as method.
#' @param cVar Calculated variation components.
#' @param model Model used for calculating power. Options, so far, are
#' 'single.factor' and 'nested.symmetric'.
#' @param completePlot Logical. Is the plot to be drawn complete? If FALSE the plot will
#' be trimmed to present a better distribution of the density plot.
#'
#' @return  A density plot for the observed pseudoF values and a line marking
#' the value of pseudoF that marks the significance level indicated in [sim_beta()].
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
#' [plot_power()]
#'
#' @aliases densityplot
#'
#' @importFrom graphics hist
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_label theme_bw
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous scale_color_manual
#' @importFrom ggplot2 scale_linetype_manual element_blank element_rect element_line
#' @importFrom ggplot2 theme guides geom_col geom_area geom_vline coord_cartesian
#' @importFrom ggplot2 annotate
#' @importFrom stats density
#' @importFrom rlang .data
#'
#' @keywords internal

## Probability density curve ----
density_plot <- function(results, powr, m = NULL, n, method, cVar, model,
                         completePlot){

  if(model == "single.factor"){
    # Intersection point (FCrit)
    xIntersect <- powr[powr$n == n, 4]

    # Subset of results to get only the values with required n
    resultsPl <- results[results$n == n, 5:6]

    # Label for the plot
    cVarLabel <- paste0("n = ", n)
  } else {
    # intersection point (Fcrit)
    xIntersect <- powr[powr$m == m & powr$n == n,][,5]

    # Subset of results to get only the values with required m and n
    resultsPl <- results[results$m == m & results$n == n, 5:6]

    # Label for the plot
    cVarLabel <- paste0("m = ", m,
                        "\nn = ", n,
                        "\nCV_B(A) = ", as.numeric(cVar[1]))
  }

  # helper values for the histogram
  bottom <- min(resultsPl, na.rm = T)
  top <- max(resultsPl, na.rm = T)
  breaks <- seq(from = bottom, to = top, by = ((top-bottom) / 13))

  # Compute histograms and then the density values according to histogram.
  # Values for Y are normalized 0-1
  histH0 <- graphics::hist(resultsPl[,1], breaks = breaks, plot = F)
  histHa <- graphics::hist(resultsPl[,2], breaks = breaks, plot = F)
  topHist <- max(histH0$density, histHa$density, na.rm = T)
  bottHist <- min(histH0$density, histHa$density, na.rm = T)
  densHistH0 <- data.frame(x = histH0$mids,
                           y = ((histH0$density - bottHist) / (topHist - bottHist)))
  densHistHa <- data.frame(x = histHa$mids,
                           y = ((histHa$density - bottHist) / (topHist - bottHist)))

  # Compute density matrices form FHa and FH0.
  # Values for Y are normalized 0-1
  denHa <- stats::density(resultsPl$pseudoFHa, n = 128, adjust = 1.5, na.rm = T) # n = 128 as it has to be a power of 2
  denH0 <- stats::density(resultsPl$pseudoFH0, n = 128, adjust = 1.5, na.rm = T)
  densHa <- data.frame(x = 1:128, y = NA)
  densHa[,1] <- denHa[["x"]]
  densHa[,2] <- denHa[["y"]]
  densH0 <- data.frame(x = 1:128, y = NA)
  densH0[,1] <- denH0[["x"]]
  densH0[,2] <- denH0[["y"]]
  topY <- max(densHa$y, densH0$y, na.rm = T)
  bottY <- min(densHa$y, densH0$y, na.rm = T)
  densHa[,2] <- (densHa$y - bottY) / (topY - bottY)
  densH0[,2] <- (densH0$y - bottY) / (topY - bottY)

  # Plot
  p1 <- ggplot()+
    ggplot2::geom_col(data = densHistHa, aes(x = densHistHa[,1], y = densHistHa[,2]),
                      color = "#56B4E9", fill = "#56B4E9", alpha = 0.1)+
    ggplot2::geom_col(data = densHistH0, aes(x = densHistH0[,1], y = densHistH0[,2]),
                      color = "#E69F00", fill = "#E69F00", alpha = 0.1,)+
    ggplot2::geom_area(data = densHa, aes(x = densHa[,1], y = densHa[,2]),
                       fill = "#56B4E9", color = "#56B4E9",
                       alpha = 0.5)+
    ggplot2::geom_area(data = densH0, aes(x = densH0[,1], y = densH0[,2]),
                       fill = "#E69F00", color = "#E69F00",
                       alpha = 0.5)+
    geom_vline(xintercept = xIntersect, linetype = 2)+
    geom_label(aes(x = xIntersect + 1, y = 0.95,
                   label = paste0("Fcrit = ",round(xIntersect, 2))),
               alpha = 0.2, nudge_x = 0.1)+
    theme_bw()+
    scale_x_continuous(name = "pseudoF",
                       breaks = function(x)seq(floor(bottom), ceiling(top)))+
    scale_y_continuous(name = "Density", limits = c(0, 1), breaks = seq(0, 1, 0.2))+
    # coord_cartesian(xlim = c(bottom,top))+
    ggplot2::theme(panel.grid.minor = element_blank(),
                   panel.border = element_rect(linewidth = 0.4),
                   axis.ticks= element_line(linewidth = 0.2))

  if(completePlot == FALSE){
    ftop <- densHistHa[8:13,]
    ftop <- as.numeric(rownames(ftop[rank(ftop[,2], ties.method="first")[1],]))

    p1 <- p1 +
      coord_cartesian(xlim = c(floor(bottom), ceiling(densHistHa[ftop,1])))
    if(method == "density"){
      p1 <- p1 +
        geom_label(aes(x = densHistHa[ftop,1] - 1, y = 0.9, label = cVarLabel))
    }
  } else {
    p1 <- p1 +
      coord_cartesian(xlim = c(0, ceiling(top)))
    if(method == "density"){
      p1 <- p1 +
        geom_label(aes(x = top - 1, y = 0.9, label = cVarLabel))
    }
  }

  # if(method == "density"){
  #   p1 <- p1 +
  #     geom_label(aes(x = top - 1, y = 0.9, label = cVarLabel))
  # }
  return(p1)
}

#' Power surface for different sampling efforts
#'
#' \code{plot_power()} can be used to visualize the power of a study as a
#' function of the sampling effort. The power curve plot shows that the
#' power of the study increases as the sample size increases, and the density
#' plot shows the overlapping areas where \eqn{\alpha} and \eqn{\beta} are
#' significant.
#'
#' @param powr Part of the object of class "ecocbo_beta" that results from
#' [sim_beta()].
#' @param model Model used for calculating power. Options, so far, are
#' 'single.factor' and 'nested.symmetric'.
#'
#' @return  A surface plot for the observed statistical power at different sampling
#' efforts, as indicated in [sim_beta()].
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
#' [plot_power()]
#'
#' @aliases surfaceplot
#'
#' @importFrom plotly plot_ly add_surface layout
#'
#' @keywords internal

## Probability density curve ----
surface_plot <- function(powr, model){

  if(model == "single.factor"){
    message("The surface plot is available only for nested factors experiments.")
  } else {
    # Getting the unique dimensions
    ms <- sort(unique(powr$m))
    ns <- sort(unique(powr$n))

    # Saving powr to a matrix
    Z <- matrix(powr$Power,
                nrow = length(ns),
                ncol = length(ms),
                byrow = FALSE)

    # Plotting the surface with plotly
    p1 <- plotly::plot_ly(x = ~ms,
                    y = ~ns,
                    z = ~Z,
                    contours = list(z = list(show = TRUE,
                                             start = 0, end = 1,
                                             size = 0.05,
                                             highlightcolor = "red"))) |>
      plotly::add_surface() |>
      plotly::layout(scene = list(
        xaxis = list(title = "m"),
        yaxis = list(title = "n"),
        zaxis = list(title = "Power", nticks = 5)
      ))
  }

  return(p1)
}

