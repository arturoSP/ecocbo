#' Power curves for different sampling efforts
#'
#' @param data List that results from beta().
#' @param m Site label to be used as basis for the plot.
#' @param n Defaults to NULL, and then the function computes the number of samples (n) that results in a sampling effort close to 95% in power. If provided, said number of samples will be used.
#' @param method The desired plot. Options are "power", "density" or "both". "power" plots the power curve, "density" plots the density distribution of pseudoF, and "both" draws both plots one next to the other.
#'
#' @return A plot
#' @export
#' @importFrom graphics hist
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_label
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 geom_col
#' @importFrom ggplot2 geom_area
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 coord_cartesian
#' @importFrom stats density
#'
#' @examples
#' plot_power(data = epiBetaR, n = 4, m = 2, method = "both")
#' plot_power(data = epiBetaR, n = NULL, m = 3, method = "both")

plot_power <- function(data, n = NULL, m, method = "both"){
# Función para graficar curvas de frecuencia de F para H0 y Ha ----

# Helper functions ----
## Power curve ----
power_curve <- function(powr, m, n){
  dummy <- data.frame(m = m, n = 1, Power = 0, Beta = NA,
                      fCrit = NA, AMSHa = NA, RMSHa = NA)
  powrPl <- powr[powr$m == m,]
  powrPl <- rbind(powrPl, dummy)
  ggplot(data = powrPl,
         aes(x = powrPl[, "n"], y = powrPl[, "Power"]))+
    geom_line()+
    geom_point(data = powrPl[powrPl$n != 1,],
               aes(x = powrPl[powrPl$n != 1, "n"],
                   y = powrPl[powrPl$n != 1, "Power"]),
               shape = 1, size = 2)+
    geom_point(data = powrPl[powrPl$n == n,],
               aes(x = powrPl[powrPl$n == n, "n"],
                   y = powrPl[powrPl$n == n, "Power"]), color = "red")+
    geom_label(aes(x = 1, y = 1, label = paste0("m = ", m)))+
    theme_bw()+
    scale_x_continuous(name = "Sample size", breaks = function(x)seq(0, max(powrPl$n, na.rm = T)))+
    scale_y_continuous(name = "Power", limits = c(0, 1))
}

## Probability density curve ----
density_plot <- function(results, powr, m, n, method){
  # intersection point (Fcrit)
  xIntersect <- powr[powr$m == m & powr$n == n,][,5]

  # Subset of results to get only the values with required m and n
  resultsPl <- results[results$m == m & results$n == n, 5:6]

  # helper values for the histogram
  bottom <- min(resultsPl, na.rm = T)
  top <- max(resultsPl, na.rm = T)
  breaks <- seq(from = bottom, to = top, by = ((top-bottom) / 13)) # 1/10 of the n used for density

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

  # Compute density matrixes form FHa and FH0.
  # Values for Y are normalized 0-1
  denHa <- density(resultsPl$pseudoFHa, n = 128, adjust = 1.5, na.rm = T) # n = 128 as it has to be a power of 2
  denH0 <- density(resultsPl$pseudoFH0, n = 128, adjust = 1.5, na.rm = T)
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
    geom_label(aes(x = xIntersect, y = 0,
                   label = paste0("Fcrit = ",round(xIntersect, 2))),
               alpha = 0.2, nudge_x = 0.1)+
    theme_bw()+
    scale_x_continuous(name = "pseudoF",
                       breaks = function(x)seq(floor(bottom), ceiling(top)))+
    scale_y_continuous(name = "Density")+
    coord_cartesian(xlim = c(bottom,top))
  if(method == "density"){
    p1 <- p1 +
      geom_label(aes(x = bottom, y = 1, label = paste0("m = ", m)))
  }
  return(p1)
}

# Main function ----
# plot_power <- function(data, n = NULL, m, method = "both")
  ## Reading data ----
  if(class(data) != "ecocbo_beta")
    stop("data is not the right class(\"ecocbo_beta\")")

  powr <- data[["Power"]]
  results <- data[["Results"]]

  ## Validating data  ----
  if(ceiling(m) != floor(m)){stop("m must be integer")}
  if(m <= 1){stop("m must be larger than 1")}
  if(m > max(powr$m)){stop("m is langer than the simulated m value")}

  if(is.null(n)){n <- max(powr[powr$m == m & powr$Power<1,][2])}
  if(ceiling(n) != floor(n)){stop("n must be integer")}
  if(n <= 1){stop("n must be larger than 1")}
  if(n > max(powr$n)){stop("n is larger than the simulated n value")}

  if(method != "both" & method != "power" & method != "density")
    stop("Available methods are \"power\", \"density\" and \"both\"")

  ## Plot according to the parameters ----
  if(method == "both") {
    p1 <- power_curve(powr, m, n)
    p2 <- density_plot(results, powr, m, n, method)
    plotF <- ggpubr::ggarrange(p1, p2)
  } else if(method == "power") {
    plotF <- power_curve(powr, m, n)
  } else if(method == "density") {
    plotF <- density_plot(results, powr, m, n, method)
  }
  return(plotF)
}
