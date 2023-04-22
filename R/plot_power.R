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
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 geom_area
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_label
#' @importFrom ggplot2 coord_cartesian
#' @importFrom stats density
#'
#' @examples
#' library("SSP")
#' # Load data and adjust it.
#' epiH0 <- epiDat
#' epiH0[,"site"] <- as.factor("T0")
#' epiHa <- epiDat
#' epiHa[,"site"] <- as.factor(epiHa[,"site"])
#'
#' # Calculate simulation parameters.
#' parH0 <- SSP::assempar(data = epiH0, type = "counts", Sest.method = "average")
#' parHa <- SSP::assempar(data = epiHa, type = "counts", Sest.method = "average")
#'
#' # Simulation.
#' simH0 <- SSP::simdata(parH0, cases = 3, N = 1000, sites = 1)
#' simHa <- SSP::simdata(parHa, cases = 3, N = 100, sites = 10)
#'
#' betaResults <- beta(simH0, simHa, n = 10, m = 3, k = 50, alpha = 0.05)
#'
#' plot_power(data = betaResults, m = 2, n = 7, method = "both")
#' plot_power(data = betaResults, m = 3, n = NULL, method = "both")

plot_power <- function(data, m, n = NULL, method = "both"){
# Función para graficar curvas de frecuencia de F para H0 y Ha ----

# Helper functions ----
## Power curve ----
power_curve <- function(powr, m, n){
  dummy <- data.frame(m = m, n = 1, Power = 0, Beta = NA,
                      fCrit = NA, AMSHa = NA, RMSHa = NA)
  powr <- rbind(powr, dummy)
  ggplot(data = powr[powr$m == m,],
         aes(x = powr[powr$m == m, "n"], y = powr[powr$m == m, "Power"]))+
    geom_line()+
    geom_point(shape = 1, size = 2)+
    geom_point(data = powr[powr$m == m & powr$n == n,],
               aes(x = powr[powr$m == m & powr$n == n, "n"],
                   y = powr[powr$m == m & powr$n == n, "Power"]), color = "red")+
    geom_label(aes(x = 1, y = 1, label = paste0("m = ", m)))+
    theme_bw()+
    scale_x_continuous(name = "Sample size", breaks = function(x)
      seq(0, max(powr$n)))+
    scale_y_continuous(name = "Power")
}

## Probability density curve ----
density_plot <- function(results, powr, m, n){
  # intersection point (Fcrit)
  xIntersect <- powr[powr$m == m & powr$n == n,][,5]

  # Subset of results to get only the values with required m and n
  resultsPl <- results[results$m == m & results$n == n, 5:6]

  # helper values for the histogram
  bottom <- min(resultsPl)
  top <- max(resultsPl)
  breaks <- seq(from = bottom, to = top, by = ((top-bottom) / 13))

  # Compute histograms and then the density values according to histogram. Values in Y are normalized 0-1
  histH0 <- graphics::hist(resultsPl[,1], breaks = breaks, plot = F)
  histHa <- graphics::hist(resultsPl[,2], breaks = breaks, plot = F)
  densHistH0 <- data.frame(x = histH0$mids, y = ((histH0$density - min(histH0$density, histHa$density))/(max(histH0$density, histHa$density) - min(histH0$density, histHa$density))))
  densHistHa <- data.frame(x = histHa$mids, y = ((histHa$density - min(histH0$density, histHa$density))/(max(histH0$density, histHa$density) - min(histH0$density, histHa$density))))

  # Compute density matrixes form FHa and FH0. Values for Y are normalized 0-1
  denHa <- density(resultsPl$pseudoFHa, n = 128, adjust = 1.5)
  denH0 <- density(resultsPl$pseudoFH0, n = 128, adjust = 1.5)
  densHa <- data.frame(x = 1:128, y = NA, yn = NA)
  densHa[,1] <- denHa[["x"]]
  densHa[,2] <- denHa[["y"]]
  densH0 <- data.frame(x = 1:128, y = NA, yn = NA)
  densH0[,1] <- denH0[["x"]]
  densH0[,2] <- denH0[["y"]]
  densHa[,3] <- (densHa$y - min(densHa$y, densH0$y)) / (max(densHa$y, densH0$y) - min(densHa$y, densH0$y))
  densH0[,3] <- (densH0$y - min(densHa$y, densH0$y)) / (max(densHa$y, densH0$y) - min(densHa$y, densH0$y))

  # Plot
  ggplot()+
    geom_histogram(data = densHistHa, aes(x = densHistHa$x, y = densHistHa$y),
                   stat = "identity", color = "#56B4E9", fill = "#56B4E9", alpha = 0.1,
                   bins = 50)+
    geom_histogram(data = densHistH0, aes(x = densHistH0$x, y = densHistH0$y),
                   stat = "identity", color = "#E69F00", fill = "#E69F00", alpha = 0.1,
                   bins = 50)+
    ggplot2::geom_area(data = densHa, aes(x = densHa$x, y = densHa$yn),
                       fill = "#56B4E9", color = "#56B4E9",
                       alpha = 0.5)+
    ggplot2::geom_area(data = densH0, aes(x = densH0$x, y = densH0$yn),
                       fill = "#E69F00", color = "#E69F00",
                       alpha = 0.5)+
    geom_vline(xintercept = xIntersect, linetype = 2)+
    geom_label(aes(x = xIntersect, y = 0,
                   label = paste0("Fcrit = ",round(xIntersect, 2))),
               alpha = 0.2, nudge_x = 0.1)+
    theme_bw()+
    scale_x_continuous(name = "pseudoF", breaks = function(x)
      seq(min(floor(resultsPl[,1])),
          max(ceiling(resultsPl[,2]))))+
    scale_y_continuous(name = "Density")+
    coord_cartesian(xlim = c(bottom,top))
}

# Main function ----
#plot_power <- function(data, m, n = NULL, method = "both"){
  # Reading data ----
  powr <- data[["Power"]]
  results <- data[["Results"]]

  # Plot according to the parameters ----
  if(is.null(n)){
    n <- max(powr[powr$m == m & powr$Power<1,][2])
    if(method == "both") {
      p1 <- power_curve(powr, m, n)
      p2 <- density_plot(results, powr, m, n)
      plotF <- ggpubr::ggarrange(p1, p2)
    } else if(method == "power") {
      plotF <- power_curve(powr, m, n)
    } else if(method == "density") {
      plotF <- density_plot(results, powr, m, n)
    }
  } else {
    if(method == "both") {
      p1 <- power_curve(powr, m, n)
      p2 <- density_plot(results, powr, m, n)
      plotF <- ggpubr::ggarrange(p1, p2)
    } else if(method == "power") {
      plotF <- power_curve(powr, m, n)
    } else if(method == "density") {
      plotF <- density_plot(results, powr, m, n)
    }
  }
  return(plotF)
}
