# Función para graficar curvas de frecuencia de F para H0 y Ha ----
library(ggplot2)

# Helper functions ----
## Power curve ----
power_curve <- function(powr, m, n){
  dummy <- data.frame(m = m, n = 1, power = 0, beta = NA,
                      fCrit = NA, AMSHa = NA, RMSHa = NA)
  powr <- rbind(powr, dummy)
  powr[powr$m == m,] |>
    ggplot(mapping = aes(x = n, y = power))+
    geom_line()+
    geom_point(shape = 1, size = 2)+
    geom_point(data = powr[powr$m == m & powr$n == n,], color = "red")+
    geom_label(aes(x = 1, y = 1, label = paste0("m = ", m)))+
    theme_bw()+
    scale_x_continuous(name = "Sample size", breaks = function(x)
      seq(0, max(powr$n)))+
    scale_y_continuous(name = "Power")
}

## Probability density curve ----
density_plot <- function(results, powr, m, n){
  xIntersect <- powr[powr$m == m & powr$n == n,][,5]
  resultsPl <- results[results$m == m & results$n == n, 5:6]
  bottom <- min(resultsPl)
  top <- max(resultsPl)
  breaks <- seq(from = bottom, to = top, by = ((top-bottom) /50))
  histH0 <- hist(resultsPl[,1], breaks = breaks, plot = F)
  histH0 <- hist(resultsPl[,1], breaks = 10, plot = F)
  densH0 <- data.frame(x = histH0$mids, y = histH0$density / sum(histH0$density))
  histHa <- hist(resultsPl[,2], breaks = breaks, plot = F)
  histHa <- hist(resultsPl[,2], breaks = 10, plot = F)
  densHa <- data.frame(x = histHa$mids, y = histHa$density / sum(histHa$density))

  ggplot2::ggplot()+
    geom_histogram(data = densHa, aes(x = x, y = y),
                   stat = "identity", color = "#56B4E9", fill = "#56B4E9", alpha = 0.1,
                   bins = 50)+
    geom_histogram(data = densH0, aes(x = x, y = y),
                   stat = "identity", color = "#E69F00", fill = "#E69F00", alpha = 0.1,
                   bins = 50)+
    geom_density(data = resultsPl, aes(x = pseudoFHa, y = ..density../sum(histHa$density)),
                 fill = "#56B4E9", color = "#56B4E9",
                 alpha = 0.5, adjust = 1.5)+
    geom_density(data = resultsPl, aes(x = pseudoFH0, y = ..density../sum(histH0$density)),
                 fill = "#E69F00", color = "#E69F00",
                 alpha = 0.5, adjust = 1.5)+
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
plot_power <- function(data, m, n = NULL, method = "both"){
  # Reading data ----
  powr <- data[["Power"]]
  results <- data[["Results"]]

  # Plot according to the parameters ----
  if(is.null(n)){
    n <- max(powr[powr$m == m & powr$power<1,][2])
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
