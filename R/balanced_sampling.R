balanced_sampling <- function(i, Y, mm, nn, YPU, H0Sim, HaSim, resultsHa, transformation, method){
  # Get the samples index
  sel <- sampling::balancedtwostage(Y, selection = 1, m = mm[i],
                                    n = nn[i], PU = YPU, FALSE)
  ones <- which(sel[,1] %in% 1)
  y0 <- H0Sim[ones,,resultsHa[i,1]]
  ya <- HaSim[ones,,resultsHa[i,1]]
  yHa <- dim(y0)[2] - 2

  # Apply PERMANOVA to get F and mean squares
  result1 <- permanova_oneway(y0[, 1:yHa], y0[,yHa+2],
                              transformation = transformation, method = method)
  result2 <- permanova_oneway(ya[, 1:yHa], y0[,(yHa+2)],
                              transformation = transformation, method = method)
  result0 <- matrix(nrow = 1, ncol = 4)
  colnames(result0) <- c("Fobs", "Fobs", "AMS", "RMS")

  # Gather the results and return
  result0[,1] <- result1[,3]
  result0[,2] <- result2[,3]
  result0[,3] <- result2[,1]
  result0[,4] <- result2[,2]
  return(result0)
}
