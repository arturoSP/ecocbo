# Función para determinar los componentes de variación previo al modelo de
# optimización de costo-beneficio

scompvar <- function(data, m = NULL, n = NULL){
  # lee la tabla de resultados desde la lista de Beta para poder usar los MeanSquares de Ha
  resultsBeta <- data$Results
  if(is.null(m)){m <- max(resultsBeta$m)}
  if(is.null(n)){n <- max(resultsBeta$n)}
  resultsBeta <- resultsBeta[resultsBeta$m == m & resultsBeta$n == n,]
  
  # Crea un data frame con tamaño: número juegos de datos simulados x 3
  compVar <- data.frame(compVarA = NA, compVarR = NA)
  
  # Calcula componentes de variación promedio (de la Tabla 9.3)
  # σe = RMS
  compVar[,2] <- mean(resultsBeta[,8], na.rm = T)
  # σB(A) = (AMS - σe) / n
  MSA <- mean(resultsBeta[,7], na.rm = T)
  compVar[,1] <- (MSA - compVar[,2]) / n
  
  return(compVar)
}


