# función para calcular el model de optimización de costo-beneficio

# funciones de apoyo ----
# función para determinar el número de sitios a muestrear, conociendo
# costos de muestreo y número de muestras a tomar.
cost_n <- function(n, ct, ck, cj){
  b <- data.frame(nOpt = n, bOpt = NA)
  
  # aplica la ecuación 9.19 de Underwood, 1997
  b[,2] <- round(ct / (n * ck + cj))
  b[,1] <- round(b[,1])
  
  return(b)
}
  
# función para determinar el número de sitios a muestrear, conociendo
# la variabilidad desde sampsd()
cost_v <- function(comp.var, multSE, n){
  b <- data.frame(nOpt = n, bOpt = NA)
  
  # Aplica la ecuación 9.18 de Underwood, 1997
  b[,2] <- round((comp.var$compVarR + b$nOpt * comp.var$compVarA) / (multSE * multSE * b$nOpt))
  b[,1] <- round(b[,1])
  
  return(b)
}

# función principal ----
simcbo <- function(comp.var, multSE, ct, ck, cj){
  # Cálculo de n óptimo ----
  nOpt <- sqrt((cj * comp.var$compVarR) / (ck * comp.var$compVarA))
  
  # Cálculo de b óptimo ----
  if(is.null(multSE)) {
    b <- cost_n(nOpt, ct, ck, cj)
  } else {
    b <- cost_v(comp.var, multSE, nOpt)
  }
  return(b)
}

