#' dbMANOVA / PERMANOVA de una vía (factor único)
#'
#' @param x        Matriz/data.frame de comunidad (muestras x especies).
#'                 Debe ser el subset ya muestreado.
#' @param factEnv  data.frame o vector con el factor (primera columna = factor A).
#' 
#' @param transformation  Transformación para reducir dominancia. Una de: "none", "square root", "fourth root", "Log (X+1)".
#'                 Por defecto "sqrt".
#' @param method   Distancia para vegan::vegdist. Default "bray".
#' @param model    Etiqueta del modelo; aquí "single.factor".
#' @param return   "table","list","both".
#'
#' @return data.frame estilo ANOVA y/o lista con SS, df, MS, pseudoF.
#' @importFrom vegan vegdist
dbmanova_oneway <- function(x,
                            factEnv,
                            transformation = c("none", "square root", "fourth root", "Log (X+1)"),
                            method = "bray",
                            model = c("single.factor"),
                            return = c("table","list","both")) {

  # --- Validaciones ---
  transformation <- match.arg(transformation)
  method <- match.arg(method)
  model <- match.arg(model)
  return <- match.arg(return)

  if (!is.matrix(x) && !is.data.frame(x)) stop("`x` must be matrix/data.frame.")
  Xcomm <- as.matrix(x)
  if (any(!is.finite(Xcomm))) stop("`x` contains non finite values (NA/NaN/Inf).")

  # factEnv: aceptamos vector o data.frame; tomamos la primera col como factor
  if (is.null(dim(factEnv))) {
    fac <- factor(factEnv, exclude = NULL)
    nameA <- deparse(substitute(factEnv))
  } else {
    factEnv <- as.data.frame(factEnv)
    if (ncol(factEnv) < 1L) stop("`factEnv` debe tener al menos 1 columna (factor).")
    fac <- factor(factEnv[[1]], exclude = NULL)
    nameA <- if (!is.null(colnames(factEnv))) colnames(factEnv)[1] else "A"
  }

  if (nrow(Xcomm) != length(fac)) stop("`x` y `factEnv` must have the same number of rows.")

  # --- Transformaciones para atenuar dominancia ---
  transform_comm <- function(M, how = "square root") {
    how <- match.arg(how, c("none","square root","fourth root","Log (X+1)"))
    if (how == "none") return(M)
    if (how == "square root") return(sqrt(M))
    if (how == "fourth root") return(M^(1/4))
    if (how == "Log (X+1)") return(log1p(M))
  }
  Xtr <- transform_comm(Xcomm, transformation)

  # --- Distance and PCoA (Gower) ---
  d  <- vegan::vegdist(Xtr, method = method)
  Dm <- as.matrix(d)

  A  <- -0.5 * (Dm^2)
  n  <- nrow(A)
  I  <- diag(n)
  J  <- matrix(1/n, n, n)
  G  <- (I - J) %*% A %*% (I - J)

  EG <- eigen(G, symmetric = TRUE)
  sign_custom <- function(x) ifelse(x >= 0, 1, -1)
  eig_sign <- sign_custom(EG$values)

  Lhalf <- sqrt(abs(EG$values))
  Xpcoa <- sweep(EG$vectors, 2, Lhalf, `*`)   # coordinates for PCoA (n x k)

  # --- Diseño (una vía) ---
  H    <- model.matrix(~ fac - 1)     # n x a
  a    <- nlevels(fac)
  P_A  <- H %*% solve(crossprod(H)) %*% t(H)
  P_T  <- matrix(1/n, n, n)
  I_n  <- diag(n)

  # --- SSCP in PCoA space ---
  SSAm <- t(Xpcoa) %*% (P_A - P_T) %*% Xpcoa
  SSRm <- t(Xpcoa) %*% (I_n - P_A) %*% Xpcoa
  SSTm <- t(Xpcoa) %*% (I_n - P_T) %*% Xpcoa

  # Check that lengths are compatible to use with sign ponderation
  k <- ncol(Xpcoa)
  if (length(eig_sign) > k) eig_sign <- eig_sign[seq_len(k)]

  SSA <- sum(diag(SSAm) * eig_sign)
  SSR <- sum(diag(SSRm) * eig_sign)
  SST <- sum(diag(SSTm) * eig_sign)

  # --- Degrees of freedom ---
  dfA <- a - 1L
  dfR <- n - a
  dfT <- n - 1L

  # --- MS and pseudo-F ---
  MS_A <- SSA / dfA
  MS_R <- SSR / dfR
  F_A  <- MS_A / MS_R

  anova_tbl <- data.frame(
    Source   = c(nameA, "Residuals", "Total"),
    df       = c(dfA, dfR, dfT),
    SS       = c(SSA, SSR, SST),
    MS       = c(MS_A, MS_R, NA_real_),
    pseudoF  = c(F_A, NA_real_, NA_real_),
    check.names = FALSE
  )

  out_list <- list(
    SS  = c(SSA = SSA, SSR = SSR, SST = SST),
    df  = c(dfA = dfA, dfR = dfR, dfT = dfT),
    MS  = c(MS_A = MS_A, MS_R = MS_R),
    F   = c(F_A = F_A),
    method = method,
    transformation = transformation,
    model = model
  )

  switch(return,
         table = anova_tbl,
         list  = out_list,
         both  = list(table = anova_tbl, stats = out_list))
}


#' dbMANOVA para diseño anidado B(A) (tipo PERMANOVA de 2 vías)
#'
#' @param x        Matriz/data.frame de comunidad (muestras x especies). Debe ser el subset ya muestreado (balanceado).
#' @param factEnv  data.frame con las columnas de factores correspondientes a x:
#'                 - A: factor principal (p.ej. sector)
#'                 - B: factor anidado en A (p.ej. sitio)
#'                 Se aceptan nombres arbitrarios; se detectan por posición: primera col = A, segunda col = B.
#' @param transformation  Transformación para reducir dominancia. Una de: "none", "sqrt", "fourthroot", "Log (X+1)".
#'                 Por defecto "sqrt".
#' @param method   Método de disimilitud (pasado a vegan::vegdist), p.ej. "bray" (default), "euclidean", "jaccard", etc.
#' @param model    Etiqueta del modelo. Actualmente solo "nested.symmetric" (B anidado en A). Se valida.
#' @param return   "table", "list" o "both". Qué devolver (ANOVA-like table, lista de SS/df, o ambos).
#'
#' @return data.frame y/o list con SS, df, MS y pseudo-F:
#'         F_A = MS_A / MS_B(A), F_B(A) = MS_B(A) / MS_R.
#' @importFrom stats dist rchisq
#' @importFrom vegan vegdist
dbmanova_nested <- function(x,
                            factEnv,
                            transformation = c("square root","none","fourth root","Log (X+1)"),
                            method = "bray",
                            model = c("nested.symmetric"),
                            return = c("table","list","both")) {

  # --- Validaciones básicas ---
  transformation <- match.arg(transformation)
  method <- match.arg(method)
  model <- match.arg(model)
  return <- match.arg(return)

  if (!is.data.frame(factEnv)) factEnv <- as.data.frame(factEnv)
  if (!is.matrix(x) && !is.data.frame(x)) stop("`x` must be matrix/data.frame (comunidad).")

  Xcomm <- as.matrix(x)
  if (any(!is.finite(Xcomm))) stop("`x` contains non-finite values (NA/NaN/Inf).")

  n <- nrow(Xcomm)
  if (n != nrow(factEnv)) stop("`x` y `factEnv` must have the same number of rows (muestras).")

  # Coerce factores
  facA <- factor(factEnv[[1]], exclude = NULL)
  facB <- factor(factEnv[[2]], exclude = NULL)

  # --- Transformaciones comunes para atenuar dominancia ---
  transform_comm <- function(M, how = "square root") {
    how <- match.arg(how, c("none","square root","fourth root","Log (X+1)"))
    if (how == "none") return(M)
    if (how == "square root") return(sqrt(M))
    if (how == "fourthroot") return(M^(1/4))
    if (how == "Log (X+1)") return(log1p(M))
  }

  Xtr <- transform_comm(Xcomm, transformation)

  # --- Distance and PCoA (Gower) ---
  d <- vegan::vegdist(Xtr, method = method)
  Dm <- as.matrix(d)

  # Gower A Matrix from distances (classical MDS/PCoA)
  A <- -0.5 * (Dm^2)
  nn <- nrow(A)
  I <- diag(nn)
  J <- matrix(1/nn, nn, nn)
  G <- (I - J) %*% A %*% (I - J)

  EG <- eigen(G, symmetric = TRUE) 
  # Detecting the sign (>= 0 -> 1; negativo -> -1)
  sign_custom <- function(x) ifelse(x >= 0, 1, -1)

  eig_sign <- sign_custom(EG$values)

  # PCoA proyection 
  Lhalf <- sqrt(abs(EG$values))
  X <- sweep(EG$vectors, 2, Lhalf, `*`)   # PCoA coordinates

  # --- nested design: B(A) ---
  # H_A: n x a 
  H_A <- model.matrix(~ facA - 1)
  # facBA: nested combination B(A)
  facBA <- interaction(facA, facB, drop = TRUE)
  H_BA <- model.matrix(~ facBA - 1)

  a <- nlevels(facA)
  b_tot <- nlevels(facBA)

  # Proyections
  N_A  <- crossprod(H_A)              # t(H_A) %*% H_A
  P_A  <- H_A %*% solve(N_A) %*% t(H_A)

  N_BA <- crossprod(H_BA)
  P_BA <- H_BA %*% solve(N_BA) %*% t(H_BA)

  P_T  <- matrix(1/n, n, n)
  I_n  <- diag(n)

  # --- SSCP in the PCoA space ---
  SSA_mat  <- t(X) %*% (P_A  - P_T) %*% X
  SSBA_mat <- t(X) %*% (P_BA - P_A) %*% X
  SSR_mat  <- t(X) %*% (I_n  - P_BA) %*% X
  SST_mat  <- t(X) %*% (I_n  - P_T) %*% X

  # To move from SSCP -> scalar SS, we use the autovalues signs
  # Only the diagonal is used and the ponderation by eig_sign
  # It's important to ensure compatible lenghts by X eliminating null components
  k <- ncol(X)
  if (length(eig_sign) > k) eig_sign <- eig_sign[seq_len(k)]

  SSA <- sum(diag(SSA_mat)  * eig_sign)
  SSBA <- sum(diag(SSBA_mat) * eig_sign)
  SSR <- sum(diag(SSR_mat)  * eig_sign)
  SST <- sum(diag(SST_mat)  * eig_sign)

  # --- Degrees of freedom ---
  dfA  <- a - 1L
  dfBA <- b_tot - a
  dfR  <- n - b_tot
  dfT  <- n - 1L

  # --- Mean squares and pseudo-F ---
  MS_A  <- SSA  / dfA
  MS_BA <- SSBA / dfBA
  MS_R  <- SSR  / dfR

  F_A   <- MS_A  / MS_BA    
  F_BA  <- MS_BA / MS_R     

  anova_tbl <- data.frame(
    Source = c(names(factEnv)[1], paste0(names(factEnv)[2], "(", names(factEnv)[1], ")"), "Residuals", "Total"),
    df     = c(dfA, dfBA, dfR, dfT),
    SS     = c(SSA, SSBA, SSR, SST),
    MS     = c(MS_A, MS_BA, MS_R, NA_real_),
    `pseudoF` = c(F_A, F_BA, NA_real_, NA_real_),
    check.names = FALSE
  )

  out_list <- list(
    SS  = c(SSA = SSA, SSBA = SSBA, SSR = SSR, SST = SST),
    df  = c(dfA = dfA, dfBA = dfBA, dfR = dfR, dfT = dfT),
    MS  = c(MS_A = MS_A, MS_BA = MS_BA, MS_R = MS_R),
    F   = c(F_A = F_A, F_BA = F_BA),
    method = method,
    transformation = transformation,
    model = model
  )

  switch(return,
         table = anova_tbl,
         list  = out_list,
         both  = list(table = anova_tbl, stats = out_list))
}


