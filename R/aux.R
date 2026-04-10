# Estimation of Ecological Parameters of the Ecological Variables
#
# This function extracts the main parameters of the pilot data using base R functions

# ============================================================
# ecocbo2.0 - assempar_env(): estimación de parámetros de simulación
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(lme4)
})

# --------------------------
# Utilidades Beta
# --------------------------
beta_squeeze <- function(x, eps = 1e-4) pmin(pmax(x, eps), 1 - eps)

beta_moments_phi <- function(mu, var) {
  if (is.na(mu) || is.na(var) || var <= 0) {
    return(NA_real_)
  }
  phi <- mu * (1 - mu) / var - 1
  if (phi <= 0) NA_real_ else phi
}

# --------------------------
# Ajuste estático (sitio + muestra) y extracción
# --------------------------
fit_static_lmm_extract <- function(
  df,
  y,
  dist = c("lognormal", "normal"),
  site = "sitio"
) {
  dist <- match.arg(dist)

  y_raw <- df[[y]]
  if (dist == "lognormal") {
    if (any(y_raw <= 0, na.rm = TRUE)) {
      stop("Lognormal requiere y>0. Variable: ", y)
    }
    y_model <- log(y_raw)
    scale_name <- "log"
  } else {
    y_model <- y_raw
    scale_name <- "original"
  }

  d <- df %>%
    transmute(
      y = y_model,
      sitio = factor(.data[[site]])
    ) %>%
    filter(is.finite(y), !is.na(sitio))

  # Modelo: y ~ 1 + (1|sitio)
  mod <- lmer(y ~ 1 + (1 | sitio), data = d, REML = TRUE)

  vc <- as.data.frame(VarCorr(mod))
  sd_between <- vc %>%
    filter(grp == "sitio", var1 == "(Intercept)", is.na(var2)) %>%
    pull(sdcor)
  sd_within <- vc %>% filter(grp == "Residual") %>% pull(sdcor)

  tibble(
    scale = scale_name,
    mu = unname(fixef(mod)["(Intercept)"]),
    sd_between = ifelse(length(sd_between) == 0, NA_real_, sd_between),
    sd_within = ifelse(length(sd_within) == 0, NA_real_, sd_within),
    n = nrow(d),
    n_sites = n_distinct(d$sitio)
  )
}

fit_static_beta_extract <- function(
  df,
  y = "humedad",
  site = "sitio",
  min_val = 0,
  max_val = 100,
  eps = 1e-4
) {
  # Reescalar a (0,1)
  h01 <- (df[[y]] - min_val) / (max_val - min_val)
  h01 <- beta_squeeze(h01, eps = eps)

  d <- df %>%
    transmute(
      h01 = h01,
      sitio = factor(.data[[site]])
    ) %>%
    filter(is.finite(h01), !is.na(sitio))

  # Estimar mu y phi por momentos (global)
  mu01 <- mean(d$h01, na.rm = TRUE)
  var01 <- stats::var(d$h01, na.rm = TRUE)
  phi <- beta_moments_phi(mu01, var01)

  # Para separar tau y sigma_w, usamos aproximación logit-normal con random intercept
  logit <- function(p) log(p / (1 - p))
  d$z <- logit(d$h01)

  mod <- lmer(z ~ 1 + (1 | sitio), data = d, REML = TRUE)
  vc <- as.data.frame(VarCorr(mod))
  sd_between <- vc %>%
    filter(grp == "sitio", var1 == "(Intercept)", is.na(var2)) %>%
    pull(sdcor)
  sd_within <- vc %>% filter(grp == "Residual") %>% pull(sdcor)

  tibble(
    scale = "beta(0,1) + logit-normal RE",
    mu01 = mu01,
    phi = phi,
    sd_between_logit = ifelse(length(sd_between) == 0, NA_real_, sd_between),
    sd_within_logit = ifelse(length(sd_within) == 0, NA_real_, sd_within),
    n = nrow(d),
    n_sites = n_distinct(d$sitio),
    beta_min = min_val,
    beta_max = max_val,
    beta_eps = eps
  )
}

# --------------------------
# Función principal
# --------------------------
assempar_env <- function(
  data,
  site = "sitio",
  sample = "muestra",
  vars_lognormal = c("N", "P", "K", "conductividad"),
  vars_normal = c("temperatura", "pH"),
  var_beta = "humedad",
  humedad_range = c(0, 100),
  beta_eps = 1e-4
) {
  needed <- c(site, sample, vars_lognormal, vars_normal, var_beta)
  missing <- setdiff(needed, names(data))
  if (length(missing) > 0) {
    stop("Faltan columnas: ", paste(missing, collapse = ", "))
  }

  # (Opcional) verificación mínima: >=2 muestras por sitio para estimar residual con sentido
  chk <- data %>%
    count(.data[[site]], name = "n_muestras") %>%
    summarise(min_muestras = min(n_muestras), .groups = "drop")
  if (chk$min_muestras < 2) {
    warning("Hay sitios con <2 muestras; sd_within puede ser inestable.")
  }

  out_ln <- map_dfr(vars_lognormal, function(v) {
    fit_static_lmm_extract(data, y = v, dist = "lognormal", site = site) %>%
      mutate(var = v, dist = "lognormal") %>%
      relocate(var, dist)
  })

  out_n <- map_dfr(vars_normal, function(v) {
    fit_static_lmm_extract(data, y = v, dist = "normal", site = site) %>%
      mutate(var = v, dist = "normal") %>%
      relocate(var, dist)
  })

  out_b <- fit_static_beta_extract(
    data,
    y = var_beta,
    site = site,
    min_val = humedad_range[1],
    max_val = humedad_range[2],
    eps = beta_eps
  ) %>%
    mutate(var = var_beta, dist = "beta") %>%
    relocate(var, dist)

  bind_rows(out_ln, out_n, out_b) %>%
    arrange(var)
}

# --------------------------
# Uso típico
# --------------------------
# df <- read.csv("resultados.csv")
# params_static <- assempar_env_static(
#   df,
#   site = "mesocosmo",   # o "sitio"
#   sample = "lectura"    # o "muestra"
# )
# print(params_static)

# ============================================================
# simdata_env(): Simulación de datasets ambientales estáticos
# (sitio + muestra) usando parámetros estimados con assempar_env()
# ============================================================

simdata_env <- function(
  Par,
  cases,
  N,
  sites,
  jitter.base = 0,
  site_col = "site",
  sample_col = "sample",
  seed = NULL,
  return_params = FALSE
) {
  # ---- Validación básica ----
  if (!is.data.frame(Par)) {
    stop("Par debe ser un data.frame/tibble.")
  }
  req <- c("var", "dist")
  miss <- setdiff(req, names(Par))
  if (length(miss) > 0) {
    stop("Par no contiene columnas requeridas: ", paste(miss, collapse = ", "))
  }

  if (!is.numeric(cases) || cases < 1) {
    stop("cases debe ser >= 1.")
  }
  if (!is.numeric(N) || N < 1) {
    stop("N debe ser >= 1 (muestras por sitio).")
  }
  if (!is.numeric(sites) || sites < 1) {
    stop("sites debe ser >= 1 (sitios por dataset).")
  }
  if (!is.numeric(jitter.base) || jitter.base < 0) {
    stop("jitter.base debe ser >= 0.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Normalizar tipos
  Par <- Par %>%
    dplyr::mutate(
      var = as.character(.data$var),
      dist = tolower(as.character(.data$dist))
    )

  # ---- Helpers ----

  # Jitter multiplicativo aplicado a sd (como en SSP, pero aquí no hay fs/fw)
  jitter_sd <- function(sd, jitter.base) {
    if (is.na(sd) || sd == 0 || jitter.base == 0) {
      return(sd)
    }
    # sd' = sd * (1 + N(0, jitter.base)), truncado a >=0
    out <- sd * (1 + stats::rnorm(1, mean = 0, sd = jitter.base))
    pmax(out, 0)
  }

  inv_logit <- function(x) 1 / (1 + exp(-x))

  # Simulación por variable para un dataset
  simulate_variable <- function(p, sites, N) {
    dist <- p$dist[1]
    vname <- p$var[1]

    # Identificadores
    site_ids <- rep(seq_len(sites), each = N)
    sample_ids <- rep(seq_len(N), times = sites)

    if (dist %in% c("lognormal", "lnorm", "log-normal")) {
      # Requiere mu, sd_between, sd_within
      if (!all(c("mu", "sd_between", "sd_within") %in% names(p))) {
        stop("Par insuficiente para lognormal en variable: ", vname)
      }
      mu <- p$mu[1]
      sd_b <- jitter_sd(p$sd_between[1], jitter.base)
      sd_w <- jitter_sd(p$sd_within[1], jitter.base)

      b_i <- stats::rnorm(sites, mean = 0, sd = sd_b)
      z <- mu + b_i[site_ids] + stats::rnorm(sites * N, mean = 0, sd = sd_w)
      y <- exp(z)

      return(list(site = site_ids, sample = sample_ids, value = y, var = vname))
    }

    if (dist %in% c("normal", "gaussian")) {
      # Requiere mu, sd_between, sd_within
      if (!all(c("mu", "sd_between", "sd_within") %in% names(p))) {
        stop("Par insuficiente para normal en variable: ", vname)
      }
      mu <- p$mu[1]
      sd_b <- jitter_sd(p$sd_between[1], jitter.base)
      sd_w <- jitter_sd(p$sd_within[1], jitter.base)

      b_i <- stats::rnorm(sites, mean = 0, sd = sd_b)
      y <- mu + b_i[site_ids] + stats::rnorm(sites * N, mean = 0, sd = sd_w)

      return(list(site = site_ids, sample = sample_ids, value = y, var = vname))
    }

    if (dist %in% c("beta")) {
      # Requiere mu01 y phi y rango beta_min/beta_max
      # (en su salida: mu01 y phi; y además beta_min, beta_max, beta_eps)
      need <- c("mu01", "phi", "beta_min", "beta_max")
      if (!all(need %in% names(p))) {
        stop(
          "Par insuficiente para beta en variable: ",
          vname,
          ". Se requieren: ",
          paste(need, collapse = ", ")
        )
      }

      mu01 <- p$mu01[1]
      phi <- p$phi[1]
      if (is.na(mu01) || is.na(phi) || phi <= 0) {
        stop("Par beta inválidos para variable: ", vname, " (mu01/phi).")
      }

      # Parámetros alpha/beta
      alpha <- mu01 * phi
      beta <- (1 - mu01) * phi

      # Simular en (0,1)
      x <- stats::rbeta(sites * N, shape1 = alpha, shape2 = beta)

      # Reescalar a rango original (0-100 típicamente)
      minv <- p$beta_min[1]
      maxv <- p$beta_max[1]
      y <- minv + x * (maxv - minv)

      return(list(site = site_ids, sample = sample_ids, value = y, var = vname))
    }

    stop("Distribución no soportada: ", dist, " (variable: ", vname, ")")
  }

  # ---- Simular datasets ----
  vars <- unique(Par$var)

  out <- vector("list", length = cases)
  out_params <- vector("list", length = cases)

  for (cc in seq_len(cases)) {
    # Para reproducibilidad caso-a-caso si el usuario fija seed:
    # se mantiene la secuencia global; no re-seedeo aquí.

    # Simular cada variable y ensamblar en wide
    sim_list <- lapply(vars, function(v) {
      p <- Par[Par$var == v, , drop = FALSE]
      simulate_variable(p, sites = sites, N = N)
    })

    # Base de índices
    base <- data.frame(
      site = sim_list[[1]]$site,
      sample = sim_list[[1]]$sample
    )

    # Agregar variables
    for (obj in sim_list) {
      base[[obj$var]] <- obj$value
    }

    # Renombrar columnas id según argumentos
    names(base)[names(base) == "site"] <- site_col
    names(base)[names(base) == "sample"] <- sample_col

    out[[cc]] <- base

    if (return_params) out_params[[cc]] <- Par
  }

  if (return_params) {
    return(list(data = out, Par = out_params))
  }
  out
}
