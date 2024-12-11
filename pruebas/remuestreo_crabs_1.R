library(GAD)
library(tidyr)
library(dplyr)
library(stringr)

data("crabs")

Density <- crabs$Density
Region <- as.random(crabs$Region)
Location <- as.random(crabs$Location)
Site <- as.random(crabs$Site)

# aplicar esto para dos sitios (S1 + S2, S2 + S3 y S1 + S3)
modelo <- lm(Density ~ Region + Location %in% Region +
                Site %in% Location %in% Region)
estimates(modelo)
anova.pri <- gad(modelo)
anova(modelo) # nos da valores F de A/R siempre, pero GAD demuestra que
               # los denominadores son diferentes

CV <- comp.var(modelo, anova.pri) # el componente de variación es la columna Estimate
                          # el que queremos ver al remuestrear es el de
                          # Region:Location:site en sus tres valores para que
                          # se vaya centrando sobre 0.200


# escribir esto (5:16) dentro de un loop donde crabs se reescriba para distintos
# sites y quadrats (meter todo a balancedtwostage dentro del loop) y luego graficar
# estimates respecto a m y n (residual con n y sites con m) y luego al revés
# (residuals con m y sites con m) para ver el funnel y las horizontales, respectivamente.

# Repetir modelo1 para pares de Site


# Paquetes necesarios
library(GAD)
library(sampling)
library(ggplot2)
library(tidyr)
library(dplyr)

data("crabs")

# Parámetros de sampleo
site <- unique(crabs$Site)       # Sitios para muestrear
quadrat <- unique(crabs$Quadrat) # Cuadrantes para muestrear
n_site <- length(site)
n_quadrat <- length(quadrat)

# Dataframe de resultados
results <- data.frame(m = integer(), n = integer(), residuals = numeric(), estimate = numeric())

# loop para dar volumen a los datos
for(h in 1:50){
  for(i in 2:n_site){
    for(j in 2:n_quadrat){
      # se obtiene muestra de las etiquetas de sitios
      site_sampled <- sample(site, size = i)
      # se obtiene muestra de las etiquetas de cuadrantes
      quadrat_sampled <- sample(quadrat, size = j)

      # se extrae la muestra a partir de las etiquetas
      final_sample <- crabs[crabs$Site %in% site_sampled &
                              crabs$Quadrat %in% quadrat_sampled,]

      # Se guardan los datos como random, para usar funciones de GAD
      Density <- final_sample$Density
      Region <- as.fixed(final_sample$Region)
      Location <- as.random(final_sample$Location)
      Site <- as.random(final_sample$Site)

      modelo1 <- lm(Density ~ Region + Location %in% Region +
                      Site %in% Location %in% Region)

      # estimates(modelo1)
      anova1 <- gad(modelo1)
      # anova(modelo1)

      compvar <- comp.var(modelo1, anova1)

      results <- rbind(results, cbind(m = i,
                                      n = j,
                                      CVres = compvar$comp.var$Estimate[4],
                                      CVsite = compvar$comp.var$Estimate[3]))

    }
  }
}

# Graficar resultados
# Gráfico 1: Residuales respecto a n
ggplot(results, aes(x = n, y = CVres, color = factor(n))) +
  geom_point()+
  #geom_jitter(height=0)+
  scale_x_continuous(breaks = seq(2,5,1))+
  scale_y_continuous(breaks = seq(0.0, 1.8, 0.2))+
  labs(title = "Residuales vs n", x = "n (Cuadrantes)", y = "Residuales", color = "m (Sitios)") +
  theme_bw()

# Gráfico 2: Comp. Var. sitios respecto a m
ggplot(results, aes(x = m, y = CVsite, color = factor(m)))+
  #geom_jitter(height=0)+
  scale_x_continuous(breaks = seq(2,3,1))+
  scale_y_continuous(breaks = seq(-0.5, 0.8, 0.2))+
  geom_point()+
  labs(title = "Comp. Variation Site vs m", x = "m (Sitios)", y = "Comp. Variation Site", color = "n (Cuadrantes)")+
  theme_bw()

# Gráfico 3: Residuales respecto a m
ggplot(results, aes(x = m, y = CVres, color = factor(m))) +
  #geom_jitter(height=0)+
  scale_x_continuous(breaks = seq(2,3,1))+
  scale_y_continuous(breaks = seq(0.0, 1.8, 0.2))+
  geom_point()+
  labs(title = "Residuales vs m", x = "m (Sites)", y = "Residulas", color = "n (Cuadrantes)") +
  theme_bw()

# Gráfico 4: Estimaciones respecto a n y m
ggplot(results, aes(x = n, y = CVsite, color = factor(n)))+
  #geom_jitter(height=0)+
  geom_point()+
  scale_x_continuous(breaks = seq(2,5,1))+
  scale_y_continuous(breaks = seq(-0.5, 0.8, 0.2))+
  labs(title = "Comp. Var vs n", x = "n (Cuadrantes)", y = "Comp. Variation Site", color = "m (Sitios)")+
  theme_bw()

# Reduciendo la dimensionalidad en un nivel ----
data("crabs")

crabs |>
  mutate(sites = str_c(Location, Site, sep = "-")) |>
  ggplot(aes(x = sites, y = Crabs, color = sites)) +
  geom_jitter()

crabs2 <- crabs |>
  transmute(Region, Location, Site = paste0(Location, Site), Quadrat, Crabs, Density)

# Parámetros de sampleo
site <- unique(crabs2$Site)       # Sitios para muestrear
quadrat <- unique(crabs2$Quadrat) # Cuadrantes para muestrear
n_site <- length(site)
n_quadrat <- length(quadrat)

# Dataframe de resultados
results2 <- data.frame(
  m = integer(),
  n = integer(),
  residuals = numeric(),
  estimate = numeric()
)

# loop para dar volumen a los datos
for (h in 1:50) {
  for (i in 2:n_site) {
    for (j in 2:n_quadrat) {
      # se obtiene muestra de las etiquetas de sitios
      site_sampled <- sample(site, size = i)
      # se obtiene muestra de las etiquetas de cuadrantes
      quadrat_sampled <- sample(quadrat, size = j)

      # se extrae la muestra a partir de las etiquetas
      final_sample <- crabs2[crabs2$Site %in% site_sampled &
                               crabs2$Quadrat %in% quadrat_sampled, ]

      # Se guardan los datos como random, para usar funciones de GAD
      Density <- final_sample$Density
      Region <- as.fixed(final_sample$Region)
      Location <- as.random(final_sample$Location)
      Site <- as.random(final_sample$Site)

      modelo1 <- lm(Density ~ Region + Site %in% Region)
      print(paste("combinacion", i, j))

      # estimates(modelo1)
      anova1 <- gad(modelo1)
      # anova(modelo1)

      compvar <- comp.var(modelo1, anova1)

      results2 <- rbind(
        results2,
        cbind(
          m = i,
          n = j,
          CVres = compvar$comp.var$Estimate[3],
          CVsite = compvar$comp.var$Estimate[2]
        )
      )

    }
  }
}


# Graficar resultados
# Gráfico 1: Residuales respecto a n
ggplot(results2, aes(x = n, y = CVres, color = factor(n))) +
  geom_point() +
  #geom_jitter(height=0)+
  # scale_x_continuous(breaks = seq(2, 5, 1)) +
  scale_y_continuous(breaks = seq(0.0, 2.6, 0.2)) +
  labs(
    title = "Residuales vs n",
    x = "n (Cuadrantes)",
    y = "Residuales",
    color = "m (Sitios)"
  ) +
  theme_bw()

# Gráfico 2: Comp. Var. sitios respecto a m
ggplot(results2, aes(x = m, y = CVsite, color = factor(m))) +
  #geom_jitter(height=0)+
  # scale_x_continuous(breaks = seq(2, 10, 1)) +
  scale_y_continuous(breaks = seq(-0.5, 2.4, 0.2)) +
  geom_point() +
  labs(
    title = "Comp. Variation Site vs m",
    x = "m (Sitios)",
    y = "Comp. Variation Site",
    color = "n (Cuadrantes)"
  ) +
  theme_bw()

# Gráfico 3: Residuales respecto a m
ggplot(results2, aes(x = m, y = CVres, color = factor(m))) +
  #geom_jitter(height=0)+
  scale_y_continuous(breaks = seq(0.0, 2.6, 0.2)) +
  # scale_x_continuous(breaks = seq(2, 10, 1)) +
  geom_point() +
  labs(
    title = "Residuales vs m",
    x = "m (Sites)",
    y = "Residulas",
    color = "n (Cuadrantes)"
  ) +
  theme_bw()

# Gráfico 4: Estimaciones respecto a n y m
ggplot(results2, aes(x = n, y = CVsite, color = factor(n))) +
  #geom_jitter(height=0)+
  # scale_x_continuous(breaks = seq(2, 5, 1)) +
  scale_y_continuous(breaks = seq(-0.5, 2.4, 0.2)) +
  geom_point() +
  labs(
    title = "Comp. Var vs n",
    x = "n (Cuadrantes)",
    y = "Comp. Variation Site",
    color = "m (Sitios)"
  ) +
  theme_bw()
