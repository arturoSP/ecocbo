library(GAD)

data("crabs")

Density <- crabs$Density
Region <- as.random(crabs$Region)
Location <- as.random(crabs$Location)
Site <- as.random(crabs$Site)

# aplicar esto para dos sitios (S1 + S2, S2 + S3 y S1 + S3)
modelo1 <- lm(Density ~ Region + Location %in% Region +
                Site %in% Location %in% Region)
estimates(modelo1)
anova1 <- gad(modelo1)
anova(modelo1) # nos da valores F de A/R siempre, pero GAD demuestra que
               # los denominadores son diferentes

comp.var(modelo1, anova1) # el componente de variación es la columna Estimate
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
      Region <- as.random(final_sample$Region)
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
                                      residuals = compvar$comp.var[4,2],
                                      estimate = compvar$comp.var[3,2]))

    }
  }
}

# Graficar resultados
# Gráfico 1: Residuales respecto a n y m
ggplot(results, aes(x = n, y = residuals, color = factor(n))) +
  geom_jitter(height=0)+
  labs(title = "Residuales vs n y m", x = "n (Cuadrantes)", y = "Residuales", color = "m (Sitios)") +
  theme_minimal()

# Gráfico 2: Residuales respecto a m y n
ggplot(results, aes(x = m, y = residuals, color = factor(m)))+
  geom_point(height=0)+
  labs(title = "Residuales vs m y n", x = "m (Sitios)", y = "Residuales", color = "n (Cuadrantes)")+
  theme_minimal()

# Gráfico 3: Estimaciones respecto a m y n
ggplot(results, aes(x = m, y = estimate, color = factor(m))) +
  geom_jitter(height=0)+
  # geom_point()+
  labs(title = "Estimaciones vs m y n", x = "m (Sitios)", y = "Estimación", color = "n (Cuadrantes)") +
  theme_minimal()

# Gráfico 4: Estimaciones respecto a n y m
ggplot(results, aes(x = n, y = estimate, color = factor(n)))+
  geom_point(height=0)+
  labs(title = "Estimaciones vs n y m", x = "n (Cuadrantes)", y = "Estimación", color = "m (Sitios)")+
  theme_minimal()

##################
# Reduciendo la dimensionalidad en un nivel ----
##################
crabs2 <- crabs |>
  transmute(Region, Location,
            Site = paste0(Site,Quadrat),
            Crabs, Density)

# Parámetros de sampleo
site <- unique(crabs2$Location)       # Sitios para muestrear
quadrat <- unique(crabs2$Site) # Cuadrantes para muestrear
n_site <- length(site)
n_quadrat <- length(quadrat)

# Dataframe de resultados
results <- data.frame(m = integer(), n = integer(), residuals = numeric(), estimate = numeric())

# loop para dar volumen a los datos
for(h in 1:50){
  for(i in 2:n_site){
    for(j in 3:n_quadrat){
      # se obtiene muestra de las etiquetas de sitios
      site_sampled <- sample(site, size = i)
      # se obtiene muestra de las etiquetas de cuadrantes
      quadrat_sampled <- sample(quadrat, size = j)

      # se extrae la muestra a partir de las etiquetas
      final_sample <- crabs2[crabs2$Location %in% site_sampled &
                              crabs2$Site %in% quadrat_sampled,]

      # Se guardan los datos como random, para usar funciones de GAD
      Density <- final_sample$Density
      Region <- as.random(final_sample$Region)
      Location <- as.random(final_sample$Location)
      Site <- as.random(final_sample$Site)

      modelo1 <- lm(Density ~ Region + Location %in% Region)

      # estimates(modelo1)
      anova1 <- gad(modelo1)
      # anova(modelo1)

      compvar <- comp.var(modelo1, anova1)

      results <- rbind(results, cbind(m = i,
                                      n = j,
                                      residuals = compvar$comp.var[3,2],
                                      estimate = compvar$comp.var[2,2]))

    }
  }
}


# Graficar resultados
# Gráfico 1: Residuales respecto a n y m
ggplot(results, aes(x = n, y = residuals, color = factor(m))) +
  # geom_jitter(height=0)+
  geom_point()+
  labs(title = "Residuales vs n y m", x = "n (Cuadrantes)", y = "Residuales", color = "m (Sitios)") +
  theme_minimal()

# Gráfico 2: Residuales respecto a m y n
ggplot(results, aes(x = m, y = residuals, color = factor(n)))+
  #geom_jitter(height=0)+
  geom_point()+
  labs(title = "Residuales vs m y n", x = "m (Sitios)", y = "Residuales", color = "n (Cuadrantes)")+
  theme_minimal()

# Gráfico 3: Estimaciones respecto a m y n
ggplot(results, aes(x = m, y = estimate, color = factor(n))) +
  #geom_jitter(height=0)+
  geom_point()+
  labs(title = "Estimaciones vs m y n", x = "m (Sitios)", y = "Estimación", color = "n (Cuadrantes)") +
  theme_minimal()

# Gráfico 4: Estimaciones respecto a n y m
ggplot(results, aes(x = n, y = estimate, color = factor(m)))+
  geom_jitter(height=0)+
  labs(title = "Estimaciones vs n y m", x = "n (Cuadrantes)", y = "Estimación", color = "m (Sitios)")+
  theme_minimal()


# probando con epibionts ----

epibionts <- SSP::epibionts

# Parámetros de sampleo
site <- unique(epibionts$sector)       # Sitios para muestrear
quadrat <- as.factor(unique(epibionts$site)) # Cuadrantes para muestrear
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
      final_sample <- epibionts[epibionts$sector %in% site_sampled &
                                epibionts$site %in% quadrat_sampled,]

      # Se guardan los datos como random, para usar funciones de GAD
      Density <- final_sample$Amathia.sp
      Sector <- as.fixed(final_sample$sector)
      Site <- as.random(final_sample$site)

      modelo1 <- lm(Density ~ Sector + Site %in% Sector)

      # estimates(modelo1)
      anova1 <- gad(modelo1)
      # anova(modelo1)

      compvar <- comp.var(modelo1, anova1)

      results <- rbind(results, cbind(m = i,
                                      n = j,
                                      residuals = compvar$comp.var[3,2],
                                      estimate = compvar$comp.var[2,2]))

    }
  }
}

# Graficar resultados
# Gráfico 1: Residuales respecto a n y m
ggplot(results, aes(x = n, y = residuals, color = factor(n))) +
  # geom_jitter(height=0)+
  geom_point()+
  scale_y_continuous(breaks= seq(0,120,20))+
  labs(title = "Residuales vs n y m", x = "n (Cuadrantes)", y = "Residuales", color = "m (Sitios)") +
  theme_minimal()

# Gráfico 2: Residuales respecto a m y n
ggplot(results, aes(x = m, y = residuals, color = factor(m)))+
  # geom_jitter(height=0)+
  geom_point()+
  scale_y_continuous(breaks= seq(0,120,20))+
  labs(title = "Residuales vs m y n", x = "m (Sitios)", y = "Residuales", color = "n (Cuadrantes)")+
  theme_minimal()

# Gráfico 3: Estimaciones respecto a m y n
ggplot(results, aes(x = m, y = estimate, color = factor(m))) +
  # geom_jitter(height=0)+
  geom_point()+
  scale_y_continuous(breaks=seq(-14,14,2))+
  labs(title = "Estimaciones vs m y n", x = "m (Sitios)", y = "Estimación", color = "n (Cuadrantes)") +
  theme_minimal()

# Gráfico 4: Estimaciones respecto a n y m
ggplot(results, aes(x = n, y = estimate, color = factor(n)))+
  # geom_jitter(height=0)+
  geom_point()+
  scale_y_continuous(breaks=seq(-14,14,2))+
  labs(title = "Estimaciones vs n y m", x = "n (Cuadrantes)", y = "Estimación", color = "m (Sitios)")+
  theme_minimal()
