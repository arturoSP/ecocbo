library(devtools)
build()


library(SSP)
library(ecocbo)
library(dplyr)
library(tidyr)


# Importar datos playas proyecto PAPIIT2020
# Filtrar para datos de Sisal, 1era campaña

datos <- raw_dat|>
  filter(locality == "Sisal" & campaign == 1 & zone == "Wet")|>
  select(-c(locality, zone, campaign))|>
  relocate(site, .before = "labels")|>
  select(-labels)|>
  mutate(dummy = 1)

# Simulación de datos para H0 cierta, Ha cierta

#H0
datosH0 <- as.data.frame(datos)
datosH0$site <- as.factor("T0")

#Ha
datosHa <- as.data.frame(datos)
datosHa$site <- as.factor(datos$site)


# Parámetros para simular con assempar{SSP}

parH0 <- assempar(data = datosH0, type = "counts")
parHa <- assempar(data = datosHa, type = "counts")


#Simulaciones con simdata{SSP}

simH0Dat <- simdata(Par = parH0, cases = 10, N = 1000, sites = 1)
simHaDat <- simdata(Par = parHa, cases = 10, N = 100, sites = 10)

# Estimación de potencia con ecocbo

betaResult <- sim_beta(simH0 = simH0Dat,
                       simHa = simHaDat,
                       n = 10, m = 5, k = 10,
                       alpha = 0.05) #Agregar argumentos: transformation (ver sampsd), method, dummy (agrega especie con valores 1)

# Componentes de variación
cv <- scompvar(data = betaResult, n = 10, m = 10)
cv

plot_power(data = betaResult, n = 6, m = 4)
