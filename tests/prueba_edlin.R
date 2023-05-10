library(devtools)
load_all()
build()


library(SSP)
library(ecocbo)
library(dplyr)
library(tidyr)

raw_dat <- readxl::read_xls("./data-raw/dat.Abundancias.Playas.PAPIIT.xls")

# Importar datos playas proyecto PAPIIT2020
# Filtrar para datos de Sisal, 1era campaña

datos <- raw_dat|>
  filter(locality == "Sisal" & campaign == 1 & zone == "Wet")|>
  select(-c(locality, zone, campaign))|>
  relocate(site, .before = "labels")|>
  select(-labels)

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
                       n = 15, m = 5, k = 20,
                       alpha = 0.05,
                       transformation = "fourth root",
                       method = "bray",
                       dummy = TRUE)

profvis::profvis(sim_beta(simH0 = simH0Dat,
                       simHa = simHaDat,
                       n = 10, m = 3, k = 20,
                       alpha = 0.05,
                       transformation = "square root")) #Agregar argumentos: transformation (ver sampsd), method, dummy (agrega especie con valores 1)

# Componentes de variación
cv <- scompvar(data = betaResult)
cv

plot_power(data = betaResult, m = 5)
plot_power(data = betaResult2, m = 5)
