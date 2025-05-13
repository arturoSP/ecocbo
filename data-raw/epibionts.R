epibionts <- readr::read_csv("./data-raw/epibionts.csv")[,-1]

# epiDat ----
# After checking for the sites that are less similar it was found that E-2, M-1 and I-3 will be used.
 epiDat <- epibionts[epibionts$sector=="E"&epibionts$site==2 |
             epibionts$sector=="M"&epibionts$site==1 |
             epibionts$sector=="I"&epibionts$site==3,]
 epiDat[,"site"] <- paste0(epiDat$sector, epiDat$site)
 epiDat <- epiDat[,-1]

usethis::use_data(epiDat, overwrite = TRUE)

# epiBetaR ----
# Apply `ecocbo::sim_beta()` to epiDat and save an iteration so that examples
# that use the result of said function will load faster.
## Declaración de parámetros ----
cases = 5
N = 20
M = 15
n = 12
m = 12
k = 20
transformation = "square root"
method = "bray"
dummy = F
useParallel = T
model = "single.factor"
type = "counts"
Sest.method = "average"
data <- epiDat

## Función prep_data ----
epiBetaR <- prep_data(data, type, Sest.method, cases, N, M, n, m, k,
                           transformation, method, dummy, useParallel, model)

epiBetaR <- sim_beta(epiBetaR, 0.05)
usethis::use_data(epiBetaR, overwrite=TRUE)
