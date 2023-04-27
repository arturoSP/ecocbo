epibionts <- read_csv("./data-raw/epibionts.csv")[,-1]

# epiDat ----
# After checking for the sites that are less similar it was found that E-2, M-1 and I-3 will be used.
 epiDat <- epibionts[epibionts$sector=="E"&epibionts$site==2 |
             epibionts$sector=="M"&epibionts$site==1 |
             epibionts$sector=="I"&epibionts$site==3,]
 epiDat[,"site"] <- paste0(epiDat$sector, epiDat$site)
 epiDat <- epiDat[,-1]

usethis::use_data(epiDat, overwrite = TRUE)

# epiBetaR ----
# Apply `ecocbo::sim_beta()` to epiDat and save an iteration so that examples that use the result
# of said function will load faster.

epiH0 <- epiDat
epiH0[,"site"] <- as.factor("T0")
epiHa <- epiDat
epiHa[,"site"] <- as.factor(epiHa[,"site"])

# Calculate simulation parameters.
parH0 <- SSP::assempar(data = epiH0, type = "counts", Sest.method = "average")
parHa <- SSP::assempar(data = epiHa, type = "counts", Sest.method = "average")

# Simulation.
simH0 <- SSP::simdata(parH0, cases = 3, N = 1000, sites = 1)
simHa <- SSP::simdata(parHa, cases = 3, N = 100, sites = 10)

epiBetaR <- sim_beta(simH0, simHa, n = 5, m = 4, k = 20, 0.05)
