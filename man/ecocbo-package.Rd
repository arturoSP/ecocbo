\name{ecocbo-package}

\alias{ecocbo-package}
\alias{ecocbo}

\docType{package}

\title{
ecocbo: Calculating Optimum Sampling Effort in Community Ecology
}

\description{
 A system for calculating the optimal sampling effort, based on the ideas of "Ecological cost-benefit optimization" as developed by A. Underwood (1997, ISBN 0 521 55696 1). Data is obtained from simulated ecological communities, and the optimization follows the following procedure of four functions (1) sim_beta() estimates statistical power and type 2 error by using Permutational Multivariate Analysis of Variance, (2) plot_power() represents the results of the previous function, (3) scompvar() calculates the variation components necessary for (4) sim_cbo() to calculate the optimal combination of number of sites and samples depending on either an economical budget or on a desired statistical accuracy.
}

\details{
The functions in \pkg{ecocbo} package can be used to identify the optimal number of sites and samples that must be considered in a community ecology study by using simulated data. Together with \pkg{SSP} package, \pkg{ecocbo} proposes a novel approach to the determination of he appropriate sampling effort in community ecology studies.

\pkg{ecocbo} is composed by four functions: \code{\link{sim_beta}} calculates statistical power for different sampling efforts and \code{\link{plot_power}} plots those results to help the user define the a combination of sampling effort and power to move on. \code{\link{scompvar}} calculates the components of variation for the analized dataset, and finally, \code{\link{sim_cbo}} determines an estimate of the number of sites and samples to consider to optimize the cost-benefit for an ecological sampling study.

\pkg{ecocbo} is being developed at Github(\url{https://github.com/arturoSP/ecocbo}), where up-to-date versions can be found.
}

\references{
Underwood, A. J. (1997). Experiments in ecology: their logical design and interpretation using analysis of variance. Cambridge university press.

Underwood, A. J., & Chapman, M. G. (2003). Power, precaution, Type II error and sampling design in assessment of environmental impacts. Journal of Experimental Marine Biology and Ecology, 296(1), 49-70.

Anderson, M. J. (2014). Permutational multivariate analysis of variance (PERMANOVA). Wiley statsref: statistics reference online, 1-15.

Guerra‐Castro, E. J., Cajas, J. C., Simões, N., Cruz‐Motta, J.J., & Mascaró, M. (2021). SSP: an R package to estimate sampling effort in studies of ecological communities. Ecography, 44(4), 561-573.
}

\author{
The \pkg{ecocbo} development team is Edlin Guerra-Castro and Arturo Sanchez-Porras.
}

\examples{
library(SSP)

\donttest{
# Load and adjust data.
data(epiDat)

epiH0 <- epiDat
epiH0[,"site"] <- as.factor("T0")
epiHa <- epiDat
epiHa[,"site"] <- as.factor(epiHa[,"site"])

# Calculate simulation parameters.
parH0 <- SSP::assempar(data = epiH0, type = "counts", Sest.method = "average")
parHa <- SSP::assempar(data = epiHa, type = "counts", Sest.method = "average")

# Simulation.
simH0Dat <- SSP::simdata(parH0, cases = 3, N = 1000, sites = 1)
simHaDat <- SSP::simdata(parHa, cases = 3, N = 100, sites = 10)

# Determination of power
epiBetaR <- sim_beta(simH0Dat, simHaDat,
                     n = 5, m = 4, k = 30,
                     alpha = 0.05,
                     transformation = "square root", method = "bray",
                     dummy = FALSE,
                     useParallel = FALSE)
}
epiBetaR

# Visualization of power
plot_power(data = epiBetaR, n = NULL, m = 3, method = "both")

# Computing components of variation
compVar <- scompvar(data = epiBetaR)
compVar

# Cost-benefit optimization
cboResult <- sim_cbo(comp.var = compVar, ct = 20000, ck = 100, cj = 2500)
cboResult

}

\keyword{ package }
