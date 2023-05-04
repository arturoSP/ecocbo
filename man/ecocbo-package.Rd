\name{ecocbo-package}

\alias{ecocbo-package}
\alias{ecocbo}

\docType{package}

\title{
ecocbo: Cost-benefit analysis in ecological sampling from simulated data
}

\description{
\pkg{ecocbo} is an R package that is designed to perform cost-benefit optimization for ecological sampling (U derwood & Chapman, 2003) based in the idea of analizing simulated ecological communities as a basis for this optimization (Guerra-Castro, et al., 2021).
}

\details{
The functions in \pkg{ecocbo} package can be used to identify the optimal number of sites and samples that must be considered in a community ecology study by using simulated data. Together with \pkg{SSP} package, \pkg{ecocbo} proposes a novel approach to the determination of he appropriate sampling effort in community ecology studies.

\pkg{ecocbo} is composed by four functions: \code{\link{sim_beta}} calculates statistical power for different sampling efforts and \code{\link{plot_power}} plots those results to help the user define the a combination of sampling effort and power to move on. \code{\link{scompvar}} calculates the components of variation for the analized dataset, and finally, \code{\link{sim_cbo}} determines an estimate of the number of sites and samples to consider to optimize the cost-benefit for an ecological sampling study.

\pkg{ecocbo} is being developed at Github(\url{https://github.com/edlinguerra/ecocbo}), where up-to-date versions can be found.
}

\references{
Underwood, A. J., Underwood, A. J., & Wnderwood, A. J. (1997). Experiments in ecology: their logical design and interpretation using analysis of variance. Cambridge university press.

Underwood, A. J., & Chapman, M. G. (2003). Power, precaution, Type II error and sampling design in assessment of environmental impacts. Journal of Experimental Marine Biology and Ecology, 296(1), 49-70.

Anderson, M. J. (2014). Permutational multivariate analysis of variance (PERMANOVA). Wiley statsref: statistics reference online, 1-15.

Guerra‐Castro, E. J., Cajas, J. C., Simões, N., Cruz‐Motta, J.J., & Mascaró, M. (2021). SSP: an R package to estimate sampling effort in studies of ecological communities. Ecography, 44(4), 561-573.
}

\author{
The \pkg{ecocbo} development team is Edlin Guerra-Castro and Arturo Sanchez-Porras.
}

\examples{

}

\keyword{ package }
