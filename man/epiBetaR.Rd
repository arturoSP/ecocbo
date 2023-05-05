\name{epiBetaR}
\alias{epiBetaR}
\docType{data}
\title{
Dataset containing an example of the results of applying `ecocbo::sim_beta`.
}
\description{
The dataset contains the results of applying ecocbo::sim_beta to epiDat. The result is a list with two levels. $Power A data frame with the results of statistical power according to different sampling efforts. $Results A data frame with the results of applying PERMANOVA to epiDat a number of times, it contains the values of pseudoF and the mean squares for different repeated sampling efforts. This dataset can be used to study the variability of the pseudoF-statistic, beta and the power when an experiment is applied to a varying number of samples, sampling units, or sampling sites.
}
\usage{data("epiBetaR")}
\format{
  An object of class "ecocbo_beta", also a list of two data frames. The format is:
  \describe{
  \item{code{$Power}}{a data frame with 12 obs. of 7 variables}
    \item{code{m}}{a numeric vector}
    \item{code{n}}{a numeric vector}
    \item{code{Power}}{a numeric vector}
    \item{code{Beta}}{a numeric vector}
    \item{code{fCrit}}{a numeric vector}
    \item{code{AMSHa}}{a numeric vector}
    \item{code{RMSHa}}{a numeric vector}
  \item{code{$Results}}{a data frame with 720 obs. of 8 variables}
    \item{code{dat.sim}}{a numeric vector}
    \item{code{k}}{a numeric vector}
    \item{code{m}}{a numeric vector}
    \item{code{n}}{a numeric vector}
    \item{code{pseudoFH0}}{a numeric vector}
    \item{code{pseudoFHa}}{a numeric vector}
    \item{code{AMSHa}}{a numeric vector}
    \item{code{RMSHa}}{a numeric vector}
}
}
\details{
This dataset comes from applying `ecocbo::sim_beta()` to the basic data from `ecocbo::epiDat`
}
\source{
Data available from the Dryad Digital Repository: <http://dx.doi.org/10.5061/dryad.3bk3j9kj5> (Guerra-Castro et al. 2020).
}
\references{
Guerra-Castro, E. J. et al. 2016. Scales of spatial variation in tropical benthic assemblages and their ecological relevance: epibionts on Caribbean mangrove roots as a model system. – Mar. Ecol. Prog. Ser. 548: 97–110.
}
\examples{
data(epiBetaR)
plot_power(epiBetaR, n = 3, m = 4, method = "both")
}
\keyword{datasets}
