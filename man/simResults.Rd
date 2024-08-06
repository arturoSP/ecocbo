\name{simResults}
\alias{simResults}
\docType{data}
\title{
Data set containing the results of applying ecocbo::prep_data().
}
\description{
The dataset contains the results of applying ecocbo::prep_data() to epiDat. The result is a list with one level: $Results is a data frame with the results of applying PERMANOVA to epiDat a number of times, it contains the values of pseudoF and the mean squares for different repeated sampling efforts.

This dataset can be used to study the variability of the pseudoF-statistic, beta and the power when an experiment is applied to a varying number of samples, sampling units, or sampling sites.
}
\usage{data("simResults")}
\format{
An object of class "ecocbo_data", also a list containing one data frame. The format is:
  \tabular{rll}{
  $Results \tab             \tab a data frame that contains the results of the evaluation of sim_beta. \cr
           \tab dat.sim     \tab simulation from which the results are obtained. \cr
           \tab k           \tab number of resample for the result. \cr
           \tab m           \tab number of sites considered for the result. \cr
           \tab n           \tab number of replicates within each site for the result. \cr
           \tab pseudoFH0   \tab observed F value for the experimental design, when all observations belong to one site. \cr
           \tab pseudoFHa   \tab observed F value for the experimental design, when observations belong to different sites. \cr
           \tab AMSHa       \tab calculated mean squares among sites in the experiment. \cr
           \tab RMSHa       \tab calculated mean squares for the residuals in the experiment. \cr
  }
}
\details{
This dataset comes from applying ecocbo::prep_data() to the basic data from ecocbo::epiDat.
}
\source{
Data available from the Dryad Digital Repository: <http://dx.doi.org/10.5061/dryad.3bk3j9kj5> (Guerra-Castro et al. 2020).
}
\references{
Guerra-Castro, E. J. et al. 2016. Scales of spatial variation in tropical benthic assemblages and their ecological relevance: epibionts on Caribbean mangrove roots as a model system. – Mar. Ecol. Prog. Ser. 548: 97–110.
}
\examples{
data(simResults)

sim_beta(simResults, alpha = 0.05)
}
\keyword{datasets}
