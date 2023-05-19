\name{simH0Dat}
\alias{simH0Dat}
\docType{data}
\title{
Dataset for using in the example of `sim_beta`
}
\description{
This dataset contains half of the necessary data to test \code{sim_beta}, the other half being \code{\link{simHaDat}}. These datasets are the result of applying the preparation steps mentioned in the vignette "ecocbo-guide": (1) subset `epibionts`, (2) adjust the `site` labels, (3) apply \code{assempar}, and (4) apply \code{simdata}.
}
\usage{data("simH0Dat")}
\format{
  A list formed by 3 data frames that contain 1000 rows (observations) and 95 columns (93 simulated species and 2 labels) each. The last column \code{site} is a factor with one level, to indicate how all observations belong to the same site.
}
\source{
Data available from the Dryad Digital Repository: <http://dx.doi.org/10.5061/dryad.3bk3j9kj5> (Guerra-Castro et al. 2020).
}
\references{
Guerra-Castro, E. J. et al. 2016. Scales of spatial variation in tropical benthic assemblages and their ecological relevance: epibionts on Caribbean mangrove roots as a model system. – Mar. Ecol. Prog. Ser. 548: 97–110.
}
\examples{
data(simH0Dat)
\dontrun{
sim_beta(simH0Dat, simHaDat, n = 10, m = 3, k = 20, alpha = 0.05)
}
}
\keyword{datasets}
