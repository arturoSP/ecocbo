% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{permanova_oneway}
\alias{permanova_oneway}
\title{PERMANOVA one-way}
\usage{
permanova_oneway(
  x,
  factEnv,
  type = "P",
  method = "bray",
  transformation = "none"
)
}
\arguments{
\item{x}{ecological community data.}

\item{factEnv}{label for the community data.}

\item{type}{which algorithm to use for the calculation? At the moment, the only
option is "P".}

\item{method}{appropriate distance/dissimilarity metric (e.g. Gower,
Bray–Curtis, Jaccard, etc).}

\item{transformation}{Mathematical function to reduce the weight of very
dominant species.}
}
\value{
A data frame containing the resulting PERMANOVA table.
}
\description{
Calculates observed F and mean squares for the residuals and among sites. This
function is a helper for \code{\link[=prep_data]{prep_data()}}.
}
\references{
Underwood, A. J. (1997). Experiments in ecology: their logical
design and interpretation using analysis of variance. Cambridge university
press.

Anderson, M. J. (2014). Permutational multivariate analysis of
variance (PERMANOVA). Wiley statsref: statistics reference online, 1-15.
}
\seealso{
\code{\link[vegan:vegdist]{vegan::vegdist()}}
}
\author{
Edlin Guerra-Castro (\email{edlinguerra@gmail.com}), Arturo Sanchez-Porras
}
\keyword{internal}
