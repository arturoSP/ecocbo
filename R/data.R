#' Dataset on species count of marine communities.
#'
#' This is a dataset containing a subset from the epibionts dataset from `SSP`
#' which was made by using the three local communities that differ the most.
#'
#' @format A data frame with count of individuals for 24 observations on 151 species.
#'
#' @source Data available from the Dryad Digital Repository:
#' <http://dx.doi.org/10.5061/dryad.3bk3j9kj5> (Guerra-Castro et al. 2020).
#'

"epiDat"

#' Data set containing the results of applying ecocbo::prep_data().
#'
#' The dataset contains the results of applying ecocbo::prep_data() to epiDat.
#' The result is a list with one level: $Results is a data frame with the results
#' of applying PERMANOVA to epiDat a number of times, it contains the values of
#' pseudoF and the mean squares for different repeated sampling efforts.
#'
#' This dataset can be used to study the variability of the pseudoF-statistic,
#' beta and the power when an experiment is applied to a varying number of samples,
#' sampling units, or sampling sites.
#'
#' @format An object of class "ecocbo_data", also a list containing one data frame.
#' The format is:
#' \describe{
#'   \item{$Results}{
#'     \describe{
#'       \item{dat.sim}{simulation from which the results are obtained.}
#'       \item{k}{number of resample for the result.}
#'       \item{n}{number of replicates within each site for the result.}
#'       \item{pseudoFH0}{observed F value for the experimental design, when all observations belong to one site.}
#'       \item{pseudoFHa}{observed F value for the experimental design, when observations belong to different sites.}
#'       \item{MSR}{calculated mean squares for the residuals in the experiment.}
#'   }}
#'   \item{$model}{"single.factor"}
#'   \item{attribute}{class: ecocbo_data}
#' }
#'
#' @source Data available from the Dryad Digital Repository:
#' <http://dx.doi.org/10.5061/dryad.3bk3j9kj5> (Guerra-Castro et al. 2020).
#'

"simResults"

#' Data set containing the results of applying ecocbo::sim_beta() to a single factor
#' experiment.
#'
#' The dataset contains the results of applying ecocbo::sim_beta() to an excerpt
#' from the dataset epibionts from the package SSP. The result is a list with 4
#' components.
#'
#' This dataset can be used to study the variability of the pseudoF-statistic,
#' beta and the power when an experiment is applied to a varying number of
#' samples, sampling units, or sampling sites.
#'
#' @format An object of class "ecocbo_beta", also a list containing four components.
#' The format is:
#' \describe{
#'   \item{$Power}{
#'   \describe{
#'     \item{m}{number of sites considered for the result.}
#'     \item{n}{number of replicates within each site for the result.}
#'     \item{Power}{estimated statistical power.}
#'     \item{Beta}{estimated type II error.}
#'     \item{fCrit}{estimated pseudoF value that corresponds to the 1-alpha quartile of the distribution of pseudoF.}
#'   }
#'   }
#'   \item{$Results}{
#'     \describe{
#'       \item{dat.sim}{simulation from which the results are obtained.}
#'       \item{k}{number of resample for the result.}
#'       \item{m}{number of sites considered for the result.}
#'       \item{n}{number of replicates within each site for the result.}
#'       \item{pseudoFH0}{observed F value for the experimental design, when all observations belong to one site.}
#'       \item{pseduoFHa}{observed F value for the experimental design, when observations belong to different sites.}
#'       \item{MSB(A)}{calculated mean squares among sites in the experiment.}
#'       \item{MSR}{calculated mean squares for the residuals in the experiment.}
#'     }
#'   }
#'   \item{$alpha}{usually 0.05}
#'   \item{$model}{nested.symmetric}
#'   \item{attribute}{ecocbo.beta}
#' }
#'
#' @source Data available from the GitHub Digital Repository:
#' <https://github.com/edlinguerra/SSP/tree/master/data> (Guerra-Castro et al.
#' 2022).
#'

"epiBetaR"

#' Dataset on species count of coastal macrofauna.
#'
#' This is a dataset containing a subset from the macrofauna recorded in the
#' PAPIIT experiment.
#'
#' @format A dataframe with counts of individuals for 43 observations on 34 species.
#'
#' @source Data available from the GitHub Digital Repository:
#' <https://github.com/edlinguerra/IA206320_publico/tree/main/datos> (Guerra-Castro
#' et al. 2022).
#'

"macrofDat"

#' Data set containing the results of applying ecocbo::prep_data() to a nested
#' factors experiment.
#'
#' The dataset contains the results of applying ecocbo::prep_data() to epiDat.
#' The result is a list with one level: $Results is a data frame with the results
#' of applying PERMANOVA to epiDat a number of times, it contains the values of
#' pseudoF and the mean squares for different repeated sampling efforts.
#'
#' This dataset can be used to study the variability of the pseudoF-statistic,
#' beta and the power when an experiment is applied to a varying number of samples,
#' sampling units, or sampling sites.
#'
#' @format An object of class "ecocbo_data", also a list containing one data frame.
#' The format is:
#' \describe{
#'   \item{$Results}{
#'     \describe{
#'       \item{dat.sim}{simulation from which the results are obtained.}
#'       \item{k}{number of resample for the result.}
#'       \item{m}{number of sites considered for the result.}
#'       \item{n}{number of replicates within each site for the result.}
#'       \item{pseudoFH0}{observed F value for the experimental design, when all observations belong to one site.}
#'       \item{pseudoFHa}{observed F value for the experimental design, when observations belong to different sites.}
#'       \item{MSB(A)}{calculated mean squares among sites in the experiment.}
#'       \item{MSR}{calculated mean squares for the residuals in the experiment.}
#'   }}
#'   \item{$model}{"single.factor"}
#'   \item{attribute}{class: ecocbo_data}
#' }
#'
#' @source Source data is available from
#' <https://github.com/edlinguerra/IA206320_publico/tree/main/datos>
#' (Guerra-Castro et al. 2020).
#'

"simResultsNested"


#' Data set containing the results of applying ecocbo::sim_beta() to a nested
#' factors experiment.
#'
#' The dataset contains the results of applying ecocbo::sim_beta() to the dataset
#' from PAPIIT experiment. The result is a list with 4 components.
#'
#' This dataset can be used to study the variability of the pseudoF-statistic,
#' beta and the power when an experiment is applied to a varying number of
#' samples, sampling units, or sampling sites.
#'
#' @format An object of class "ecocbo_beta", also a list containing four components.
#' The format is:
#' \describe{
#'   \item{$Power}{
#'   \describe{
#'     \item{m}{number of sites considered for the result.}
#'     \item{n}{number of replicates within each site for the result.}
#'     \item{Power}{estimated statistical power.}
#'     \item{Beta}{estimated type II error.}
#'     \item{fCrit}{estimated pseudoF value that corresponds to the 1-alpha quartile of the distribution of pseudoF.}
#'     }}
#'   \item{$Results}{
#'     \describe{
#'       \item{dat.sim}{simulation from which the results are obtained.}
#'       \item{k}{number of resample for the result.}
#'       \item{m}{number of sites considered for the result.}
#'       \item{n}{number of replicates within each site for the result.}
#'       \item{pseudoFH0}{observed F value for the experimental design, when all observations belong to one site.}
#'       \item{pseduoFHa}{observed F value for the experimental design, when observations belong to different sites.}
#'       \item{MSB(A)}{calculated mean squares among sites in the experiment.}
#'       \item{MSR}{calculated mean squares for the residuals in the experiment.}
#'     }}
#'   \item{$alpha}{usually 0.05}
#'   \item{$model}{"nested.symmetric"}
#'   \item{attribute}{"ecocbo.beta"}
#' }
#'
#' @source Data available from the GitHub Digital Repository:
#' <https://github.com/edlinguerra/IA206320_publico/tree/main/datos> (Guerra-Castro
#' et al. 2022).
#'

"betaNested"


