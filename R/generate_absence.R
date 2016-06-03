#' Simulate cancer incidence in the absence of screening.
#'
#' Generate a data frame representing a population of \code{pop.size}
#' individuals and record year of clinical diagnosis for relevant cancers that
#' develop in each year with a given sojourn time distribution under specified
#' follow-up.
#'
#' @param pop.size Number of individuals in the simulated population.
#' @param onset.rate Individual rate at which relevant preclinical cancers
#'   develop each year.
#' @param sojourn.min Minimum years of preclinical detectable period.
#' @param sojourn.max Maximum years of preclinical detectable period.
#' @param followup.years Number of years of follow-up.
#' @return A data frame of simulated cancer incidence organized by year of
#'   preclinical onset, number of preclinical cancers, sojourn time, and year of
#'   clinical diagnosis.
#' @seealso \code{\link{generate_presence}}, \code{\link{generate_overdiag}}
#' @examples
#' library(plyr)
#' library(reshape)
#' dset <- generate_absence(1000, 0.001, 0, 6, 10)
#' print(head(dset))
#' @export

generate_absence <- function(pop.size,
                             onset.rate,
                             sojourn.min,
                             sojourn.max,
                             followup.years){
    pop.size <- floor(pop.size)
    followup.years <- floor(followup.years)
    stopifnot(pop.size > 0)
    stopifnot(followup.years > 0)
    stopifnot(onset.rate >= 0)
    stopifnot(0 <= sojourn.min & sojourn.min <= sojourn.max)
    if(sojourn.min == sojourn.max){
        sojourn.denom <- 1
        sojourn.time <- sojourn.min
    } else {
        sojourn.denom <- sojourn.max-sojourn.min+1
        sojourn.time <- seq(sojourn.min, sojourn.max)
    }
    probabilities <- rep(1/sojourn.denom, sojourn.denom)
    stopifnot(sum(probabilities) == 1)
    pset <- data.frame(onset_year=seq(-sojourn.max, followup.years),
                       count_onset=pop.size*onset.rate)
    if(sojourn.max > 0)
        pset <- rbind(pset,
                      data.frame(onset_year=followup.years+seq(sojourn.max),
                                 count_onset=0))
    pset <- ddply(pset,
                  .(onset_year),
                  function(x)
                      with(x,
                           data.frame(onset_year=unique(onset_year),
                                      count_onset=probabilities*unique(count_onset),
                                      sojourn=sojourn.time)))
    pset <- transform(pset, clinical_year=onset_year+sojourn)
    return(pset)
}

