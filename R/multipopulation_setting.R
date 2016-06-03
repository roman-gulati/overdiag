#' Simulate incidence with gradual dissemination of screening in a population.
#'
#' Generate a data frame representing a population of \code{pop.size}
#' individuals and record relevant disease diagnosed with and without screening
#' and overdiagnoses under gradual dissemination of screening.
#'
#' @param pop.size Number of individuals in the simulated population.
#' @param onset.rate Annual incidence rate of relevant preclinical disease.
#' @param sojourn.min Shortest relevant preclinical duration.
#' @param sojourn.max Longest relevant preclinical duration.
#' @param sensitivity Screen test episode sensitivity.
#' @param overdiag.rate Proportion of screen detections that are overdiagnosed.
#' @param proportion Vector of proportions (that sum to 1) representing
#'   components of the population that initiate screening in specified years.
#' @param start.year Vector of years in which the corresponding component of
#'   the population initiates screening.
#' @param followup.years Number of years of follow-up.
#' @return A data frame of simulated disease incidence organized by year of
#'   preclinical onset, sojourn time, and year of diagnosis.
#' @seealso \code{\link{trial_setting}}
#' @examples
#' library(plyr)
#' library(reshape)
#' mpset <- multipopulation_setting()
#' print(head(mpset))
#' multipopulation_2ways <- function(dissemination.pattern='default'){
#'     if(dissemination.pattern == 'default'){
#'         proportion <- c(0.05, 0.1, 0.15, 0.15, 0.05, 0.5)
#'         start.year <- c(2, 3, 4, 5, 6, 28)
#'         dissemination.name <- 'Cumulative uptake (years 2-6): 5%,15%,30%,45%,50%'
#'     } else {
#'         proportion <- c(0.1, 0.3, 0.1, 0.5)
#'         start.year <- c(2, 3, 4, 28)
#'         dissemination.name <- 'Cumulative uptake (years 2-4): 10%,40%,50%'
#'     }
#'     mpset <- multipopulation_setting(proportion=proportion,
#'                                      start.year=start.year)
#'     mpset <- transform(mpset, dissemination=dissemination.name)
#'     return(mpset)
#' }
#' mpset_default <- multipopulation_2ways(dissemination.pattern='default')
#' mpset_variant <- multipopulation_2ways(dissemination.pattern='variant')
#' mpset_2ways <- rbind(mpset_default, mpset_variant)
#' print(head(mpset_2ways))
#' @export

multipopulation_setting <- function(pop.size=1e5,
                                    onset.rate=0.001,
                                    sojourn.min=0,
                                    sojourn.max=6,
                                    sensitivity=0.5,
                                    overdiag.rate=0.25,
                                    proportion=c(0.05, 0.1, 0.15, 0.15, 0.05, 0.5),
                                    start.year=c(2, 3, 4, 5, 6, 28),
                                    followup.years=30){
    dissemination <- data.frame(pop.size=pop.size,
                                proportion=proportion,
                                start.year=start.year)
    stopifnot(with(dissemination, sum(proportion)) == 1)
    stopifnot(with(dissemination, length(proportion) == length(start.year)))
    mpset <- ddply(dissemination,
                   .(proportion, start.year),
                   function(x){
                       with(x,
                            population_setting(pop.size=proportion*pop.size,
                                               onset.rate=onset.rate,
                                               sojourn.min=sojourn.min,
                                               sojourn.max=sojourn.max,
                                               sensitivity=sensitivity,
                                               overdiag.rate=overdiag.rate,
                                               screen.start.year=start.year,
                                               followup.years=followup.years))
                   })
    mpset <- ddply(mpset,
                   .(year),
                   summarize,
                   count_screen=sum(count_screen),
                   count_clinical=sum(count_clinical),
                   count_overdiag=sum(count_overdiag))
    return(mpset)
}

