#' Simulate cancer incidence in a population setting.
#'
#' Generate a data frame representing a population of \code{pop.size}
#' individuals and record relevant cancers diagnosed with and without screening
#' and overdiagnosed cancers.
#'
#' @param dset data frame containing \code{sojourn.min},
#'   \code{sojourn.max}, \code{sensitivity}, and \code{overdiag.rate}.
#' @param pop.size Number of individuals in the simulated population.
#' @param followup.years Number of years of follow-up.
#' @param onset.rate Number of relevant preclinical cancers that develop each year.
#' @param screen.start.year Year of follow-up at which screening starts.
#' @param screen.stop.year Year of follow-up at which screening stops.
#' @param verbose A logical flag to print setting information.
#' @return A data frame of simulated cancer incidence organized by year of
#'   preclinical onset, sojourn time, and year of diagnosis.
#' @seealso \code{\link{trial_setting}}
#' @examples
#' library(plyr)
#' library(reshape)
#' population_incidence_3ways <- function(pop.size=1e5,
#'                                        screen.start.year=3){
#'     dset <- data.frame(sojourn.min=c(0, 0),
#'                        sojourn.max=c(12, 6),
#'                        sensitivity=0.5,
#'                        overdiag.rate=0.25)
#'     pset <- ddply(dset,
#'                   .(sojourn.min,
#'                     sojourn.max,
#'                     sensitivity,
#'                     overdiag.rate),
#'                   population_setting,
#'                   pop.size=pop.size,
#'                   screen.start.year=screen.start.year)
#'     return(pset)
#' }
#' 
#' multipopulation_incidence_3ways <- function(pop.size=1e5,
#'                                             dissemination.pattern='default'){
#'     if(dissemination.pattern == 'default'){
#'         proportion <- c(0.05, 0.1, 0.15, 0.15, 0.05, 0.5)
#'         start.year <- c(2, 3, 4, 5, 6, 28)
#'         dissemination.name <- 'Cumulative uptake (years 2-6): 5%,15%,30%,45%,50%'
#'     } else {
#'         proportion <- c(0.1, 0.3, 0.1, 0.5)
#'         start.year <- c(2, 3, 4, 28)
#'         dissemination.name <- 'Cumulative uptake (years 2-4): 10%,40%,50%'
#'     }
#'     dissemination <- data.frame(pop.size=pop.size,
#'                                 proportion=proportion,
#'                                 start.year=start.year)
#'     stopifnot(with(dissemination, sum(proportion)) == 1)
#'     mpset <- ddply(dissemination,
#'                    .(proportion, start.year),
#'                    function(x){
#'                        cat('proportion:', with(x, proportion), '\n')
#'                        cat('start.year:', with(x, start.year), '\n')
#'                        with(x,
#'         population_incidence_3ways(pop.size=proportion*pop.size,
#'                                    screen.start.year=start.year))
#'                          })
#'     mpset <- ddply(mpset,
#'                    .(sojourn.min,
#'                      sojourn.max,
#'                      sensitivity,
#'                      overdiag.rate,
#'                      year),
#'                    summarize,
#'                    count_screen=sum(count_screen),
#'                    count_clinical=sum(count_clinical),
#'                    count_overdiag=sum(count_overdiag))
#'     mpset <- transform(mpset, dissemination=dissemination.name)
#'     return(mpset)
#' }
#' mpset_default <- multipopulation_incidence_3ways(dissemination.pattern='default')
#' mpset_variant <- multipopulation_incidence_3ways(dissemination.pattern='variant')
#' mpset_3ways <- rbind(mpset_default, mpset_variant)
#' @export

population_setting <- function(dset,
                               pop.size=1e5,
                               followup.years=30,
                               onset.rate=0.001,
                               screen.start.year=4,
                               screen.stop.year=30,
                               verbose=FALSE){
    if(verbose)
        with(dset, cat(paste(rep('-', 40), collapse=''),
                       '\nsensitivity:', as.character(unique(sensitivity)),
                       '\noverdiag.rate:', as.character(unique(overdiag.rate)),
                       '\n'))
    # generate population of pop.size individuals and
    # record year of clinical diagnosis for batches of relevant
    # cancers that develop in each year with a given sojourn time
    pset <- with(dset, generate_absence(pop.size,
                                        followup.years,
                                        onset.rate,
                                        sojourn.min,
                                        sojourn.max))
    # screen the population under assumed sensitivity by
    # counting screen diagnoses in each year of screening
    # for batches of relevant cancers that develop in each
    # year with a given sojourn time
    pset <- with(dset, generate_presence(pset,
                                         followup.years,
                                         screen.start.year,
                                         screen.stop.year,
                                         sojourn.min,
                                         sojourn.max,
                                         sensitivity,
                                         attendance=1))
    # append overdiagnoses as a constant fraction of screen
    # diagnoses in each year of screening
    pset <- with(dset, generate_overdiag(pset, overdiag.rate))
    # isolate clinical diagnoses by year and sojourn time
    clinical_sojourn <- ddply(pset,
                              .(clinical_year, sojourn),
                              summarize,
                              count_clinical=unique(count_clinical))
    # isolate screen diagnoses and overdiagnoses and sojourn time
    screened_sojourn <- ddply(pset,
                              .(screen_year, sojourn),
                              summarize,
                              count_screen=sum(count_screen),
                              count_overdiag=sum(count_overdiag))
    # count clinical diagnoses in each year
    clinical <- ddply(clinical_sojourn,
                      .(clinical_year),
                      summarize,
                      count_clinical=sum(count_clinical))
    # count screen diagnoses and overdiagnoses in each year
    screened <- ddply(screened_sojourn,
                      .(screen_year),
                      summarize,
                      count_screen=sum(count_screen),
                      count_overdiag=sum(count_overdiag))
    if(verbose)
        if(with(dset, sojourn.min != sojourn.max & overdiag.rate == 0)){
            leadtime <- with(dset,
                             sum(sapply(seq(sojourn.max-1),
                                        function(x)
                                            sensitivity*(1-sensitivity)^(sojourn.max-x-1)*x)))
            cat('sensitivity:', with(dset, sensitivity), '\n')
            cat('lead time:', round(leadtime, 2), '\n')
        }
    # pad screen diagnoses and overdiagnoses with 0
    # in all pre-screening years
    screened <- rbind(data.frame(screen_year=seq(0, screen.start.year-1),
                                 count_screen=0,
                                 count_overdiag=0),
                      screened)
    # merge screen diagnoses, overdiagnoses, and clinical diagnoses
    # from pre-screening through all years of screening
    merged <- merge(rename(screened, c('screen_year'='year')),
                    rename(clinical, c('clinical_year'='year')))
    return(merged)
}

