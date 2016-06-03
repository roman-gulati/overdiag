#' Simulate incidence in the presence of screening.
#'
#' Screen a simulated population using specified sensitivity and count screen
#' diagnoses in each year of screening for relevant disease that develops in each
#' year with a given sojourn time.
#'
#' @param dset A data frame of simulated population as produced by
#'   \code{generate_absence}.
#' @param sojourn.min Minimum years of preclinical detectable period.
#' @param sojourn.max Maximum years of preclinical detectable period.
#' @param sensitivity Proportion of relevant disease detected by screening.
#' @param attendance Proportion of individuals who attend screening tests.
#' @param screen.start.year Year of follow-up at which screening starts.
#' @param screen.stop.year Year of follow-up at which screening stops.
#' @param followup.years Number of years of follow-up.
#' @return A data frame of simulated disease incidence organized by year
#'   of preclinical onset, sojourn time, and year of clinical diagnosis.
#' @seealso \code{\link{generate_absence}}, \code{\link{generate_overdiag}}
#' @examples
#' library(plyr)
#' library(reshape)
#' dset <- generate_absence(1000, 0.001, 0, 6, 10)
#' dset <- generate_presence(dset, 0, 6, 0.5, 0.8, 0, 10, 10)
#' print(head(dset))
#' @export

generate_presence <- function(dset,
                              sojourn.min,
                              sojourn.max,
                              sensitivity,
                              attendance,
                              screen.start.year,
                              screen.stop.year,
                              followup.years){
    followup.years <- floor(followup.years)
    screen.start.year <- floor(screen.start.year)
    screen.stop.year <- floor(screen.stop.year)
    stopifnot(0 <= screen.start.year & screen.start.year <= screen.stop.year)
    stopifnot(screen.stop.year <= followup.years)
    stopifnot(0 <= attendance & attendance <= 1)
    stopifnot(0 <= sensitivity & sensitivity <= 1)
    dset <- ddply(dset,
                  .(sojourn),
                  calculate_clinical,
                  screen.start.year=screen.start.year,
                  screen.stop.year=screen.stop.year,
                  attendance=attendance,
                  sensitivity=sensitivity)
    dset <- ddply(dset,
                  .(sojourn),
                  calculate_screen,
                  screen.start.year=screen.start.year,
                  screen.stop.year=screen.stop.year,
                  attendance=attendance,
                  sensitivity=sensitivity)
    dset <- subset(dset, select=-tests_offered)
    return(dset)
}

