#' Calculate clinical incidence in the presence of screening.
#'
#' Append to input data frame number of screen tests offered and number of
#' interval cancers not detected by screening.
#'
#' @details This function expects that the input data frame \code{dset} contains
#'   a variable \code{sojourn} with a single unique value. The number of false
#'   negatives are bounded by (1) the later of onset and the start of screening
#'   and (2) the earlier of clinical presentation and the end of screening.
#'   Cancers present clinically when all screening tests are false negatives.
#'   Imperfect attendance is represented by removing the expected number of
#'   cancers in individuals who attend a sensitive test scaled by the number of
#'   tests offered. In other words, false negative tests represent either not
#'   attending a test or attending a test but it does not detect latent cancer.
#' @param dset A data frame of cancer incidence as produced by
#'   \code{generate_absence}.
#' @param screen.start.year Year of follow-up at which screening starts.
#' @param screen.stop.year Year of follow-up at which screening stops.
#' @param attendance Proportion of individuals who attend screening tests.
#' @param sensitivity Proportion of relevant cancers detected by screening.
#' @return A data frame of simulated cancer incidence organized by year of
#'   preclinical onset, sojourn time, and year of clinical diagnosis.
#' @seealso \code{\link{calculate_screen}}
#' @examples
#' library(plyr)
#' library(reshape)
#' dset <- generate_absence(1000, 0.001, 0, 6, 10)
#' cset <- ddply(dset,
#'               .(sojourn),
#'               calculate_clinical,
#'               sensitivity=0.5,
#'               attendance=0.8,
#'               screen.start.year=0,
#'               screen.stop.year=10)
#' print(head(cset))
#' @export

calculate_clinical <- function(dset,
                               sensitivity,
                               attendance,
                               screen.start.year,
                               screen.stop.year){
    stopifnot(with(dset, length(unique(sojourn)) == 1))
    screen.start.year <- floor(screen.start.year)
    screen.stop.year <- floor(screen.stop.year)
    stopifnot(0 <= screen.start.year & screen.start.year <= screen.stop.year)
    stopifnot(0 <= attendance & attendance <= 1)
    stopifnot(0 <= sensitivity & sensitivity <= 1)
    dset <- transform(dset, lower_year=pmax(screen.start.year, onset_year))
    dset <- transform(dset, upper_year=pmin(screen.stop.year, clinical_year))
    dset <- transform(dset, tests_offered=pmax(upper_year-lower_year, 0))
    dset <- transform(dset, count_clinical=count_onset*(attendance*(1-sensitivity)+(1-attendance))^tests_offered)
    dset <- subset(dset, select=-c(lower_year, upper_year))
    return(dset)
}

