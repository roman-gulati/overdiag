#' Calculate screen incidence in the presence of screening.
#'
#' Expand input data frame to indicate year of diagnosis and number of cancers
#' detected by screening.
#'
#' @details This function expects that the input data frame \code{dset} contains
#'   a variable \code{sojourn} with a single unique value. Create new variables
#'   to record year of each screening round (exclusive of stopping year). The
#'   number of false negatives are bounded (1) below by 0 and screen round - 1
#'   and (2) above by number of tests - 1 and sojourn time - 1 Check that
#'   incidence of onset matches diagnoses. Reshape dataset to indicate number of
#'   screen detections in each screen year.
#' @param dset A data frame of cancer incidence as produced by
#'   \code{generate_absence} after processing by \code{calculate_clinical}.
#' @param screen.start.year Year of follow-up at which screening starts.
#' @param screen.stop.year Year of follow-up at which screening stops.
#' @param attendance Proportion of individuals who attend screening tests.
#' @param sensitivity Proportion of relevant cancers detected by screening.
#' @return A data frame of simulated cancer incidence organized by year of
#'   preclinical onset, sojourn time, and year of clinical diagnosis.
#' @seealso \code{\link{calculate_clinical}}
#' @examples
#' library(plyr)
#' library(reshape)
#' dset <- generate_absence(1000, 10, 0.001, 0, 6)
#' cset <- ddply(dset,
#'               .(sojourn),
#'               calculate_clinical,
#'               screen.start.year=0,
#'               screen.stop.year=10,
#'               attendance=0.8,
#'               sensitivity=0.5)
#' sset <- ddply(cset,
#'               .(sojourn),
#'               calculate_screen,
#'               screen.start.year=0,
#'               screen.stop.year=10,
#'               attendance=0.8,
#'               sensitivity=0.5)
#' print(head(sset))
#' @export

calculate_screen <- function(dset,
                             screen.start.year,
                             screen.stop.year,
                             attendance,
                             sensitivity){
    stopifnot(with(dset, length(unique(sojourn)) == 1))
    screen.start.year <- floor(screen.start.year)
    screen.stop.year <- floor(screen.stop.year)
    stopifnot(0 <= screen.start.year & screen.start.year <= screen.stop.year)
    stopifnot(0 <= attendance & attendance <= 1)
    stopifnot(0 <= sensitivity & sensitivity <= 1)
    sojourn_time <- with(dset, unique(sojourn))
    for(screen_year in seq(screen.start.year, screen.stop.year-1))
        dset[[paste('count_screen', screen_year, sep='_')]] <- 0
    if(sojourn_time > 0){
        for(screen_year in seq(screen.start.year, screen.stop.year-1)){
            prev_years <- seq(screen_year-sojourn_time+1, screen_year)
            latent_index <- with(dset, which(onset_year %in% prev_years))
            latent_dset <- dset[latent_index, ]
            latent_times <- with(latent_dset, seq(sojourn_time-1, 0))
            latent_tests_offered <- with(latent_dset,
                                         pmin(latent_times, tests_offered-1))
            latent_tests_offered_and_valid <- with(latent_dset,
                                           pmin(latent_tests_offered,
                                                screen_year-screen.start.year))
            latent_tests_offered_and_valid <- with(latent_dset,
                                       pmax(latent_tests_offered_and_valid, 0))
            prob_onetp <- with(latent_dset,
                               sensitivity*attendance*(attendance*(1-sensitivity)+(1-attendance))^latent_tests_offered_and_valid)
            cases <- with(latent_dset, prob_onetp*count_onset)
            screen_index <- paste('count_screen', screen_year, sep='_')
            dset[latent_index, screen_index] <- cases
        }
    }
    count_all <- grep('count_screen|count_clinical', names(dset), value=TRUE)
    diag_counts <- rowSums(dset[, count_all])
    onset_counts <- dset[, 'count_onset']
    #index <- which(diag_counts != onset_counts)
    #if(any(diag_counts != onset_counts))
    #    browser()
    stopifnot(isTRUE(all.equal(diag_counts, onset_counts)))
    count_screen <- grep('count_screen', names(dset), value=TRUE)
    melted <- melt(dset, measure.vars=count_screen)
    melted <- rename(melted, c('variable'='screen_year', 'value'='count_screen'))
    melted <- transform(melted, screen_year=as.integer(sub('count_screen_', '', screen_year)))
    return(melted)
}

