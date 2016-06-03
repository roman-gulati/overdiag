#' Simulate incidence of overdiagnosis.
#'
#' Append to input data frame the number of overdiagnosed cases given the
#' specified frequency of overdiagnosis.
#'
#' @param dset data frame of simulated population as produced by
#'   \code{generate_absence} after processing by \code{generate_presence}.
#' @param overdiag.rate Proportion of screen-detected cases that are
#'   overdiagnosed.
#' @return A data frame of simulated disease incidence organized by year of
#'   preclinical onset, sojourn time, and year of clinical diagnosis.
#' @seealso \code{\link{generate_absence}}, \code{\link{generate_presence}}
#' @examples
#' library(plyr)
#' library(reshape)
#' dset <- generate_absence(1000, 0.001, 0, 6, 10)
#' dset <- generate_presence(dset, 0, 6, 0.5, 0.8, 0, 10, 10)
#' dset <- generate_overdiag(dset, 0.25)
#' print(head(dset))
#' @export

generate_overdiag <- function(dset, overdiag.rate){
    stopifnot('count_screen' %in% names(dset))
    stopifnot(0 <= overdiag.rate & overdiag.rate <= 1)
    dset <- transform(dset, count_overdiag=overdiag.rate*count_screen/(1-overdiag.rate))
    return(dset)
}

