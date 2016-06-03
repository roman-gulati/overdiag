#' Simulate cancer incidence in a trial setting.
#'
#' Generate a data frame representing a screening trial of \code{arm.size}
#' individuals in each arm and record relevant cancers diagnosed with and
#' without screening and overdiagnosed cancers.
#' @param arm.size Number of individuals in each trial arm.
#' @param onset.rate Annual incidence rate of relevant preclinical cancers.
#' @param sojourn.min Shortest relevant preclinical duration.
#' @param sojourn.max Longest relevant preclinical duration.
#' @param sensitivity Screen test episode sensitivity.
#' @param attendance Proportion of participants who attend a screen test.
#' @param overdiag.rate Proportion of screen-detected cancers that are
#'   overdiagnosed.
#' @param screen.start.year Year of follow-up at which screening starts.
#' @param screen.stop.year Year of follow-up at which screening stops.
#' @param followup.years Number of years of follow-up.
#' @return A data frame of simulated cancer incidence organized by year of
#'   preclinical onset, sojourn time, and year of diagnosis.
#' @seealso \code{\link{population_setting}}
#' @examples
#' library(plyr)
#' library(reshape)
#' trial_incidence_3types <- function(){
#'     tset <- data.frame(sojourn.min=0,
#'                        sojourn.max=6,
#'                        sensitivity=0.5,
#'                        attendance=0.8,
#'                        overdiag.rate=0.25)
#'     tset <- rbind(transform(tset, screen.start.year=1,
#'                                   screen.stop.year=5,
#'                                   type='screen-stop'),
#'                   transform(tset, screen.start.year=1,
#'                                   screen.stop.year=30,
#'                                   type='screen-continue'),
#'                   transform(tset, screen.start.year=5,
#'                                   screen.stop.year=30,
#'                                   type='screen-continue-delay'))
#'     tset <- ddply(tset,
#'                   .(sojourn.min,
#'                     sojourn.max,
#'                     sensitivity,
#'                     attendance,
#'                     overdiag.rate,
#'                     screen.start.year,
#'                     screen.stop.year,
#'                     type),
#'                   function(x)
#'                   with(x, trial_setting(sojourn.min=sojourn.min,
#'                                         sojourn.max=sojourn.max,
#'                                         sensitivity=sensitivity,
#'                                         attendance=attendance,
#'                                         overdiag.rate=overdiag.rate,
#'                                         screen.start.year=screen.start.year,
#'                                         screen.stop.year=screen.stop.year))
#'                   )
#'     tset <- transform(tset, count_clinical=ifelse(year == 0, 0, count_clinical))
#'     newscreen <- subset(tset, arm == 'screen' & type == 'screen-continue')
#'     newcontrol <- subset(tset, arm == 'screen' & type == 'screen-continue-delay')
#'     newscreen <- transform(newscreen, type='screen-continue-control-start')
#'     newcontrol <- transform(newcontrol, arm='control', type='screen-continue-control-start')
#'     tset <- subset(tset, type != 'screen-continue-delay')
#'     tset <- rbind(tset, newscreen, newcontrol)
#'     return(tset)
#' }
#' tset_3types <- trial_incidence_3types()
#' @export

trial_setting <- function(arm.size=50000,
                          onset.rate=0.001,
                          sojourn.min=0,
                          sojourn.max=6,
                          sensitivity=0.5,
                          attendance=0.8,
                          overdiag.rate=0.25,
                          screen.start.year=1,
                          screen.stop.year=30,
                          followup.years=30){
    # generate trial population of arm.size individuals and
    # record year of clinical diagnosis for batches of relevant
    # cancers that develop in each year with a given sojourn time
    # to serve as a basis for either arm
    tset <- generate_absence(arm.size,
                             onset.rate,
                             sojourn.min,
                             sojourn.max,
                             followup.years)
    # construct the control arm by "screening" the trial population
    # under 0 sensitivity and counting screen diagnoses in each year
    # of screening for batches of relevant cancers that develop in each
    # year with a given sojourn time
    cset <- generate_presence(tset,
                              sojourn.min,
                              sojourn.max,
                              0,
                              0,
                              screen.start.year,
                              screen.stop.year,
                              followup.years)
    # append overdiagnoses as a constant fraction of screen
    # diagnoses in each year of screening
    cset <- generate_overdiag(cset, overdiag.rate)
    # count control arm screen diagnoses in each year of screening
    control_screen <- ddply(cset,
                           .(screen_year),
                           summarize,
                           arm='control',
                           count_screen=sum(count_screen),
                           count_overdiag=sum(count_overdiag))
    # count clinical diagnoses in each year after first collapsing
    # across screening years
    control_clinical <- ddply(cset,
                             .(clinical_year, sojourn),
                             summarize,
                             count_clinical=unique(count_clinical))
    control_clinical <- ddply(control_clinical,
                             .(clinical_year),
                             summarize,
                             arm='control',
                             count_clinical=sum(count_clinical))
    # merge control arm screen and clinical diagnoses
    control_arm <- merge(rename(control_clinical, c('clinical_year'='year')),
                         rename(control_screen, c('screen_year'='year')),
                         all=TRUE)
    # coerce missing control arm screen diagnoses and overdiagnoses to 0
    control_arm[is.na(control_arm)] <- 0
    # process the trial population under assumed sensitivity
    # to construct the screen arm
    # construct the screen arm by screening the trial population under
    # given sensitivity and counting screen diagnoses in each year of
    # screening for batches of relevant cancers that develop in each
    # year with a given sojourn time
    sset <- generate_presence(tset,
                              sojourn.min,
                              sojourn.max,
                              sensitivity,
                              attendance,
                              screen.start.year,
                              screen.stop.year,
                              followup.years)
    # append overdiagnoses as a constant fraction of screen
    # diagnoses in each year of screening
    sset <- generate_overdiag(sset, overdiag.rate)
    # count screen arm screen diagnoses in each year of screening
    screen_screen <- ddply(sset,
                           .(screen_year),
                           summarize,
                           arm='screen',
                           count_screen=sum(count_screen),
                           count_overdiag=sum(count_overdiag))
    # count clinical diagnoses in each year after first collapsing
    # across screening years
    screen_clinical <- ddply(sset,
                             .(clinical_year, sojourn),
                             summarize,
                             count_clinical=unique(count_clinical))
    screen_clinical <- ddply(screen_clinical,
                             .(clinical_year),
                             summarize,
                             arm='screen',
                             count_clinical=sum(count_clinical))
    # merge screen arm screen and clinical diagnoses
    screen_arm <- merge(rename(screen_clinical, c('clinical_year'='year')),
                        rename(screen_screen, c('screen_year'='year')),
                        all=TRUE)
    # coerce missing screen arm screen diagnoses and overdiagnoses to 0
    screen_arm[is.na(screen_arm)] <- 0
    # merge screen and control arms
    merged <- rbind(control_arm, screen_arm)
    # restrict to given years of follow-up
    merged <- subset(merged, 0 <= year & year <= followup.years)
    return(merged)
}

