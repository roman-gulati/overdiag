#' Plot incidence with gradual dissemination of screening in a population.
#'
#' Illustrate incidence projections produced as by
#' \code{multipopulation_setting}.
#'
#' @param dset A data frame as produced by \code{multipopulation_setting}.
#' @param minyear Minimum year to display projections.
#' @param maxyear Maximum year to display projections.
#' @param mininc Minimum incidence to display projections.
#' @param maxinc Maximum incidence to display projections.
#' @param line.size Width of lines to display projections.
#' @param minyear.text Minimum year to display text annotation.
#' @param maxyear.text Maximum year to display text annotation.
#' @param text.size Scale of text annotation.
#' @param text.offset Scale of gap between projections and text annotation.
#' @param text.angle Angle of text annotation.
#' @return A ggplot object.
#' @seealso \code{\link{trial_plot}}
#' @examples
#' library(plyr)
#' library(reshape)
#' library(ggplot2)
#' mpset <- multipopulation_setting()
#' \dontrun{
#' mp <- multipopulation_plot(mpset)
#' print(mp)
#' }
#' @export

multipopulation_plot <- function(dset,
                                 minyear=0,
                                 maxyear=20,
                                 mininc=0,
                                 maxinc=200,
                                 line.size=0.1,
                                 minyear.text=2,
                                 maxyear.text=maxyear-1,
                                 text.size=3,
                                 text.offset=12,
                                 text.angle=30){
    dset <- droplevels(subset(dset, minyear <= year & year <= maxyear))
    # extract background incidence
    baseline <- dset[with(dset, year == 0), 'count_clinical']
    # extract incidence under screening
    overdiag_inc <- with(dset, count_overdiag)
    screen_inc <- with(dset, count_screen)
    clinical_inc <- with(dset, count_clinical)
    # calculate empirical excess
    inc_total <- clinical_inc+screen_inc+overdiag_inc
    excess_total <- inc_total-baseline
    # specify valid years for evaluating empirical excess
    year_offset <- screen_inc > 0
    # calculate first year excess matches true overdiagnosis
    exact <- year_offset & abs(overdiag_inc/excess_total-1) < 0.000001
    exact_year <- min(which(exact))
    # append exact year
    dset <- transform(dset,
                      exact=year == year[exact_year],
                      exact_tag='Unbiased')
    dset.text <- subset(dset, year %in% seq(minyear.text, maxyear.text))
    dset.text <- transform(dset.text,
                           text_size=text.size,
                           text_offset=text.offset,
                           text_angle=text.angle)
    theme_set(theme_bw())
    theme_update(panel.grid.major=element_blank(),
                 panel.grid.minor=element_blank(),
                 panel.border=element_rect(colour='black', fill='transparent'),
                 axis.text=element_text(size=14),
                 axis.title=element_text(angle=0, size=18),
                 legend.position='none')
    gg <- ggplot(dset)
    gg <- gg+geom_blank(aes(y=maxinc))
    gg <- gg+geom_ribbon(aes(x=year,
                             ymin=0,
                             ymax=count_clinical),
                         fill='skyblue')
    gg <- gg+geom_ribbon(aes(x=year,
                             ymin=count_clinical,
                             ymax=count_clinical+count_screen),
                         fill='orange')
    gg <- gg+geom_ribbon(aes(x=year,
                             ymin=count_clinical+count_screen,
                             ymax=count_clinical+count_screen+count_overdiag),
                         fill='seagreen')
    gg <- gg+geom_line(aes(x=year, y=count_clinical), size=line.size)
    gg <- gg+geom_line(aes(x=year, y=count_clinical+count_screen), size=line.size)
    gg <- gg+geom_line(aes(x=year, y=count_clinical+count_screen+count_overdiag), size=line.size)
    gg <- gg+scale_x_continuous(name='\nYear',
                                breaks=seq(minyear, maxyear, by=2),
                                limits=c(minyear, maxyear),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name='Annual incidence per 100,000 individuals\n',
                                expand=c(0, 0))
    gg <- gg+geom_text(data=dset.text,
                       aes(x=year,
                           y=count_clinical+count_screen-text_offset,
                           label=paste0('+', round(count_screen)),
                           size=text_size,
                           angle=text_angle))
    gg <- gg+geom_text(data=dset.text,
                       aes(x=year,
                           y=count_clinical+count_screen+count_overdiag+text_offset,
                           label=paste0('+', round(count_overdiag)),
                           size=text_size,
                           angle=text_angle))
    if(nrow(subset(dset.text, count_clinical <= text_offset)) > 0)
        gg <- gg+geom_text(data=subset(dset.text, count_clinical <= text_offset),
                           aes(x=year,
                               y=count_clinical+text_offset,
                               label=paste0('-', 100-round(count_clinical)),
                               size=text_size,
                               angle=text_angle))
    gg <- gg+geom_text(data=subset(dset.text, count_clinical > text_offset),
                       aes(x=year,
                           y=count_clinical-text_offset,
                           label=paste0('-', 100-round(count_clinical)),
                           size=text_size,
                           angle=text_angle))
    gg <- gg+geom_segment(data=subset(dset, exact),
                          aes(x=year,
                              xend=year,
                              y=maxinc*0.85,
                              yend=maxinc*0.7),
                          arrow=arrow(length=unit(0.3, 'cm'), type='closed'))
    gg <- gg+geom_text(data=subset(dset, exact),
                       aes(x=year,
                           y=maxinc*0.9,
                           label=exact_tag),
                       vjust=0,
                       size=4)
    return(gg)
}

