#' Plot incidence in a randomized trial.
#'
#' Illustrate incidence projections produced as by \code{trial_setting}.
#'
#' @param dset A data frame as produced by \code{multipopulation_setting}.
#' @param minyear Minimum year to display projections.
#' @param maxyear Maximum year to display projections.
#' @return A ggplot object.
#' @seealso \code{\link{multipopulation_plot}}
#' @examples
#' library(plyr)
#' library(reshape)
#' library(ggplot2)
#' tset <- trial_setting()
#' \dontrun{
#' tp <- trial_plot(tset)
#' print(tp)
#' }
#' @export

trial_plot <- function(dset, minyear=0, maxyear=12){
    dset <- droplevels(subset(dset, minyear <= year & year <= maxyear))
    dset <- ddply(dset,
                  .(arm),
                  transform,
                  count_nonoverdiag=count_clinical+count_screen,
                  count_total=count_clinical+count_screen+count_overdiag)
    dset <- transform(dset, measure='Annual')
    cset <- transform(dset, measure='Cumulative')
    cset <- ddply(cset,
                  .(arm),
                  transform,
                  count_clinical=cumsum(count_clinical),
                  count_screen=cumsum(count_screen),
                  count_overdiag=cumsum(count_overdiag),
                  count_nonoverdiag=cumsum(count_nonoverdiag),
                  count_total=cumsum(count_total))
    dset <- rbind(dset, cset)
    dset <- ddply(dset,
                  .(measure),
                  function(x){
                      # extract years for each arm
                      screen_year <- with(x, year[arm == 'screen'])
                      control_year <- with(x, year[arm == 'control'])
                      # extract screen diagnoses for each arm
                      screen_screen <- with(x, count_screen[arm == 'screen'])
                      control_screen <- with(x, count_screen[arm == 'control'])
                      # extract overdiagnoses for each arm
                      screen_overdiag <- with(x, count_overdiag[arm == 'screen'])
                      control_overdiag <- with(x, count_overdiag[arm == 'control'])
                      # extract total diagnoses for each arm
                      screen_total <- with(x, count_total[arm == 'screen'])
                      control_total <- with(x, count_total[arm == 'control'])
                      # calculate empirical excess and true overdiagnosis
                      excess_overdiag <- screen_overdiag-control_overdiag
                      excess_total <- screen_total-control_total
                      # specify valid years for evaluating empirical excess
                      if(any(control_screen > 0))
                          year_offset <- FALSE # control_screen > 0
                      else
                          year_offset <- screen_screen > 0
                      # calculate first year excess matches true overdiagnosis
                      exact <- year_offset & abs(excess_overdiag/excess_total-1) < 0.000001
                      exact_year <- min(which(exact))
                      stopifnot(isTRUE(all.equal(screen_year[exact_year],
                                                 control_year[exact_year])))
                      # append exact years to stratum
                      x <- transform(x,
                         exact=arm == 'screen' & year == screen_year[exact_year],
                         exact_tag='Unbiased')
                      return(x)
                  })
    maxanninc <- with(subset(dset, measure == 'Annual'), max(count_total))*1.34
    maxcuminc <- with(subset(dset, measure == 'Cumulative'), max(count_total))*1.34
    theme_set(theme_bw())
    theme_update(panel.grid.major=element_blank(),
                 panel.grid.minor=element_blank(),
                 panel.border=element_rect(colour='black', fill='transparent'),
                 axis.text=element_text(size=14),
                 axis.title=element_text(angle=0, size=18),
                 strip.background=element_rect(fill=NA, colour=NA),
                 strip.text=element_text(angle=0, size=12),
                 legend.position='none')
    gg <- ggplot(dset)
    # highlight overdiagnosis in the screen arm
    {
        dset.screen <- subset(dset, arm == 'screen')
        gg <- gg+geom_ribbon(data=dset.screen,
                             aes(x=year,
                                 ymin=count_nonoverdiag,
                                 ymax=count_total),
                             size=0.5,
                             fill=alpha('darkgreen', 0.4),
                             colour=alpha('darkgreen', 0.4))
    }
    # highlight overdiagnosis in the control arm
    {
        dset.control <- subset(dset, arm == 'control')
        gg <- gg+geom_ribbon(data=dset.control,
                             aes(x=year,
                                 ymin=count_nonoverdiag,
                                 ymax=count_total),
                             size=0.5,
                             fill=alpha('purple', 0.4),
                             colour=alpha('purple', 0.4))
    }
    gg <- gg+geom_line(aes(x=year, y=count_total, colour=arm),
                       alpha=0.75,
                       size=1)
    gg <- gg+geom_point(aes(x=year,
                            y=count_total,
                            group=arm,
                            colour=arm),
                        alpha=1,
                        size=2,
                        shape=19)
    gg <- gg+facet_grid(measure~., scales='free_y')
    gg <- gg+scale_x_continuous(name='\nYears of follow-up',
                                breaks=seq(minyear, maxyear, by=2),
                                limits=c(minyear, maxyear),
                                expand=c(0, 0.5))
    gg <- gg+scale_y_continuous(name='Number of cases\n',
                                expand=c(0.05, 0))
    gg <- gg+geom_blank(aes(y=0))
    gg <- gg+geom_blank(data=subset(dset, measure == 'Annual'),
                        aes(y=maxanninc))
    gg <- gg+geom_blank(data=subset(dset, measure == 'Cumulative'),
                        aes(y=maxcuminc))
    gg <- gg+scale_colour_manual(values=c('control'=alpha('purple', 0.4),
                                          'screen'=alpha('darkgreen', 0.4)))
    gg <- gg+scale_fill_manual(values=c('control'=alpha('purple', 0.4),
                                        'screen'=alpha('darkgreen', 0.4)))
    {# highlight screen events in each arm
        sset <- subset(dset, arm == 'screen' &
                             count_screen > 0 &
                             measure == 'Annual')
        if(nrow(sset) > 0)
            gg <- gg+geom_segment(data=sset,
                                  aes(x=year,
                                      xend=year,
                                      y=-Inf,
                                      yend=0),
                                  size=1,
                                  colour=alpha('darkgreen', 0.4))
        cset <- subset(dset, arm == 'control' &
                             count_screen > 0 &
                             measure == 'Annual')
        if(nrow(cset) > 0)
            gg <- gg+geom_segment(data=cset,
                                  aes(x=year,
                                      xend=year,
                                      y=0,
                                      yend=10),
                                  size=1,
                                  colour=alpha('purple', 0.4))
    }
    if(nrow(subset(dset, exact) > 0)){
        maxinc <- ifelse(subset(dset, exact)[['measure']] == 'Annual',
                         maxanninc,
                         maxcuminc)
        gg <- gg+geom_vline(data=subset(dset, exact),
                            aes(xintercept=year),
                            colour='red')
        gg <- gg+geom_text(data=subset(dset, exact),
                           aes(x=year,
                               y=maxinc,
                               label=exact_tag),
                           colour='red',
                           angle=90,
                           hjust=1.1,
                           vjust=1.5,
                           size=5)
    }
    return(gg)
}

