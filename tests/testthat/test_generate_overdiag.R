context('Simulated incidence in the presence of screening and overdiagnosis')

pop.size <- 1000
onset.rate <- 0.001
sojourn.min <- 0
sojourn.max <- 6
followup.years <- 10

dset <- generate_absence(pop.size,
                         onset.rate,
                         sojourn.min,
                         sojourn.max,
                         followup.years)

sensitivity <- 0.5
attendance <- 0.8
screen.start.year <- 2
screen.stop.year <- 8

dset <- generate_presence(dset,
                          sojourn.min,
                          sojourn.max,
                          sensitivity,
                          attendance,
                          screen.start.year,
                          screen.stop.year,
                          followup.years)

overdiag.rate <- 0.25

dset <- generate_overdiag(dset, overdiag.rate)

test_that('Overdiagnosed cancers is correct fraction of screen cancers', {
          expect_equal(sum(dset$count_screen)*overdiag.rate/(1-overdiag.rate),
                       sum(dset$count_overdiag))
})

test_that('Years of onset have correct range', {
          onset.range <- seq(-sojourn.max, followup.years+sojourn.max)
          expect_equal(unique(dset$onset_year), onset.range)
})

test_that('Years of clinical diagnosis have correct range', {
          if(sojourn.min == sojourn.max)
              clinical.range <- seq(0, followup.years+2*sojourn.max)
          else
              clinical.range <- seq(-sojourn.max,
                                    followup.years+2*sojourn.max)
          expect_equal(unique(dset$clinical_year), clinical.range)
})

test_that('Data frame has correct number of rows', {
          sojourn.time.range <- seq(sojourn.min, sojourn.max)
          onset.range <- seq(-sojourn.max, followup.years+sojourn.max)
          screen.range <- seq(screen.start.year, screen.stop.year-1)
          expect_equal(nrow(dset), length(sojourn.time.range)*
                                   length(onset.range)*
                                   length(screen.range))
})

