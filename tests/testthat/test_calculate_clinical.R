context('Calculated clinical incidence in the presence of screening')

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

dset <- ddply(dset,
              .(sojourn),
              calculate_clinical,
              sensitivity=sensitivity,
              attendance=attendance,
              screen.start.year=screen.start.year,
              screen.stop.year=screen.stop.year)

test_that('Number of tests offered is sane', {
          expect_true(max(dset$tests_offered) <= sojourn.max)
          expect_true(max(dset$tests_offered) <= followup.years)
})

