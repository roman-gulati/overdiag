context('Calculated clinical incidence in the presence of screening')

pop.size <- 1000
followup.years <- 10
onset.rate <- 0.001
sojourn.min <- 0
sojourn.max <- 6

dset <- generate_absence(pop.size,
                         followup.years,
                         onset.rate,
                         sojourn.min,
                         sojourn.max)

screen.start.year <- 2
screen.stop.year <- 8
attendance <- 0.8
sensitivity <- 0.5

dset <- ddply(dset,
              .(sojourn),
              calculate_clinical,
              screen.start.year=screen.start.year,
              screen.stop.year=screen.stop.year,
              attendance=attendance,
              sensitivity=sensitivity)

test_that('Number of tests offered is sane', {
          expect_true(max(dset$tests_offered) <= sojourn.max)
          expect_true(max(dset$tests_offered) <= followup.years)
})

