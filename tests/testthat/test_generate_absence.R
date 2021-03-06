context('Simulated incidence in the absence of screening')

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

test_that('Sojourn times have correct range', {
          sojourn.time.range <- seq(sojourn.min, sojourn.max)
          expect_equal(unique(dset$sojourn), sojourn.time.range)
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
          expect_equal(nrow(dset), length(sojourn.time.range)*length(onset.range))
})

