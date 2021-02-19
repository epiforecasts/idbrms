context("id_stancode.idbrms_convolution")

if (!testthat:::on_cran()) {
  files <- c("discretised_lognormal_pmf.stan",
             "calc_pmf.stan",
             "calc_unique_pmfs.stan")
  suppressMessages(
    expose_idbrms_stan_fns(
      files, dir =  system.file("stan/functions", package = "idbrms")
      )
    )
}


test_that("calc_unique_pmfs: Successfully calculates a singe PMF", {
  expect_equal(
    round(calc_unique_pmfs(1.5, 0.2, 10)[[1]], 2),
    c(0.00, 0.00, 0.01, 0.06, 0.22, 0.42, 0.26, 0.02, 0.00, 0.00)
  )
  expect_equal(
    round(calc_unique_pmfs(1.5, 0.2, 6)[[1]], 2),
    c(0.24, 0.46, 0.28, 0.02, 0.00, 0.00)
  )
})

test_that("calc_unique_pmfs: Successfully calculates multiple PMFs with unique
         summary statistics", {
  pmfs <- calc_unique_pmfs(c(1.5, 1.2, 1.9, 0.2), c(0.2, 0.1, 0.4, 0.6), 10)
  expect_equal(length(unique(pmfs)), 4)
})

test_that("calc_unique_pmfs: Successfully calculates multiple PMFs with
           duplicate summary statistics", {
           pmfs <- calc_unique_pmfs(c(1.5, 1.2, 1.5, 1.2),
                                   c(0.2, 0.1, 0.2, 0.1), 10)
           expect_equal(length(unique(pmfs)), 2)
         })
