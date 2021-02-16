
<!-- README.md is generated from README.Rmd. Please edit that file -->

# idbrms: Infectious Disease Modelling using brms

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Codecov test
coverage](https://codecov.io/gh/epiforecasts/brms.ide/branch/master/graph/badge.svg)](https://codecov.io/gh/epiforecasts/brms.ide?branch=master)

Provides population-level infectious disease models as an extension of
`brms`.

## Installation

You can install the unstable development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("epiforecasts/idbrms")
```

## Example

  - Load `idbrms`, `brms` and `data.table` (used for manipulating data).

<!-- end list -->

``` r
library(idbrms)
library(brms)
library(data.table)
```

  - Define toy data with a fixed 5 day lag between observations.

<!-- end list -->

``` r
dt <- data.table(
  region = "France", cases = seq(10, 500, by = 10),
  date = seq(as.Date("2020-10-01"), by = "days", length.out = 50))
dt[, deaths := as.integer(shift(cases, 5) * 0.1)]
dt[is.na(deaths), deaths := 0]
```

  - Prepare the data to be fit using the convolution model.

<!-- end list -->

``` r
dt <- prepare(dt, model = "convolution", location = "region",
              primary = "cases", secondary = "deaths")
head(dt, 10)
#>     location       date time index init_obs cstart cmax primary secondary
#>  1:   France 2020-10-01    0     1        1      1    1      10         0
#>  2:   France 2020-10-02    1     2        1      1    2      20         0
#>  3:   France 2020-10-03    2     3        1      1    3      30         0
#>  4:   France 2020-10-04    3     4        1      1    4      40         0
#>  5:   France 2020-10-05    4     5        1      1    5      50         0
#>  6:   France 2020-10-06    5     6        1      1    6      60         1
#>  7:   France 2020-10-07    6     7        1      1    7      70         2
#>  8:   France 2020-10-08    7     8        1      1    8      80         3
#>  9:   France 2020-10-09    8     9        1      1    9      90         4
#> 10:   France 2020-10-10    9    10        1      1   10     100         5
```

  - Fit the model.

<!-- end list -->

``` r
fit <- idbrm(data = dt)
summary(fit)
```

  - Expose model stan functions.

<!-- end list -->

``` r
expose_functions(fit)
```

  - Posterior predictive check. Runs without error but looks unlikely to
    be correct.

<!-- end list -->

``` r
pp_check(fit)
#> Using 10 posterior samples for ppc type 'dens_overlay' by default.
```

<img src="man/figures/README-pred-check-1.png" width="100%" />

  - Plot conditional effects (fails with due to mismatched vector sizes
    in the stan dot product + there are no variables in this model so
    should always fail).

<!-- end list -->

``` r
plot(conditional_effects(fit), ask = FALSE)
```
