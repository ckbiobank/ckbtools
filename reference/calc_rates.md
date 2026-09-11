# Calculate rates and their standard errors using various methods

Calculate rates and their standard errors using various methods

## Usage

``` r
calc_rates(data, per_pyar, method = c("norm", "gamma"))
```

## Arguments

- data:

  A data frame.

- per_pyar:

  Crude and standardised incidence rates are given per `per_pyar`
  person-years at risk. (Default: 1000)

- method:

  Method to calculate confidence intervals, one of "norm" or "gamma".
  Default "norm" is Normal approximation.

## Value

A data frame
