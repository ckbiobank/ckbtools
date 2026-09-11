# Calculate standardised incidence rates

Calculate standardised incidence rates

## Usage

``` r
incidence_rates(
  data,
  col.dob = "dob_anon",
  col.entry_date = "study_date",
  col.exit_date = "endpoint_date",
  col.status = "endpoint",
  ages = c(0, 100),
  group = NULL,
  adjust_for = NULL,
  per_pyar = 1000,
  weights = NULL,
  method = c("norm", "gamma")
)
```

## Arguments

- data:

  A data frame.

- col.dob:

  Name of column containing date of birth. (Must be of class `Date`.)

- col.entry_date:

  Name of column containing entry date. (Must be of class `Date`.)

- col.exit_date:

  Name of column containing exit date. (Must be of class `Date`.)

- col.status:

  Name of column containing status. This should be 1 if the event has
  occurred at the exit date, or 0 if it has not.

- ages:

  Vector of ages (years) to create age-at-risk strata. (Must be a
  numeric vector with distinct, increasing values. Default: `c(0, 100)`)

- group:

  Name of column containing groups to separately calculate incidence
  rates.

- adjust_for:

  Vector of names of columns containing adjustment variables. Include
  "agegrp" to adjust for age-at-risk strata.

- per_pyar:

  Crude and standardised incidence rates are given per `per_pyar`
  person-years at risk. (Default: 1000)

- weights:

  A data frame with columns that match the `adjust_for` variables, and
  one row for each combination of levels of those variables, and a
  column named ‘std_pyar’ with the weights (person-years at risk for
  that stratum in the reference/standard population).

- method:

  Method to calculate confidence intervals, one of "norm" or "gamma".
  Default "norm" is Normal approximation.

## Value

A named list:

- incidence_rates_in_all_strata

- standardised_incidence_rates

## Methods

See
[`vignette("incidence_rates")`](https://ckbiobank.github.io/ckbtools/articles/incidence_rates.md)
for details.
