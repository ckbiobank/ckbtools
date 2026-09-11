# Test for heterogeneity

A chi-squared test (on \\k - 1\\ degrees of freedom) for heterogeneity
given a vector of \\k\\ independent estimates and their standard errors.

## Usage

``` r
heterogeneity(beta, se)
```

## Arguments

- beta:

  vector of independent estimates

- se:

  vector of their standard errors

## Value

A named list:

- test statistic (chi-squared test statistic)

- degrees of freedom (number of degrees of freedom used for the test,
  equal to `length(beta)-1`)

- p (p-value)

## Details

See
[`vignette("trend_heterogeneity")`](https://ckbiobank.github.io/ckbtools/articles/trend_heterogeneity.md)
for details.
