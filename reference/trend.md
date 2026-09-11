# Test for trend

A chi-squared test for trend (on 1 degree of freedom) given a vector of
independent estimates and their standard errors.

## Usage

``` r
trend(beta, se)
```

## Arguments

- beta:

  vector of independent estimates

- se:

  vector of their standard errors

## Value

A named list:

- test statistic (chi-squared test statistic)

- p (p-value)

## Details

See
[`vignette("trend_heterogeneity")`](https://ckbiobank.github.io/ckbtools/articles/trend_heterogeneity.md)
for details.
