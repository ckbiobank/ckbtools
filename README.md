
# ckbtools <img src="man/figures/logo.png" align="right" width = "120" />

`ckbtools` provides tools for statistical analyses in epidemiology.

It is developed by [China Kadoorie Biobank](http://www.ckbiobank.org) researchers for users of CKB data and other researchers in epidemiology.

Features of the package include:

- Standardised incidence rates
- Age-at-risk expansion
- Variance of the log risk
- Trend and heterogeneity tests
- Preparation of CKB proteomics data
- Regression model helper functions


## Installation

Install the latest version of `ckbtools` from the ckbiobank R-universe:

``` r
install.packages('ckbtools',
                 repos = c('https://ckbiobank.r-universe.dev/', getOption('repos')))
```

This will also install required dependencies using your existing package repositories.

## Get started

Read `vignette("ckbtools")` to get started.

## Citation

If you find this package useful, please cite as:

China Kadoorie Biobank Collaborative Group (2026). "ckbtools: Tools for analyses in epidemiology." <https://ckbiobank.github.io/ckbtools>.


## Development and support

This package is under development and will be updated with bug fixes and additional functions related to statistical analyses in epidemiology. Please refer to the package documentation for guidance on usage and use GitHub issues to report bugs. This package is provided “as is” for research purposes, without warranty of any kind. The authors accept no liability for its use or for any results derived from it.
