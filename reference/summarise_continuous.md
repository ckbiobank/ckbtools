# Summarise continuous variables

Summarise continuous variables

## Usage

``` r
summarise_continuous(data, cols, by = NULL)
```

## Arguments

- data:

  Data set (e.g. produced by
  [`prepare_proteomics_analysis_dataset()`](https://ckbiobank.github.io/ckbtools/reference/prepare_proteomics_analysis_dataset.md))

- cols:

  Names of columns containing continuous variables to be summarised.

- by:

  Names of columns that contain categorical participant characteristics
  to produce subgroup summaries.
