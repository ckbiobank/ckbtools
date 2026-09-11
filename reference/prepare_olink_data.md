# Preparation of Olink protein data

Preparation of Olink protein data

## Usage

``` r
prepare_olink_data(
  data,
  meta_data,
  transform = c("none"),
  qc = c("default", "none"),
  columns = NULL
)
```

## Arguments

- data:

  Raw Olink protein data (imported from data_baseline_olink_explore.rds
  or data_baseline_olink_explore.csv)

- meta_data:

  Olink meta data produced by
  [`prepare_olink_meta_data()`](https://ckbiobank.github.io/ckbtools/reference/prepare_olink_meta_data.md).

- transform:

  Transformation to apply to protein values. Currently only "none" is
  available.

- qc:

  QC steps to apply:

  - "default" (default): replace protein values with a QC warning with a
    missing value

  - "none"

- columns:

  Other columns in raw data to include in output. For example, set to
  "lod" to include the limit of detection for each protein.
