# Preparation of a proteomics analysis data set

Preparation of a proteomics analysis data set

## Usage

``` r
prepare_proteomics_analysis_dataset(
  participant_data,
  protein_data,
  protein_meta_data,
  plate_correction = c("none", "median-subcohort", "mean-subcohort"),
  scaling = c("none", "center-subcohort", "scale-subcohort", "center", "scale", "rint"),
  id = "csid"
)
```

## Arguments

- participant_data:

  A data frame of participant data

- protein_data:

  Protein data produced by
  [`prepare_olink_data()`](https://ckbiobank.github.io/ckbtools/reference/prepare_olink_data.md)
  or
  [`prepare_somascan_data()`](https://ckbiobank.github.io/ckbtools/reference/prepare_somascan_data.md).

- protein_meta_data:

  Corresponding protein meta data produed by
  [`prepare_olink_meta_data()`](https://ckbiobank.github.io/ckbtools/reference/prepare_olink_meta_data.md)
  or
  [`prepare_somascan_meta_data()`](https://ckbiobank.github.io/ckbtools/reference/prepare_somascan_meta_data.md).

- plate_correction:

  Type of plate correction to apply:

  - "none" (default)

  - "median-subcohort": For each protein, subtract the plate median (of
    subcohort values) and add the overall median (of subcohort values).

  - "mean-subcohort": For each protein, subtract the plate mean (of
    subcohort values) and add the overall mean (of subcohort values).

- scaling:

  Scaling to apply to values for each protein:

  - "none" (default)

  - "center-subcohort": Subtract the overall mean of subcohort values.

  - "scale-subcohort": Subtract the overall mean and divide by the
    sample standard deviation of subcohort values.

  - "center": Subtract the overall mean.

  - "scale": Subtract the overall mean and divide by the sample standard
    deviation.

  - "rint": Rank inverse normal transformation.

- id:

  Name of column used to join participant_data and protein_data.
  (Default: "csid")
