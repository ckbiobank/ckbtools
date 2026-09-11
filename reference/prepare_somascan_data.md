# Prepare SomaScan protein data

This functions prepares a data set of SomaScan protein data given the
raw CKB SomaScan data.

## Usage

``` r
prepare_somascan_data(
  data,
  sample_data,
  meta_data,
  transform = c("log2", "log", "log10", "none"),
  qc = c("default", "flagged", "none"),
  repeat_samples = c("remove", "retain"),
  columns = NULL
)
```

## Arguments

- data:

  Raw CKB SomaScan protein data (imported from
  data_baseline_somalogic.rds or data_baseline_somalogic.csv)

- sample_data:

  Raw CKB SomaScan protein sample information data (imported from
  data_baseline_somalogic_samples.rds or
  data_baseline_somalogic_samples.csv)

- meta_data:

  SomaScan meta data produced by
  [`prepare_somascan_meta_data()`](https://ckbiobank.github.io/ckbtools/reference/prepare_somascan_meta_data.md).

- transform:

  Transformation to apply to protein values:

  - "log2" (default): log base 2

  - "log10": log base 10

  - "none"

- qc:

  QC steps to apply.

  - "default" (default): Exclude participant samples with hybridization
    control scale factor out of range.

  - "flagged": Exclude participants samples flagged for any
    normalization acceptance criteria.

  - "none": Do not exclude any participant samples.

- repeat_samples:

  Default "remove" will remove data for the second sample for
  participants who had a sample run twice. Use "retain" to keep this
  data; the second samples will have sample_source=1.

- columns:

  Names of any additional columns from sample_data to include in the
  output.

## QC acceptance criteria

Hybridization control scale factor (`hybcontrolnormscale`) and median
signal normalization scale factors (`normscale_20`, `normscale_0_5`,
`normscale_0_005`) have an acceptance range of 0.4-2.5.

Percentage of measurements used for ANML scale factor calculation
(`anmlfractionused_20`, `anmlfractionused_0_5`,
`anmlfractionused_0_005`) has a criteria of \> 30% (0.3).

Samples that meet these criteria are given a PASS for normalization
acceptance criteria for all row scale factors (`rowcheck`), while
samples that are outside the criteria are given a FLAG.

Hybridization control scale factor out of range suggests a technical
issue. Set `qc = "default"` to exclude participant samples with
hybridization control scale factor (`hybcontrolnormscale`) out of range
(i.e. less than 0.4 or greater than 2.5).

Other criteria out of range suggests that the sample is fairly different
from a typical population. Set `qc = "flagged"` to exclude participant
samples with a FLAG for normalization acceptance criteria for all row
scale factors (`rowcheck`).
