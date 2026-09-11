# Preparation of Olink meta data

Preparation of Olink meta data

## Usage

``` r
prepare_olink_meta_data(
  data,
  raw_olink_protein_data,
  names = c("protein", "id"),
  columns = c("uniprot", "olinkid")
)
```

## Arguments

- data:

  Raw CKB Olink meta data (imported from olink_explore_meta.rds or
  olink_explore_meta.csv)

- raw_olink_protein_data:

  Raw Olink protein data (imported from data_baseline_olink_explore.rds
  or data_baseline_olink_explore.csv)

- names:

  Naming scheme to use for proteins: "protein" (default) to use assay
  names or "id" to use Olink IDs.

- columns:

  Optional columns to include in the meta data. (Default: c("uniprot",
  "olinkid"))
