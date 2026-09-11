# Preparation of SomaScan meta data

Preparation of SomaScan meta data

## Usage

``` r
prepare_somascan_meta_data(
  data,
  names = c("protein", "id"),
  columns = c("seqid", "somaid", "target_full_name", "target", "uniprot", "organism",
    "type", "dilution"),
  filter = "HumanProtein"
)
```

## Arguments

- data:

  Raw CKB SomaScan meta data (imported from somalogic_meta.rds or
  somalogic_meta.csv)

- names:

  Naming scheme to use for proteins: "protein" (default) to use target
  name or "id" to use SomaLogic sequence IDs.

- columns:

  Optional columns to include in the meta data. (Default: c("seqid",
  "somaid", "target_full_name", "target", "uniprot", "organism", "type",
  "dilution"))

- filter:

  Filter to apply to rows:

  - "HumanProtein" (default): Only return meta data for aptamers where
    organism is "Human" and type is "Protein".

  - "All": Return all rows of meta data.
