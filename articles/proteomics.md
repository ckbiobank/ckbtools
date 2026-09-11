# Proteomics

## Introduction to CKB proteomics data

The [China Kadoorie Biobank](http://www.ckbiobank.org) includes
proteomics data from a nested case-cohort study of ~2k MI cases and ~2k
subcohort participants. Plasma samples were analysed using two
proteomics platforms: [Olink Explore (~3k
proteins)](https://olink.com/products/olink-explore-series) and
[SomaScan (~7k proteins)](https://somalogic.com/somascan-platform/).

Details of the CKB Data Sharing Policy, data release schedules and data
request application procedures are available at
[www.ckbiobank.org/data-access](https://www.ckbiobank.org/data-access).

## Observational proteomics analyses

The functions in this package can be used to prepare analysis-ready data
sets from raw CKB proteomics data and perform some observational
proteomics analyses as shown in these articles:

- **[Associations of proteins with participant
  characteristics](https://ckbiobank.github.io/ckbtools/articles/proteomics_characteristics.md)**
  demonstrates how to use this package to estimate associations between
  proteins and participant characteristics. The methods used are based
  on the papers [“Associations of 2923 Olink proteins with demographic,
  lifestyle, environmental and health characteristics in middle-aged
  Chinese adults”](https://doi.org/10.1007/s10654-025-01311-z) and
  [“Associations of 6600 SomaScan proteins with demographic, lifestyle,
  environmental and health characteristics in Chinese
  adults”](https://doi.org/10.1038/s41598-026-41444-z).

- **[Associations of proteins with incident
  disease](https://ckbiobank.github.io/ckbtools/articles/proteomics_diseases.md)**
  demonstrates how to use this package to estimate associations between
  proteins and risk of incident disease using Cox proportional hazards
  models. The methods used are from [“Proteome-wide assessment of ~9500
  plasma proteins with major diseases and multi-morbidity in a Chinese
  population”](https://doi.org/10.21203/rs.3.rs-5356315/v1).

- **[Associations of proteins in a case-subcohort
  study](https://ckbiobank.github.io/ckbtools/articles/proteomics_case-subcohort.md)**
  demonstrates how to use this package to estimate associations between
  proteins and risk of incidence disease using Cox proportional hazards
  models with methods appropriate for a case-subcohort study design. The
  methods used are from [“Plasma Proteomics to Identify Drug Targets for
  Ischemic Heart Disease”](https://doi.org/10.1016/j.jacc.2023.09.804).

## Preparing proteomics data

This package includes several functions to assist with preparing the CKB
proteomics data for analysis.

### Import data

First, all relevant protein data must be imported into R. For example,
use [`readRDS()`](https://rdrr.io/r/base/readRDS.html) for RDS format
data or [`read.csv()`](https://rdrr.io/r/utils/read.table.html) for csv
data.

``` r

raw_olink_meta_data <- readRDS("olink_explore_meta.rds")
raw_olink_protein_data <- readRDS("data_baseline_olink_explore.rds")

raw_somascan_meta_data <- readRDS("somalogic_meta.rds")
raw_somascan_protein_data <- readRDS("data_baseline_somalogic.rds")
raw_somascan_sample_data  <- readRDS("data_baseline_somalogic_samples.rds")
```

For demonstration, the synthetic data sets available in this package can
be used:

1.  ckbtools_raw_olink_meta_data
2.  ckbtools_raw_olink_protein_data
3.  ckbtools_raw_somascan_meta_data
4.  ckbtools_raw_somascan_protein_data
5.  ckbtools_raw_somascan_sample_data

### Prepare protein data

To use the proteomics data preparation functions, the CKB data release
number must be set using
[`options()`](https://rdrr.io/r/base/options.html). The release number
for CKB data can be found in the release.txt file. This ensures the
correct code is used for the format of the supplied data.

``` r

options(ckb.data.release = "19.02")
```

To prepare protein meta data, use the
[`prepare_olink_meta_data()`](https://ckbiobank.github.io/ckbtools/reference/prepare_olink_meta_data.md)
and
[`prepare_somascan_meta_data()`](https://ckbiobank.github.io/ckbtools/reference/prepare_somascan_meta_data.md)
functions. These return data frames with one row per assay/aptamer and
columns for protein details such as UniProt ID. The `names` argument
specifies the naming system used for proteins: “protein” for names based
on assay (for Olink) or target (for SomaScan); or “id” for names based
on OlinkID (Olink ID) or SeqID (SomaLogic sequence identifier).

``` r

olink_meta_data <- prepare_olink_meta_data(ckbtools_raw_olink_meta_data,
                                           ckbtools_raw_olink_protein_data,
                                           names = "id")
somascan_meta_data <- prepare_somascan_meta_data(ckbtools_raw_somascan_meta_data,
                                                 names = "id")
```

To prepare protein data, use the
[`prepare_olink_data()`](https://ckbiobank.github.io/ckbtools/reference/prepare_olink_data.md)
and
[`prepare_somascan_data()`](https://ckbiobank.github.io/ckbtools/reference/prepare_somascan_data.md)
functions. The functions each create a ‘wide’ data set with one row per
sample/participant and columns for each protein measure.

``` r

olink_protein_data <- prepare_olink_data(ckbtools_raw_olink_protein_data,
                                         olink_meta_data)
```

``` outp
# ✔ Protein values with QC warnings have been replaced with missing values.
# [Default]
```

The `transform` argument can be used to apply a transformation to all
SomaScan protein values (for example `"log2"` will log₂ transform the
RFU values).

``` r

somascan_protein_data <- prepare_somascan_data(ckbtools_raw_somascan_protein_data,
                                               ckbtools_raw_somascan_sample_data,
                                               somascan_meta_data,
                                               transform = "log2")
```

``` outp
# ✔ RFU values have been log transformed (base 2). [Default]
```

``` outp
# ✔ Participant samples with hybcontrolnormscale out of range (i.e. less than 0.4
# or greater than 2.5) have been removed. [Default]
```

``` outp
# ✔ Data for the second sample (for participants who had a sample run twice) have
# been removed. [Default]
```

To merge participant data with protein data, use the
[`prepare_proteomics_analysis_dataset()`](https://ckbiobank.github.io/ckbtools/reference/prepare_proteomics_analysis_dataset.md)
function. Only observations in both the participant data and proteomics
data sets are retained. There must be a column in the participant data
named “subcohort” that identifies observations that are part of the
subcohort.

``` r

ckbtools_participant_data$subcohort <- ckbtools_participant_data$proteomics_subcohort
```

The `plate_correction` argument can be used to correct protein values
for plate effects. Set this to `"median-subcohort"`/`"mean-subcohort"`
to centre on medians/medians across plates based on subcohort values.

The `scaling` argument can be used to scale protein values. For example,
`"scale-subcohort"` will scale values (as in the R function
[`scale()`](https://rdrr.io/r/base/scale.html)) using the mean and
standard deviation of subcohort values. `"rint"` will apply a rank
inverse normal transformation to all values.

``` r

olink_analysis_data <- prepare_proteomics_analysis_dataset(
  participant_data  = ckbtools_participant_data,
  protein_data      = olink_protein_data,
  protein_meta_data = olink_meta_data,
  plate_correction  = "median-subcohort",
  scaling           = "none"
)
```

``` outp
# ✔ Plate correction 'median-subcohort' has been applied to protein values.
```

``` outp
# ✔ No scaling of protein values has been applied. [Default]
```

``` r

somascan_analysis_data <- prepare_proteomics_analysis_dataset(
  participant_data  = ckbtools_participant_data,
  protein_data      = somascan_protein_data,
  protein_meta_data = somascan_meta_data,
  plate_correction  = "none",
  scaling           = "none"
)
```

``` outp
# ✔ No plate correction has been applied to protein values. [Default]
# ✔ No scaling of protein values has been applied. [Default]
```

In summary, the proteomics data preparation steps in order are:

1.  Log transformation (using the `transform` argument in
    [`prepare_somascan_data()`](https://ckbiobank.github.io/ckbtools/reference/prepare_somascan_data.md)).
2.  Removing values/samples with QC warnings (using the `qc` argument in
    [`prepare_olink_data()`](https://ckbiobank.github.io/ckbtools/reference/prepare_olink_data.md)
    and
    [`prepare_somascan_data()`](https://ckbiobank.github.io/ckbtools/reference/prepare_somascan_data.md)).
3.  Removing repeated samples (using the repeat_samples argument in
    [`prepare_somascan_data()`](https://ckbiobank.github.io/ckbtools/reference/prepare_somascan_data.md)).
4.  Merging proteomics and participant data (using the
    `participant_data` and `protein_data arguments` in
    [`prepare_proteomics_analysis_dataset()`](https://ckbiobank.github.io/ckbtools/reference/prepare_proteomics_analysis_dataset.md)).
    Note that this is an “inner join” so that only matching observations
    in both data sets are kept.
5.  Plate correction (using the `plate_correction` argument in
    [`prepare_proteomics_analysis_dataset()`](https://ckbiobank.github.io/ckbtools/reference/prepare_proteomics_analysis_dataset.md)).
6.  Scaling protein values (using the `scaling` argument in
    [`prepare_proteomics_analysis_dataset()`](https://ckbiobank.github.io/ckbtools/reference/prepare_proteomics_analysis_dataset.md)).

### Descriptive summary statistics

Use
[`summarise_continuous()`](https://ckbiobank.github.io/ckbtools/reference/summarise_continuous.md)
to create basic descriptive summary statistics for proteins.

``` r

protein_summary_stats <- summarise_continuous(olink_analysis_data,
                                              cols = olink_meta_data$name,
                                              by = "sex")

protein_summary_stats[1:12, 1:7]
```

``` outp
#         var    sex n_nonmissing n_missing percent_missing       mean        sd
# 1  OID13950   <NA>          987        13       1.3000000 0.17825730 0.7761187
# 2  OID13950   Male          413         6       1.4319809 0.18644455 0.7509019
# 3  OID13950 Female          574         7       1.2048193 0.17236647 0.7943648
# 4  OID01935   <NA>          988        12       1.2000000 0.16706032 0.8281963
# 5  OID01935   Male          413         6       1.4319809 0.19642842 0.8623041
# 6  OID01935 Female          575         6       1.0327022 0.14596636 0.8029020
# 7  OID12648   <NA>          993         7       0.7000000 0.13701380 0.8096633
# 8  OID12648   Male          416         3       0.7159905 0.10770473 0.7944150
# 9  OID12648 Female          577         4       0.6884682 0.15814478 0.8205168
# 10 OID91452   <NA>          994         6       0.6000000 0.07626006 0.8198431
# 11 OID91452   Male          415         4       0.9546539 0.09757931 0.8470904
# 12 OID91452 Female          579         2       0.3442341 0.06097942 0.8001332
```

### Olink: Limit of detection

The number and percentages of Olink protein values below the limit of
detection can be obtained as follows:

``` r

olink_meta_data <- prepare_olink_meta_data(ckbtools_raw_olink_meta_data,
                                           ckbtools_raw_olink_protein_data)
olink_protein_data <- prepare_olink_data(ckbtools_raw_olink_protein_data,
                                         olink_meta_data,
                                         columns = "lod")
```

``` outp
# ✔ Protein values with QC warnings have been replaced with missing values.
# [Default]
```

``` r

olink_meta_data %>%
  dplyr::arrange(name) %>%
  dplyr::select(name) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(n = sum(!is.na(olink_protein_data[name])),
                below_lod = sum(olink_protein_data[name] <
                                  olink_protein_data[paste0("lod_", name)], na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p = 100 * below_lod / n)
```

``` outp
# # A tibble: 48 × 4
#    name       n below_lod     p
#    <chr>  <int>     <int> <dbl>
#  1 ACAA1    991       205  20.7
#  2 BAMBI    994       187  18.8
#  3 BCR      994       206  20.7
#  4 CHEK2    988       213  21.6
#  5 CLASP1   991       188  19.0
#  6 CNGB3    986       205  20.8
#  7 CNTF     988       223  22.6
#  8 COMT     988       199  20.1
#  9 DAND5    983       204  20.8
# 10 DAPP1    984       199  20.2
# # ℹ 38 more rows
```
