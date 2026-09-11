# Associations of proteins in a case-subcohort study

This article demonstrates how to use this package to estimate
associations between proteins and risk of incidence disease using Cox
proportional hazards models, where the study uses a case-subcohort
design. For a real application, see [“Plasma Proteomics to Identify Drug
Targets for Ischemic Heart
Disease”](https://doi.org/10.1016/j.jacc.2023.09.804).

## Import data

First, all relevant participant and protein data must be imported into
R. For example, use [`readRDS()`](https://rdrr.io/r/base/readRDS.html)
for RDS format data or
[`read.csv()`](https://rdrr.io/r/utils/read.table.html) for csv data.

``` r

baseline       <- readRDS("baseline.rds")
ascertainments <- readRDS("ascertainments.rds")
endpoints      <- readRDS("endpoints.rds")

raw_olink_meta_data <- readRDS("olink_explore_meta.rds")
raw_olink_protein_data <- readRDS("data_baseline_olink_explore.rds")

raw_somascan_meta_data <- readRDS("somalogic_meta.rds")
raw_somascan_protein_data <- readRDS("data_baseline_somalogic.rds")
raw_somascan_sample_data  <- readRDS("data_baseline_somalogic_samples.rds")
```

For demonstration, the synthetic data sets available in this package can
be used:

1.  ckbtools_participant_data
2.  ckbtools_ascertainments
3.  ckbtools_raw_olink_meta_data
4.  ckbtools_raw_olink_protein_data
5.  ckbtools_raw_somascan_meta_data
6.  ckbtools_raw_somascan_protein_data
7.  ckbtools_raw_somascan_sample_data

### Baseline and endpoint data

A data set containing relevant participant baseline and endpoint data
must be prepared manually. See [CKB
data](https://ckbiobank.github.io/ckbtools/articles/ckb_data.md) for an
example.

In this example, the outcome is IHD.

### Identify case-subcohort participants

``` r

ckbtools_participant_data$subcohort <- ckbtools_participant_data$proteomics_subcohort

## keep only the subcohort and IHD cases outside the subcohort
ckbtools_participant_data <- ckbtools_participant_data[ckbtools_participant_data$subcohort == 1 |
                                                         ckbtools_participant_data$ihd == 1, ]
```

### Prentice weighting

To use the [Prentice weighting method for case-subcohort
designs](https://doi.org/10.1093/biomet/73.1.1), follow up time starts
at zero for participants in the subcohort but starts just before outcome
for cases outside subcohort.

``` r

ckbtools_participant_data$time_out <- as.numeric(ckbtools_participant_data$ihd_date - ckbtools_participant_data$study_date)
ckbtools_participant_data$time_in <- ckbtools_participant_data$time_out - 0.001
ckbtools_participant_data$time_in[ckbtools_participant_data$subcohort == 1] <- 0
```

## Prepare protein data

Prepare the protein data as described in
[Proteomics](https://ckbiobank.github.io/ckbtools/articles/proteomics.html#preparing-protein-data).

For this analysis, no plate correction is applied to the protein values
but we will adjust for “plateid” in the analysis models. SomaScan values
are log₂ transformed. Protein values are then scaled to have standard
deviations of one.

``` r

## set CKB data release
options(ckb.data.release = "19.02")

## meta data
olink_meta_data <- prepare_olink_meta_data(ckbtools_raw_olink_meta_data,
                                           ckbtools_raw_olink_protein_data)

somascan_meta_data <- prepare_somascan_meta_data(ckbtools_raw_somascan_meta_data)

## protein data
olink_protein_data <- prepare_olink_data(ckbtools_raw_olink_protein_data,
                                         olink_meta_data)

somascan_protein_data <- prepare_somascan_data(ckbtools_raw_somascan_protein_data,
                                               ckbtools_raw_somascan_sample_data,
                                               somascan_meta_data,
                                               transform = "log2")

## proteomics analysis data
olink_analysis_data <- prepare_proteomics_analysis_dataset(
  ckbtools_participant_data,
  olink_protein_data,
  olink_meta_data,
  plate_correction = "none",
  scaling = "scale"
)

somascan_analysis_data <- prepare_proteomics_analysis_dataset(
  ckbtools_participant_data,
  somascan_protein_data,
  somascan_meta_data,
  plate_correction = "none",
  scaling = "scale"
)
```

## Analysis

The
[`multiple_models()`](https://ckbiobank.github.io/ckbtools/reference/multiple_models.md)
function can be used to run multiple regression models.

*It is highly recommended to first run some individual analyses to check
that the models and output are as intended, before running analyses for
many proteins.*

The code below fits models for all Olink and SomaScan proteins. To set
the correct entry and exit times for the Prentice weighting use the
`time_in` and `time_out` arguments. Robust standard errors should also
be used by setting `clusted = "csid"`.

``` r

fully_adjusted <- c("strata(sex)",
                    "strata(region)",
                    "poly(age_at_study_date, 2)",
                    "poly(hours_since_last_ate, 2)",
                    "poly(region_mean_temp, 2)",
                    "education",
                    "smoking",
                    "alcohol",
                    "met",
                    "sbp_mean",
                    "has_diabetes",
                    "bmi_calc",
                    "plateid")


olink_results <- multiple_models(
  olink_analysis_data,
  model      = "cox",
  exposures  = olink_meta_data$name,
  outcomes   = "ihd",
  time_in    = "time_in",
  time_out   = "time_out",
  cluster    = "csid",
  adjust_for = list(full = fully_adjusted)
)

somascan_results <- multiple_models(
  somascan_analysis_data,
  model      = "cox",
  exposures  = somascan_meta_data$name,
  outcomes   = "ihd",
  time_in    = "time_in",
  time_out   = "time_out",
  cluster    = "csid",
  adjust_for = list(full = fully_adjusted)
)
```

*Note that education, smoking and alcohol are formatted as factors in
the data, so they are correctly used as categorical variables in the
model. If you want to include a categorical variable that is coded
numerically then make sure to use `"as.factor(variable)"`.*

## Summarise results

### Check for errors and warnings

After fitting multiple models, always check for warnings and errors. If
there are any, decide how these will be handled when summarising results
or for further analyses.

``` r

olink_results %>%
  dplyr::filter(!is.na(error) | !is.na(warnings)) %>%
  dplyr::distinct(outcome, exposure, adjustment, error, warnings)
```

``` outp
# # A tibble: 0 × 5
# # ℹ 5 variables: outcome <chr>, exposure <chr>, adjustment <chr>, error <chr>,
# #   warnings <chr>
```

``` r

somascan_results %>%
  dplyr::filter(!is.na(error) | !is.na(warnings)) %>%
  dplyr::distinct(outcome, exposure, adjustment, error, warnings)
```

``` outp
# # A tibble: 0 × 5
# # ℹ 5 variables: outcome <chr>, exposure <chr>, adjustment <chr>, error <chr>,
# #   warnings <chr>
```

### Select one result per protein

Some Olink assays are repeated across panels. To summarise the results
as in the referenced paper, select the assay by panel preference
Cardiometabolic \> Inflammation \> Neurology \> Oncology, and prefer
models with no error or warning.

``` r

olink_summary_results <- olink_results %>%
  dplyr::filter(exposure %in% olink_meta_data$name) %>%
  dplyr::inner_join(olink_meta_data, by = c("exposure" = "name")) %>%
  dplyr::group_by(outcome, adjustment, uniprot) %>%
  dplyr::arrange(!is.na(error) | !is.na(warnings), panel) %>%
  dplyr::filter(dplyr::row_number() == 1) %>%
  dplyr::ungroup()
```

Some SomaScan aptamers target the same protein. To summarise the results
as in the referenced paper, select one aptamer per target choosing the
result with the smallest p-value and no error or warning (or just the
smallest p-value if all have error or warning)

``` r

somascan_summary_results <- somascan_results %>%
  dplyr::inner_join(somascan_meta_data, by = c("exposure" = "name")) %>%
  dplyr::group_by(adjustment, target) %>%
  dplyr::arrange(!is.na(error) | !is.na(warnings), p.value) %>%
  dplyr::filter(dplyr::row_number() == 1) %>%
  dplyr::ungroup()
```

### Multiple testing

In this example, apply a false discovery rate correction separately in
the Olink and SomaScan analyses.

``` r

olink_summary_results <- olink_summary_results %>%
  ## handle errors and warnings
  dplyr::mutate(p.value = dplyr::if_else(!is.na(error) | !is.na(warnings), 1, p.value)) %>%
  ## multiple testing correction
  dplyr::mutate(p.value.fdr = p.adjust(p.value, method = "fdr"))

somascan_summary_results <- somascan_summary_results %>%
  ## handle errors and warnings
  dplyr::mutate(p.value = dplyr::if_else(!is.na(error) | !is.na(warnings), 1, p.value)) %>%
  ## multiple testing correction
  dplyr::mutate(p.value.fdr = p.adjust(p.value, method = "fdr"))
```

## Present results

The number of proteins significantly associated (FDR-corrected p \<
0.05) for Olink and SomaScan proteins can then be summarised.

``` r

olink_summary_results %>%
  dplyr::summarise(n = dplyr::n(),
                   error_or_warning = sum(!is.na(error) | !is.na(warnings)),
                   significant = sum(p.value.fdr < 0.05),
                   significant_inverse = sum(p.value.fdr < 0.05 & estimate < 0),
                   .groups = "drop")
```

``` outp
# # A tibble: 1 × 4
#       n error_or_warning significant significant_inverse
#   <int>            <int>       <int>               <int>
# 1    48                0           1                   0
```

``` r

somascan_summary_results %>%
  dplyr::summarise(n = dplyr::n(),
                   error_or_warning = sum(!is.na(error) | !is.na(warnings)),
                   significant = sum(p.value.fdr < 0.05),
                   significant_inverse = sum(p.value.fdr < 0.05 & estimate < 0),
                   .groups = "drop")
```

``` outp
# # A tibble: 1 × 4
#       n error_or_warning significant significant_inverse
#   <int>            <int>       <int>               <int>
# 1    48                0          48                   0
```

``` r

olink_summary_results %>%
  dplyr::arrange(p.value.fdr) %>%
  dplyr::select(assay, uniprot, estimate, robust.se, p.value, p.value.fdr) %>%
  head(5)
```

``` outp
# # A tibble: 5 × 6
#   assay  uniprot estimate robust.se   p.value p.value.fdr
#   <chr>  <chr>      <dbl>     <dbl>     <dbl>       <dbl>
# 1 MYL4   P12829     0.304    0.0744 0.0000449     0.00216
# 2 FOLH1  Q04609    -0.166    0.0780 0.0336        0.403  
# 3 CNTF   P26441     0.190    0.0803 0.0177        0.403  
# 4 GABRA4 P48169    -0.160    0.0741 0.0311        0.403  
# 5 TIA1   P31483    -0.140    0.0761 0.0654        0.539
```

``` r

somascan_summary_results %>%
  dplyr::arrange(p.value.fdr) %>%
  dplyr::select(target, uniprot, estimate, robust.se, p.value, p.value.fdr) %>%
  head(5)
```

``` outp
# # A tibble: 5 × 6
#   target uniprot estimate robust.se  p.value p.value.fdr
#   <chr>  <chr>      <dbl>     <dbl>    <dbl>       <dbl>
# 1 PIN1   Q13526      26.3      2.87 4.64e-20    1.91e-19
# 2 EI2BA  Q14232      26.0      2.84 5.29e-20    1.91e-19
# 3 SGK3   Q96BR1      26.2      2.86 5.89e-20    1.91e-19
# 4 SLUG   O43623      26.2      2.87 6.23e-20    1.91e-19
# 5 BRD2   P25440      25.9      2.84 6.66e-20    1.91e-19
```
