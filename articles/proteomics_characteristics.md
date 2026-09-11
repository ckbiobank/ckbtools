# Associations of proteins with participant characteristics

This article demonstrates how to use this package to estimate
cross-sectional associations between proteins and participant
characteristics. The methods used are based on the papers [“Associations
of 2923 Olink proteins with demographic, lifestyle, environmental and
health characteristics in middle-aged Chinese
adults”](https://doi.org/10.1007/s10654-025-01311-z) and [“Associations
of 6600 SomaScan proteins with demographic, lifestyle, environmental and
health characteristics in Chinese
adults”](https://doi.org/10.1038/s41598-026-41444-z).

## Import data

First, all relevant participant and protein data must be imported into
R. For example, use [`readRDS()`](https://rdrr.io/r/base/readRDS.html)
for RDS format data or
[`read.csv()`](https://rdrr.io/r/utils/read.table.html) for csv data.

``` r

baseline       <- readRDS("baseline.rds")
ascertainments <- readRDS("ascertainments.rds")

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

## Preparation of participant data

A data set containing relevant participant baseline data must be
prepared manually. For these analyses, the data must contain all
participant characteristics used for exclusions or as adjustment
variables. See [CKB
data](https://ckbiobank.github.io/ckbtools/articles/ckb_data.md) for an
example.

### Select subcohort participants

These analyses of baseline characteristics will use only participants in
the subcohort of the nested case-subcohort proteomics study.

``` r

subcohort_participant_data <- ckbtools_participant_data[ckbtools_participant_data$proteomics_subcohort == 1, ]
subcohort_participant_data$subcohort <- 1
```

## Prepare protein data

Prepare the protein data as described in
[Proteomics](https://ckbiobank.github.io/ckbtools/articles/proteomics.html#preparing-protein-data).

``` r

options(ckb.data.release = "19.02")

## meta data
olink_meta_data <- prepare_olink_meta_data(ckbtools_raw_olink_meta_data,
                                           ckbtools_raw_olink_protein_data)

somascan_meta_data <- prepare_somascan_meta_data(ckbtools_raw_somascan_meta_data)

## protein data
olink_protein_data <- prepare_olink_data(ckbtools_raw_olink_protein_data,
                                         olink_meta_data)
```

``` outp
# ✔ Protein values with QC warnings have been replaced with missing values.
# [Default]
```

``` r

somascan_protein_data <- prepare_somascan_data(ckbtools_raw_somascan_protein_data,
                                               ckbtools_raw_somascan_sample_data,
                                               somascan_meta_data,
                                               transform = "log")
```

``` outp
# ✔ RFU values have been (natural) log transformed.
```

``` outp
# ✔ Participant samples with hybcontrolnormscale out of range (i.e. less than 0.4
# or greater than 2.5) have been removed. [Default]
```

``` outp
# ✔ Data for the second sample (for participants who had a sample run twice) have
# been removed. [Default]
```

``` r

## proteomics analysis data
olink_analysis_data <- prepare_proteomics_analysis_dataset(
  subcohort_participant_data,
  olink_protein_data,
  olink_meta_data,
  plate_correction = "none",
  scaling = "scale"
)
```

``` outp
# ✔ No plate correction has been applied to protein values. [Default]
```

``` outp
# ✔ Protein values have been centred and scaled using overall means and standard
# deviations.
```

``` r

somascan_analysis_data <- prepare_proteomics_analysis_dataset(
  subcohort_participant_data,
  somascan_protein_data,
  somascan_meta_data,
  plate_correction = "none",
  scaling = "scale"
)
```

``` outp
# ✔ No plate correction has been applied to protein values. [Default]
# ✔ Protein values have been centred and scaled using overall means and standard
# deviations.
```

## One model

The
[`multiple_models()`](https://ckbiobank.github.io/ckbtools/reference/multiple_models.md)
function can be used to fit a single or multiple linear regression
models. The code below will return a the results of a linear regression
of protein ACAA1 (Olink) on BMI (kg/m²). The model is adjusted for age
(and its square), sex, region, fasting time (and its square), outdoor
temperature (and its square) and plate.

``` r

basic_adjustment_terms <- c("age_at_study_date",
                            "I(age_at_study_date^2)",
                            "sex",
                            "region",
                            "hours_since_last_ate",
                            "I(hours_since_last_ate^2)",
                            "region_mean_temp",
                            "I(region_mean_temp^2)",
                            "plateid")

result <- multiple_models(
  olink_analysis_data,
  model = "linear",
  exposures = "bmi_calc",
  outcomes = "ACAA1",
  adjust_for = list(basic = basic_adjustment_terms)
)

result[, c("outcome", "term", "estimate", "std.error", "p.value", "n")]
```

``` outp
# # A tibble: 43 × 6
#    outcome term                    estimate std.error p.value     n
#    <chr>   <chr>                      <dbl>     <dbl>   <dbl> <int>
#  1 ACAA1   bmi_calc                0.0231    0.0152   0.129     496
#  2 ACAA1   (Intercept)            -2.73      1.26     0.0303     NA
#  3 ACAA1   age_at_study_date       0.0407    0.0433   0.348      NA
#  4 ACAA1   I(age_at_study_date^2) -0.000357  0.000404 0.377      NA
#  5 ACAA1   sexFemale              -0.0578    0.0960   0.547      NA
#  6 ACAA1   regionHarbin            0.237     0.223    0.288      NA
#  7 ACAA1   regionHaikou            0.275     0.255    0.282      NA
#  8 ACAA1   regionSuzhou            0.464     0.235    0.0487     NA
#  9 ACAA1   regionLiuzhou           0.670     0.242    0.00579    NA
# 10 ACAA1   regionSichuan           0.450     0.226    0.0473     NA
# # ℹ 33 more rows
```

By default, the function returns one row for a continuous exposure or
one row per group for a categorical exposure. Binary or categorical
exposures stored as numeric are treated as continuous. These should be
changed to factors to treat as categorical.

The code below returns the results of a linear regression of protein
CRBB2 (SomaScan) on BMI (grouped) adjusted for age, sex and region.

``` r

result <- multiple_models(
  somascan_analysis_data,
  model = "linear",
  exposures = "bmi_grp",
  outcomes = "CRBB2",
  adjust_for = list(basic = basic_adjustment_terms)
)

result[, c("outcome", "term", "estimate", "std.error", "p.value", "n", "n_group")]
```

``` outp
# # A tibble: 44 × 7
#    outcome term                     estimate std.error p.value     n n_group
#    <chr>   <chr>                       <dbl>     <dbl>   <dbl> <int>   <dbl>
#  1 CRBB2   bmi_grp<18.5           NA         NA         NA       495      44
#  2 CRBB2   bmi_grp18.5-25          0.0277     0.169      0.870    NA     328
#  3 CRBB2   bmi_grp25-30            0.0588     0.190      0.757    NA     119
#  4 CRBB2   bmi_grp30+              0.0811     0.548      0.883    NA       4
#  5 CRBB2   (Intercept)            -0.532      1.22       0.663    NA      NA
#  6 CRBB2   age_at_study_date      -0.00444    0.0432     0.918    NA      NA
#  7 CRBB2   I(age_at_study_date^2)  0.0000684  0.000402   0.865    NA      NA
#  8 CRBB2   sexFemale              -0.0226     0.0969     0.816    NA      NA
#  9 CRBB2   regionHarbin           -0.190      0.227      0.403    NA      NA
# 10 CRBB2   regionHaikou           -0.323      0.259      0.214    NA      NA
# # ℹ 34 more rows
```

## Analysis of multiple participant charateristics

Multiple models with different exposures can be specified by setting
`exposures` as a character vector.

Use the `remove_adjust_for` argument to remove terms from the adjustment
for a specific exposure (or outcome). Here, the age terms are removed
from adjustment when the exposure is `age_at_study_date`.

``` r

results <- multiple_models(
  olink_analysis_data,
  model = "linear",
  exposures = c("age_at_study_date",
                "sbp_mean",
                "has_diabetes",
                "current_regular_smoker"),
  outcomes = "ACAA1",
  adjust_for = list(basic = basic_adjustment_terms),
  remove_adjust_for = list(age_at_study_date = c("age_at_study_date",
                                                 "I(age_at_study_date^2)"))
)

results[, c("exposure", "term", "estimate", "std.error", "p.value", "n", "n_group")]
```

``` outp
# # A tibble: 171 × 7
#    exposure          term              estimate std.error p.value     n n_group
#    <chr>             <chr>                <dbl>     <dbl>   <dbl> <int>   <dbl>
#  1 age_at_study_date age_at_study_date  0.00184   0.00446 0.679     496      NA
#  2 age_at_study_date (Intercept)       -1.14      0.423   0.00724    NA      NA
#  3 age_at_study_date sexFemale         -0.0670    0.0960  0.485      NA      NA
#  4 age_at_study_date regionHarbin       0.197     0.222   0.374      NA      NA
#  5 age_at_study_date regionHaikou       0.237     0.254   0.351      NA      NA
#  6 age_at_study_date regionSuzhou       0.447     0.235   0.0571     NA      NA
#  7 age_at_study_date regionLiuzhou      0.636     0.241   0.00858    NA      NA
#  8 age_at_study_date regionSichuan      0.405     0.225   0.0721     NA      NA
#  9 age_at_study_date regionGansu        0.102     0.238   0.670      NA      NA
# 10 age_at_study_date regionHenan        0.515     0.239   0.0319     NA      NA
# # ℹ 161 more rows
```

Note: `has_diabetes` is numeric in the data so is treated as a
continuous variable and reported in a single row. Create a factor to
treat it is a categorical exposure.

``` r

olink_analysis_data$has_diabetes <- factor(olink_analysis_data$has_diabetes,
                                           levels = c(0, 1),
                                           labels = c("No", "Yes"))

results <- multiple_models(
  olink_analysis_data,
  model = "linear",
  exposures = c("age_at_study_date",
                "sbp_mean",
                "has_diabetes",
                "current_regular_smoker"),
  outcomes = "ACAA1",
  adjust_for = list(basic = basic_adjustment_terms),
  remove_adjust_for = list(age_at_study_date = c("age_at_study_date",
                                                 "I(age_at_study_date^2)"))
)

results[, c("exposure", "term", "estimate", "std.error", "p.value", "n", "n_group")]
```

``` outp
# # A tibble: 172 × 7
#    exposure          term              estimate std.error p.value     n n_group
#    <chr>             <chr>                <dbl>     <dbl>   <dbl> <int>   <dbl>
#  1 age_at_study_date age_at_study_date  0.00184   0.00446 0.679     496      NA
#  2 age_at_study_date (Intercept)       -1.14      0.423   0.00724    NA      NA
#  3 age_at_study_date sexFemale         -0.0670    0.0960  0.485      NA      NA
#  4 age_at_study_date regionHarbin       0.197     0.222   0.374      NA      NA
#  5 age_at_study_date regionHaikou       0.237     0.254   0.351      NA      NA
#  6 age_at_study_date regionSuzhou       0.447     0.235   0.0571     NA      NA
#  7 age_at_study_date regionLiuzhou      0.636     0.241   0.00858    NA      NA
#  8 age_at_study_date regionSichuan      0.405     0.225   0.0721     NA      NA
#  9 age_at_study_date regionGansu        0.102     0.238   0.670      NA      NA
# 10 age_at_study_date regionHenan        0.515     0.239   0.0319     NA      NA
# # ℹ 162 more rows
```

## Analysis of multiple proteins

Multiple models with different outcomes can be specified by setting
`outcomes` as a character vector.

``` r

somascan_results <- multiple_models(
  somascan_analysis_data,
  model = "linear",
  exposures = c("sbp_mean", "has_diabetes"),
  outcomes = somascan_meta_data$name,
  adjust_for = list(basic = basic_adjustment_terms)
)
```

## Full analyses

Multiple levels of adjustment are specified using a named list.

For comparison of nested models the adjustment variables of the simpler
model should be included as adjustment factors in the more complex
model.

``` r

basic_lifestyle_adjustment_terms <- c(basic_adjustment_terms,
                                      "smoking_3_groups",
                                      "alcohol_3_groups")

olink_results <- multiple_models(
  olink_analysis_data,
  model = "linear",
  exposures = c("sbp_mean", "has_diabetes"),
  outcomes = olink_meta_data$name,
  adjust_for = list(basic = basic_adjustment_terms,
                    basic_lifestyle = basic_lifestyle_adjustment_terms)
)
```

## Check for errors and warnings

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

## Select one result per protein

To summarise the somascan results as in the referenced papers, one
aptamer per Uniprot ID was selected choosing the result with the
smallest p-value and no error or warning (or just the smallest p-value
if all have error or warning).

``` r

somascan_results <- somascan_results %>%
  dplyr::filter(!is.na(estimate))

somascan_results_per_protein <- somascan_results %>%
  dplyr::inner_join(somascan_meta_data, by = c("outcome" = "name")) %>%
  dplyr::group_by(adjustment, exposure, uniprot) %>%
  dplyr::filter(p.value == min(p.value)) %>%
  dplyr::distinct() %>%
  dplyr::ungroup()
```

In the case of Olink there are a number of assays that are associated
with more than one panel. These are identified in the table below and

To summarise the Olink results as in the referenced papers, where more
than one panel was associated with one assay a single assay was selected
across all 4 panels choosing the result with the smallest p-value and no
error or warning (or just the smallest p-value if all have error or
warning)

``` r

olink_results_per_protein <- olink_results %>%
  dplyr::filter(!is.na(estimate)) %>%
  dplyr::inner_join(olink_meta_data, by = c("outcome" = "name")) %>%
  dplyr::group_by(adjustment, exposure, uniprot) %>%
  dplyr::filter(p.value == min(p.value)) %>%
  dplyr::distinct() %>%
  dplyr::ungroup()
```
