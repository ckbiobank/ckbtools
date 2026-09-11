# Associations of proteins with incident disease

This article demonstrates how to use this package to estimate
associations between proteins and risk of incidence disease using Cox
proportional hazards models. The methods used are from the paper
[“Proteome-wide assessment of ~9500 plasma proteins with major diseases
and multi-morbidity in a Chinese
population”](https://doi.org/10.21203/rs.3.rs-5356315/v1).

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

## Prepare participant data

### Baseline and endpoint data

A data set containing relevant participant baseline and endpoint data
must be prepared manually. For these analyses, the data must contain
columns for disease endpoints (status and date) and all participant
characteristics used for exclusions or as adjustment variables. See [CKB
data](https://ckbiobank.github.io/ckbtools/articles/ckb_data.md) for an
example.

### Select subcohort participants

These analyses of incident disease will use only participants in the
subcohort of the nested case-subcohort proteomics study.

``` r

subcohort_participant_data <- ckbtools_participant_data[ckbtools_participant_data$proteomics_subcohort == 1, ]
subcohort_participant_data$subcohort <- 1
```

### Participant exclusions

Where exclusions are needed for particular analyses, create columns that
identify addition participants to be excluded. These will be used later
when fitting the survival models. For example, create columns that
identify participants to be excluded in analyses of cancer (if they have
prior cancer diagnosis at baseline) or analyses of respiratory diseases
(if they have prior TB, asthma or COPD at baseline).

``` r

subcohort_participant_data$exclude_for_cancer <- subcohort_participant_data$cancer_diag == "Yes"
subcohort_participant_data$exclude_for_respiratory <- subcohort_participant_data$tb_diag == "Yes" |
  subcohort_participant_data$asthma_diag == "Yes" |
  subcohort_participant_data$has_copd == "Yes"
```

## Prepare protein data

Prepare the protein data as described in
[Proteomics](https://ckbiobank.github.io/ckbtools/articles/proteomics.html#preparing-proteomics-data).

For this analysis, Olink values are centred using the plat medians
(plate correction “median-subcohort”) and SomaScan values are log₂
transformed. Protein values are then scaled to have standard deviations
of one.

``` r

## set CKB data release
options(ckb.data.release = "19.02")

## meta data
olink_meta_data <- prepare_olink_meta_data(ckbtools_raw_olink_meta_data,
                                           ckbtools_raw_olink_protein_data,
                                           names = "id")

somascan_meta_data <- prepare_somascan_meta_data(ckbtools_raw_somascan_meta_data,
                                                 names = "id")

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

``` r

## proteomics analysis data
olink_analysis_data <- prepare_proteomics_analysis_dataset(
  subcohort_participant_data,
  olink_protein_data,
  olink_meta_data,
  plate_correction = "median-subcohort",
  scaling = "scale"
)
```

``` outp
# ✔ Plate correction 'median-subcohort' has been applied to protein values.
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

## Analysis of one incident disease and protein

The
[`cox_model()`](https://ckbiobank.github.io/ckbtools/reference/cox_model.md)
function is used to fit one Cox proportional hazards model. The code
below will return the fitted proportional hazard models object for an
analysis of protein OID13950 with all-cause mortality, with
time-on-study as the time scale. The model has baseline hazard
stratified by sex and region, and is adjusted for age at baseline and
its square, hours since last ate and its square, temperature and its
square, and education.

``` r

cox_model(olink_analysis_data,
          exposure = "OID13950",
          outcome = "all_cause_mortality",
          adjust_for = c("strata(sex)",
                         "strata(region)",
                         "poly(age_at_study_date, 2)",
                         "poly(hours_since_last_ate, 2)",
                         "poly(region_mean_temp, 2)",
                         "education_4_groups"),
          exclude = NULL,
          return = "model")
```

``` outp
# Call:
# survival::coxph(formula = model_formula, data = data, init = init_values_aligned, 
#     y = FALSE, cluster = if (is.null(cluster)) NULL else eval(as.name(cluster)), 
#     iter.max = 40)
# 
#                                               coef  exp(coef)   se(coef)      z
# OID13950                                -1.074e-01  8.981e-01  1.533e-01 -0.701
# poly(age_at_study_date, 2)1              2.223e+01  4.512e+09  5.148e+00  4.318
# poly(age_at_study_date, 2)2              4.338e+00  7.652e+01  3.833e+00  1.132
# poly(hours_since_last_ate, 2)1          -6.704e-01  5.115e-01  3.580e+00 -0.187
# poly(hours_since_last_ate, 2)2          -1.235e+00  2.908e-01  3.544e+00 -0.349
# poly(region_mean_temp, 2)1              -3.427e+00  3.247e-02  3.908e+00 -0.877
# poly(region_mean_temp, 2)2              -9.686e+00  6.215e-05  4.422e+00 -2.190
# education_4_groupsPrimary School         6.184e-01  1.856e+00  5.180e-01  1.194
# education_4_groupsMiddle School          8.075e-01  2.242e+00  5.479e-01  1.474
# education_4_groupsHigh School or higher  1.461e-02  1.015e+00  5.386e-01  0.027
#                                                p
# OID13950                                  0.4833
# poly(age_at_study_date, 2)1             1.57e-05
# poly(age_at_study_date, 2)2               0.2578
# poly(hours_since_last_ate, 2)1            0.8514
# poly(hours_since_last_ate, 2)2            0.7275
# poly(region_mean_temp, 2)1                0.3805
# poly(region_mean_temp, 2)2                0.0285
# education_4_groupsPrimary School          0.2325
# education_4_groupsMiddle School           0.1405
# education_4_groupsHigh School or higher   0.9784
# 
# Likelihood ratio test=51.7  on 10 df, p=1.296e-07
# n= 491, number of events= 48 
#    (9 observations deleted due to missingness)
```

*Note that education_4_groups is formatted as a factor in the data, so
it is correctly used as a categorical variable in the model. If you want
to include a categorical variable that is coded numerically then make
sure to use `"as.factor(variable)"`.*

The fitted model object can be used for further analyses, such as
inspection of model residuals, fitted values and fit statistics.

## Analysis of many incident diseases and proteins

The
[`multiple_models()`](https://ckbiobank.github.io/ckbtools/reference/multiple_models.md)
function can be used to run multiple regression models.

*It is highly recommended to run some individual analyses as described
in the previous section to check that the models and output are as
intended, before running analyses for many proteins.*

The code below fits models for all Olink and SomaScan proteins, for two
outcomes (`all_cause_mortality` and `respiratory_diseases`) and two sets
of adjustments.

For all `respiratory_diseases` models, participants with
`exclude_for_respiratory == TRUE` are excluded.

``` r

minimally_adjusted <- c("strata(sex)",
                        "strata(region)",
                        "poly(age_at_study_date, 2)",
                        "poly(hours_since_last_ate, 2)",
                        "poly(region_mean_temp, 2)",
                        "education_4_groups")

fully_adjusted <- c("strata(sex)",
                    "strata(region)",
                    "poly(age_at_study_date, 2)",
                    "poly(hours_since_last_ate, 2)",
                    "poly(region_mean_temp, 2)",
                    "education_4_groups",
                    "smoking_3_groups",
                    "alcohol_3_groups",
                    "met",
                    "sbp_mean",
                    "has_diabetes",
                    "bmi_calc")

olink_results <- multiple_models(
  olink_analysis_data,
  model = "cox",
  exposures = olink_meta_data$name,
  outcomes = c("all_cause_mortality",
               "respiratory_diseases"),
  exclude_by_outcome = list(respiratory_diseases = "exclude_for_respiratory"),
  adjust_for = list(minimal = minimally_adjusted,
                    full = fully_adjusted)
)
```

``` outp
# ✔ Fitting 96 models for all_cause_mortality
```

``` outp
# ✔ Fitting 96 models for respiratory_diseases
```

``` r

somascan_results <- multiple_models(
  somascan_analysis_data,
  model = "cox",
  exposures = somascan_meta_data$name,
  outcomes = c("all_cause_mortality",
               "respiratory_diseases"),
  exclude_by_outcome = list(respiratory_diseases = "exclude_for_respiratory"),
  adjust_for = list(minimal = minimally_adjusted,
                    full = fully_adjusted)
)
```

``` outp
# ✔ Fitting 96 models for all_cause_mortality 
# ✔ Fitting 96 models for respiratory_diseases
```

*Note that education, smoking and alcohol are formatted as factors in
the data, so they are correctly used as categorical variables in the
model. If you want to include a categorical variable that is coded
numerically then make sure to use `"as.factor(variable)"`.*

When your analysis includes a very large number of models you can split
your analyses into separate uses of the function (e.g. one for each
disease endpoint) so that each can be run and saved to a results file
before moving onto the next. This `multiple_model()` function or the
individual model functions can also be used if you wish to implement
another approach to iteration yourself.

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
as in the referenced paper, select one aptamer per Uniprot ID choosing
the result with the smallest p-value and no error or warning (or just
the smallest p-value if all have error or warning)

``` r

somascan_summary_results <- somascan_results %>%
  dplyr::inner_join(somascan_meta_data, by = c("exposure" = "name")) %>%
  dplyr::group_by(outcome, adjustment, uniprot) %>%
  dplyr::arrange(!is.na(error) | !is.na(warnings), p.value) %>%
  dplyr::filter(dplyr::row_number() == 1) %>%
  dplyr::ungroup()
```

### Multiple testing

In the referenced paper, a single multiple testing correction is applied
across Olink and SomaScan protein measures for each disease and
adjustment model. Where a result has an error or warning, the nominal
p-value is first set to one, as the result is not significant but should
still be counted as a significance test. In the referenced paper, the
false discovery rate correction of [Benjamini & Hochberg
(1995)](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x) is used.

``` r

## combine results from Olink and SomaScan analyses
all_results <- dplyr::bind_rows(olink_summary_results,
                                somascan_summary_results)  %>%
  ## handle errors and warnings
  dplyr::mutate(p.value = dplyr::if_else(!is.na(error) | !is.na(warnings), 1, p.value)) %>%
  ## multiple testing correction per outcome and adjustment
  dplyr::group_by(outcome, adjustment) %>%
  dplyr::mutate(p.value.fdr = p.adjust(p.value, method = "fdr")) %>%
  dplyr::ungroup()
```

## Present results

The number of proteins significantly associated (FDR-corrected p \<
0.05) with each disease for each adjustment model can then be
summarised.

``` r

all_results %>%
  dplyr::group_by(outcome, adjustment) %>%
  dplyr::summarise(n = dplyr::n(),
                   error_or_warning = sum(!is.na(error) | !is.na(warnings)),
                   significant = sum(p.value.fdr < 0.05),
                   significant_inverse = sum(p.value.fdr < 0.05 & estimate < 0),
                   .groups = "drop")
```

``` outp
# # A tibble: 4 × 6
#   outcome      adjustment     n error_or_warning significant significant_inverse
#   <chr>        <chr>      <int>            <int>       <int>               <int>
# 1 all_cause_m… full          96                0          51                   2
# 2 all_cause_m… minimal       96                0          49                   1
# 3 respiratory… full          96                0          50                   1
# 4 respiratory… minimal       96                0          50                   1
```
