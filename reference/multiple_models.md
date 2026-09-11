# Multiple linear, logistic or Cox proportional hazards regression models

Multiple linear, logistic or Cox proportional hazards regression models

## Usage

``` r
multiple_models(
  data,
  model,
  exposures,
  outcomes,
  outcome_dates = paste0(outcomes, "_date"),
  start_date = "study_date",
  time_in = NULL,
  time_out = NULL,
  age_at_risk = NULL,
  by = NULL,
  exclude = NULL,
  exclude_controls = NULL,
  exclude_by_exposure = NULL,
  exclude_by_outcome = NULL,
  exclude_controls_by_outcome = NULL,
  adjust_for = list(none = NULL),
  remove_adjust_for = NULL,
  add_adjust_for = NULL,
  cluster = NULL,
  init_values = NULL,
  return = c("coefs", "exposure", "reduced-model", "model"),
  verbose = if (rlang::is_interactive()) "progress" else "updates",
  dryrun = FALSE,
  ...
)
```

## Arguments

- data:

  A data frame.

- model:

  Type of model to use:

  - "linear": Linear regression using
    [`linear_model_safely()`](https://ckbiobank.github.io/ckbtools/reference/safemodels.md).

  - "logistic": Logistic regression using
    [`logistic_model_safely()`](https://ckbiobank.github.io/ckbtools/reference/safemodels.md).

  - "cox": Cox proportional hazards model using
    [`cox_model_safely()`](https://ckbiobank.github.io/ckbtools/reference/safemodels.md).

- exposures:

  Names of columns containing exposure variables.

- outcomes:

  Names of columns containing outcome variables.

- outcome_dates:

  Names of columns containing outcome dates (for Cox model only).

- start_date:

  Name of column containing date of start of follow up. (Default:
  "study_date")

- time_in:

  Names of columns containing starting time for interval data (for Cox
  model only).

- time_out:

  Names of columns containing follow up time (or ending time for
  interval data) (for Cox model only).

- age_at_risk:

  A list of arguments for `age_at_risk()` to expand data by age-at-risk
  for each outcome before fitting models for that outcome. If this is
  used, "strata(agegrp)" is automatically added to each element of
  `adjust_for`.

- by:

  Name of column (or multiple columns) to create subgroups.

- exclude:

  Names of columns that identify observations that should be excluded.

- exclude_controls:

  Names of columns that identify observations that are to be exluded
  when outcome is equal to zero.

- exclude_by_exposure:

  A named list of character vectors of column names. When the name of an
  element matches the exposure, the column names will be added to
  `exclude`.

- exclude_by_outcome:

  A named list of character vectors of column names. When the name of an
  element matches the outcome, the column names will be added to
  `exclude`.

- exclude_controls_by_outcome:

  A named list of character vectors of column names. When the name of an
  element matches the outcome, the column names will be added to
  `exclude_controls`.

- adjust_for:

  A named list of character vectors of terms to include in model
  formulas.

- remove_adjust_for:

  A named list of character vectors of terms. When the name of an
  element matches either the outcome or exposure, terms given in the
  character vector will be removed from every element of `adjust_for`.

- add_adjust_for:

  A named list of character vectors of terms. When the name of an
  element matches either the outcome or exposure, terms given in the
  character vector will be added to every element of `adjust_for`.

- cluster:

  Name of column that clusters the observations, for the purposes of a
  robust variance.

- init_values:

  For Cox models, initial values for model coefficients. This must be a
  data frame with columns "outcome", "subgroup", "adjustment", "term"
  and "estimate", where term and estimate are the coefficient names and
  initial values.

- return:

  What should the model function return:

  - "exposure" (default): A data frame of estimated coefficient(s) for
    the exposure and other key model statistics.

  - "coefs": A data frame of all estimated coefficients and other key
    model statistics.

  - "model": The model object.

  - "reduced-model": The model object with some large components
    removed.

  Note that R model objects can be very large (they can contain the
  model data, residuals, fitted values, etc.) so it is not recommended
  to set `return="model"` for a large number of models.

- verbose:

  "progress" (default) to show progress bars, "updates" for messages
  without progress bars, or "quiet" to hide.

- dryrun:

  Return only the data frame of models to be run. Columns (except
  'subgroup', 'adjustment' and 'outcome_name') are the arguments that
  will be used in the model function. (Default: FALSE)

- ...:

  Arguments passed to the model function.

## Value

A "tibble" data frame. If `return` is "model" or "reduced-model" then it
will contain a list-column called "model".

## Iteration

This function provides a convenient way to iterate over models using
[`linear_model()`](https://ckbiobank.github.io/ckbtools/reference/linear_model.md),
[`logistic_model()`](https://ckbiobank.github.io/ckbtools/reference/logistic_model.md)
or
[`cox_model()`](https://ckbiobank.github.io/ckbtools/reference/cox_model.md).
First a data frame defining all models is created. If the number of
exposures is greater than or equal to the number of outcomes (or Cox
models with age-at-risk expansion are used) then the models are split by
outcome, otherwise they are split by exposure. Within each set of
models,
[`purrr::pmap()`](https://purrr.tidyverse.org/reference/pmap.html) is
used to iterate over the fitting of models. Results are then combined
into a single data frame.

For iteration of other models or analyses, you must write your own code.
The Iteration chapter of the first edition of “R for Data Science”
provides an introduction to iteration using R:
<https://r4ds.had.co.nz/iteration.html>

## See also

[`linear_model()`](https://ckbiobank.github.io/ckbtools/reference/linear_model.md),
[`logistic_model()`](https://ckbiobank.github.io/ckbtools/reference/logistic_model.md),
[`cox_model()`](https://ckbiobank.github.io/ckbtools/reference/cox_model.md),
[`linear_model_safely()`](https://ckbiobank.github.io/ckbtools/reference/safemodels.md),
[`logistic_model_safely()`](https://ckbiobank.github.io/ckbtools/reference/safemodels.md),
[`cox_model_safely()`](https://ckbiobank.github.io/ckbtools/reference/safemodels.md)

## Examples

``` r

# linear regression models
multiple_models(
  ckbtools_participant_data,
  model = "linear",
  exposures = c("age_at_study_date", "sex", "bmi_calc", "bmi_grp"),
  outcomes = c("sbp_mean")
)
#> ✔ Fitting 4 models for sbp_mean 
#> # A tibble: 12 × 17
#>    model_id outcome  exposure subgroup adjustment error warnings messages output
#>       <int> <chr>    <chr>    <chr>    <chr>      <chr> <chr>    <chr>    <chr> 
#>  1        1 sbp_mean age_at_… all      none       NA    NA       NA       NA    
#>  2        1 sbp_mean age_at_… all      none       NA    NA       NA       NA    
#>  3        2 sbp_mean sex      all      none       NA    NA       NA       NA    
#>  4        2 sbp_mean sex      all      none       NA    NA       NA       NA    
#>  5        2 sbp_mean sex      all      none       NA    NA       NA       NA    
#>  6        3 sbp_mean bmi_calc all      none       NA    NA       NA       NA    
#>  7        3 sbp_mean bmi_calc all      none       NA    NA       NA       NA    
#>  8        4 sbp_mean bmi_grp  all      none       NA    NA       NA       NA    
#>  9        4 sbp_mean bmi_grp  all      none       NA    NA       NA       NA    
#> 10        4 sbp_mean bmi_grp  all      none       NA    NA       NA       NA    
#> 11        4 sbp_mean bmi_grp  all      none       NA    NA       NA       NA    
#> 12        4 sbp_mean bmi_grp  all      none       NA    NA       NA       NA    
#> # ℹ 8 more variables: term <chr>, estimate <dbl>, std.error <dbl>,
#> #   statistic <dbl>, p.value <dbl>, n <int>, formula <chr>, n_group <dbl>

# logistic regression models
multiple_models(
  ckbtools_participant_data,
  model = "logistic",
  exposures = c("age_at_study_date", "sex", "bmi_calc", "bmi_grp"),
  outcomes = c("has_diabetes")
)
#> ✔ Fitting 4 models for has_diabetes 
#> # A tibble: 12 × 21
#>    model_id outcome  exposure subgroup adjustment error warnings messages output
#>       <int> <chr>    <chr>    <chr>    <chr>      <chr> <chr>    <chr>    <chr> 
#>  1        1 has_dia… age_at_… all      none       NA    NA       NA       NA    
#>  2        1 has_dia… age_at_… all      none       NA    NA       NA       NA    
#>  3        2 has_dia… sex      all      none       NA    NA       NA       NA    
#>  4        2 has_dia… sex      all      none       NA    NA       NA       NA    
#>  5        2 has_dia… sex      all      none       NA    NA       NA       NA    
#>  6        3 has_dia… bmi_calc all      none       NA    NA       NA       NA    
#>  7        3 has_dia… bmi_calc all      none       NA    NA       NA       NA    
#>  8        4 has_dia… bmi_grp  all      none       NA    NA       NA       NA    
#>  9        4 has_dia… bmi_grp  all      none       NA    NA       NA       NA    
#> 10        4 has_dia… bmi_grp  all      none       NA    NA       NA       NA    
#> 11        4 has_dia… bmi_grp  all      none       NA    NA       NA       NA    
#> 12        4 has_dia… bmi_grp  all      none       NA    NA       NA       NA    
#> # ℹ 12 more variables: term <chr>, estimate <dbl>, std.error <dbl>,
#> #   statistic <dbl>, p.value <dbl>, n <int>, ncases <int>, ncontrols <int>,
#> #   formula <chr>, n_group <dbl>, ncases_group <dbl>, ncontrols_group <dbl>

# Cox proportional hazards models
multiple_models(
  ckbtools_participant_data,
  model = "cox",
  exposures = c("age_at_study_date", "sex", "bmi_calc", "bmi_grp"),
  outcomes = c("all_cause_mortality", "ihd")
)
#> ✔ Fitting 4 models for all_cause_mortality 
#> ✔ Fitting 4 models for ihd 
#> # A tibble: 16 × 19
#>    model_id outcome  exposure subgroup adjustment error warnings messages output
#>       <int> <chr>    <chr>    <chr>    <chr>      <chr> <chr>    <chr>    <chr> 
#>  1        1 all_cau… age_at_… all      none       NA    NA       NA       NA    
#>  2        2 ihd      age_at_… all      none       NA    NA       NA       NA    
#>  3        3 all_cau… sex      all      none       NA    NA       NA       NA    
#>  4        3 all_cau… sex      all      none       NA    NA       NA       NA    
#>  5        4 ihd      sex      all      none       NA    NA       NA       NA    
#>  6        4 ihd      sex      all      none       NA    NA       NA       NA    
#>  7        5 all_cau… bmi_calc all      none       NA    NA       NA       NA    
#>  8        6 ihd      bmi_calc all      none       NA    NA       NA       NA    
#>  9        7 all_cau… bmi_grp  all      none       NA    NA       NA       NA    
#> 10        7 all_cau… bmi_grp  all      none       NA    NA       NA       NA    
#> 11        7 all_cau… bmi_grp  all      none       NA    NA       NA       NA    
#> 12        7 all_cau… bmi_grp  all      none       NA    NA       NA       NA    
#> 13        8 ihd      bmi_grp  all      none       NA    NA       NA       NA    
#> 14        8 ihd      bmi_grp  all      none       NA    NA       NA       NA    
#> 15        8 ihd      bmi_grp  all      none       NA    NA       NA       NA    
#> 16        8 ihd      bmi_grp  all      none       NA    NA       NA       NA    
#> # ℹ 10 more variables: term <chr>, estimate <dbl>, std.error <dbl>,
#> #   statistic <dbl>, p.value <dbl>, n <int>, nevent <int>, formula <chr>,
#> #   n_group <dbl>, nevent_group <dbl>

```
