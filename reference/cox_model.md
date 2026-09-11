# Estimate the association between an exposure and incident disease

This function is used to estimate the (adjusted) association between an
exposure variable and an incident disease endpoint using a Cox
proportional hazards model. It is a wrapper for the
[coxph](https://rdrr.io/pkg/survival/man/coxph.html) function for
convenience and provides results from some additional analysis steps if
requested.

## Usage

``` r
cox_model(
  data,
  exposure,
  outcome,
  outcome_date = paste0(outcome, "_date"),
  start_date = "study_date",
  time_in = NULL,
  time_out = NULL,
  adjust_for = NULL,
  exclude = NULL,
  exposure_term = exposure,
  cluster = NULL,
  init_values = c(0),
  additional = "",
  return = c("model", "reduced-model", "coefs", "exposure"),
  return_fn = identity,
  ...
)
```

## Arguments

- data:

  A data frame.

- exposure:

  Name of column containing exposure variable.

- outcome:

  Name of the column containing status indicator (normally 0 = no
  disease, 1 = disease).

- outcome_date:

  Name of column containing date of event or censoring date. (Default:
  `paste0(outcome, "_date)`)

- start_date:

  Name of column containing date of start of follow up. (Default:
  "study_date")

- time_in:

  Name of column containing starting time for interval data.

- time_out:

  Name of column containing follow up time (or ending time for interval
  data).

- adjust_for:

  Character vector of terms to include in the model formula.

- exclude:

  Names of columns that identify observations that should be excluded.

- exposure_term:

  The term(s) to include in the model for the exposure. This can be a
  character vector or a function that take the exposure name as the
  argument and returns a character vector. (Default: exposure)

- cluster:

  Name of column that clusters the observations, for the purposes of a
  robust variance.

- init_values:

  A named vector of initial values for model coefficients.

- additional:

  Extra results to include. See "Additional results" for details.

- return:

  What should the function return:

  - "model" (default): The model object (see ?survival::coxph.object).

  - "reduced-model": The model object with some large components
    removed.

  - "coefs": A data frame of all estimated coefficients and other key
    model statistics.

  - "exposure": A data frame of estimated coefficient(s) for the
    exposure and other key model statistics.

- return_fn:

  If return="model", name of function to apply to model before it's
  returned. (Default: identity)

- ...:

  Arguments passed to
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)

## Value

Either a model object, reduced model object or data frame of
coefficients.

## Time scale

For analysis of incident disease using time-on-study as the time scale,
`outcome`, `outcome_date` and `start_date` are used to define follow up
time (right censored data).

For other anlyses, the `time_in` and `time_out` arguments can be used.
If `time_in` and `time_out` are given, then these are starting and
ending time for interval survival data. If `time_out` (but not
`time_in`) is given, then this is the follow up time (right censored
data).

## Additional results

The `additional` argument is used to add results from additional
analyses. It must be a character vector containing one or more of:

- "lrt": Likehood ratio test comparing the full model to one without the
  exposure. Output will have additional columns "Chisq" (chi-square test
  statistics), "Df" (degrees of freedom), "p_Chisq" (p-value).

- "float": Standard errors for the exposure calculated by
  [float](https://rdrr.io/pkg/Epi/man/float.html). Output will have
  additional columns "float_estimate", "float_se" and "limits".

- "qvcalc": Standard errors for the exposure calculated by
  [qvcalc](https://davidfirth.github.io/qvcalc/reference/qvcalc.html).
  Output will have additional columns "qvcalc_estimate", "quasi_se" and
  "worstErrors".

- "" (default): Nothing.

These will only be included when relevant e.g. standard errors from
[float](https://rdrr.io/pkg/Epi/man/float.html) when the exposure has
more than two levels. This only has an effect if `return` is "coefs" or
"exposure".

See [float](https://rdrr.io/pkg/Epi/man/float.html) and
[qvcalc](https://davidfirth.github.io/qvcalc/reference/qvcalc.html) for
details about these results. Note that "limits" and "worstErrors" do not
correspond to particular groups but accuracy other all contrasts.

## Examples

``` r
cox_model(ckbtools_participant_data,
          exposure = "bmi_calc",
          outcome = "all_cause_mortality",
          adjust_for = c("age_at_study_date", "strata(sex)"),
          return = "coefs")
#> # A tibble: 2 × 8
#>   term              estimate std.error statistic p.value     n nevent formula   
#>   <chr>                <dbl>     <dbl>     <dbl>   <dbl> <int>  <int> <chr>     
#> 1 bmi_calc          -0.00269   0.00588    -0.458   0.647 25000   2513 Surv(stud…
#> 2 age_at_study_date  0.101     0.00208    48.4     0        NA     NA NA        

cox_model(ckbtools_participant_data,
          exposure = "bmi_grp",
          outcome = "all_cause_mortality",
          adjust_for = c("age_at_study_date", "strata(sex)"),
          return = "coefs")
#> # A tibble: 5 × 10
#>   term    estimate std.error statistic p.value     n nevent n_group nevent_group
#>   <chr>      <dbl>     <dbl>     <dbl>   <dbl> <int>  <int>   <dbl>        <dbl>
#> 1 bmi_gr… NA        NA         NA       NA     25000   2513    2154          261
#> 2 bmi_gr… -0.00641   0.0669    -0.0958   0.924    NA     NA   15903         1597
#> 3 bmi_gr… -0.0135    0.0744    -0.181    0.856    NA     NA    6513          601
#> 4 bmi_gr…  0.214     0.150      1.43     0.152    NA     NA     430           54
#> 5 age_at…  0.101     0.00208   48.4      0        NA     NA      NA           NA
#> # ℹ 1 more variable: formula <chr>
```
