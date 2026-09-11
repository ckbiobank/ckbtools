# Estimate the association between an exposure and a binary outcome

This function is used to estimate the (adjusted) association (log odds
ratio) between an exposure variable and a binary outcome variable using
a logistic regression model. It is a wrapper for the
[glm](https://rdrr.io/r/stats/glm.html) function with
`family = "binomial"` for convenience and provides results from some
additional analysis steps if requested.

## Usage

``` r
logistic_model(
  data,
  exposure,
  outcome,
  adjust_for = NULL,
  exclude = NULL,
  exclude_controls = NULL,
  exposure_term = exposure,
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

  Name of the column containing binary outcome variable.

- adjust_for:

  Character vector of terms to include in the model formula.

- exclude:

  Names of columns that identify observations that are to be excluded.

- exclude_controls:

  Names of columns that identify observations that are to be exluded
  when outcome is equal to zero.

- exposure_term:

  The term(s) to include in the model for the exposure. This can be a
  character vector or a function that take the exposure name as the
  argument and returns a character vector. (Default: exposure)

- additional:

  Extra results to include. See "Additional results" for details.

- return:

  What should the function return:

  - "model" (default): The model object (object of class "glm").

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

  Arguments passed to [`stats::glm()`](https://rdrr.io/r/stats/glm.html)

## Value

Either a model object, reduced model object or data frame of
coefficients.

## Additional results

The `additional` argument is used to add results from additional
analyses. It must be a character vector containing one or more of:

- "lrt": Likehood ratio test comparing the full model to one without the
  exposure. Output will have additional columns "Chisq" (chi-square test
  statistics), "Df" (degrees of freedom), "p_Chisq" (p-value).

- "float": Standard errors for the exposure calculated by
  [float](https://rdrr.io/pkg/Epi/man/float.html) or by
  [ftrend](https://rdrr.io/pkg/Epi/man/ftrend.html) if the model has no
  intercept term. Output will have additional columns "float_estimate",
  "float_se" and "limits".

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
logistic_model(ckbtools_participant_data,
               exposure = "bmi_calc",
               outcome = "has_diabetes",
               adjust_for = c("age_at_study_date", "sex"),
               return = "coefs")
#> # A tibble: 4 × 9
#>   term     estimate std.error statistic   p.value     n ncases ncontrols formula
#>   <chr>       <dbl>     <dbl>     <dbl>     <dbl> <int>  <int>     <int> <chr>  
#> 1 bmi_calc   0.109    0.00738     14.8  2.60e- 49 25000   1893     23107 has_di…
#> 2 (Interc…  -8.60     0.233      -36.9  1.20e-298    NA     NA        NA NA     
#> 3 age_at_…   0.0623   0.00232     26.9  2.14e-159    NA     NA        NA NA     
#> 4 sexFema…   0.205    0.0505       4.06 4.88e-  5    NA     NA        NA NA     

logistic_model(ckbtools_participant_data,
               exposure = "bmi_grp",
               outcome = "has_diabetes",
               adjust_for = c("age_at_study_date", "sex"),
               return = "coefs")
#> # A tibble: 7 × 12
#>   term    estimate std.error statistic    p.value     n ncases ncontrols n_group
#>   <chr>      <dbl>     <dbl>     <dbl>      <dbl> <int>  <int>     <int>   <dbl>
#> 1 bmi_gr…  NA       NA           NA    NA         25000   1893     23107    2154
#> 2 bmi_gr…   0.629    0.113        5.54  2.94e-  8    NA     NA        NA   15903
#> 3 bmi_gr…   1.09     0.117        9.32  1.15e- 20    NA     NA        NA    6513
#> 4 bmi_gr…   1.74     0.171       10.2   2.10e- 24    NA     NA        NA     430
#> 5 (Inter…  -6.77     0.179      -37.8   2.87e-312    NA     NA        NA      NA
#> 6 age_at…   0.0616   0.00231     26.6   2.58e-156    NA     NA        NA      NA
#> 7 sexFem…   0.215    0.0504       4.27  1.95e-  5    NA     NA        NA      NA
#> # ℹ 3 more variables: ncases_group <dbl>, ncontrols_group <dbl>, formula <chr>
```
