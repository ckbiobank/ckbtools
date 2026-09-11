# Estimate the association between an exposure and a continuous outcome

This function is used to estimate the (adjusted) linear association
between an exposure variable and a continuous outcome variable using a
linear regression model. It is a wrapper for the
[lm](https://rdrr.io/r/stats/lm.html) function for convenience and
provides results from some additional analysis steps if requested.

## Usage

``` r
linear_model(
  data,
  exposure,
  outcome,
  adjust_for = NULL,
  exclude = NULL,
  exposure_term = exposure,
  additional = "",
  return = c("model", "reduced-model", "coefs", "exposure"),
  return_fn = identity,
  nuisance = character(0),
  ...
)
```

## Arguments

- data:

  A data frame.

- exposure:

  Name of column containing exposure variable.

- outcome:

  Name of column containing outcome variable.

- adjust_for:

  Character vector of terms to include in the model formula.

- exclude:

  Names of columns that identify observations that should be excluded.

- exposure_term:

  The term(s) to include in the model for the exposure. This can be a
  character vector or a function that take the exposure name as the
  argument and returns a character vector. (Default: exposure)

- additional:

  Extra results to include. See "Additional results" for details.

- return:

  What should the function return:

  - "model" (default): The model object (object of class "lm").

  - "reduced-model": The model object with some large components
    removed.

  - "coefs": A data frame of all estimated coefficients and other key
    model statistics.

  - "exposure": A data frame of estimated coefficient(s) for the
    exposure and other key model statistics.

- return_fn:

  If return="model", name of function to apply to model before it's
  returned. (Default: identity)

- nuisance:

  If "emmeans" is requested via the `additional` argument, use this to
  specify predictors to omit from the reference grid. See the help for
  [emmeans](https://rvlenth.github.io/emmeans/reference/emmeans.html)
  and
  [ref_grid](https://rvlenth.github.io/emmeans/reference/ref_grid.html).

- ...:

  Arguments passed to [`stats::lm()`](https://rdrr.io/r/stats/lm.html)

## Value

Either a model object, reduced model object or data frame of
coefficients.

## Additional results

The `additional` argument is used to add results from additional
analyses. It must be a character vector containing one or more of:

- "anova": ANOVA F-test comparing the full model to one without the
  exposure. Output will have additional columns "F" (F test statistics),
  "Df" (degrees of freedom), "p_F" (p-value).

- "float": Standard errors for the exposure calculated by
  [float](https://rdrr.io/pkg/Epi/man/float.html). Output will have
  additional columns "float_estimate", "float_se" and "limits".

- "qvcalc": Standard errors for the exposure calculated by
  [qvcalc](https://davidfirth.github.io/qvcalc/reference/qvcalc.html).
  Output will have additional columns "qvcalc_estimate", "quasi_se" and
  "worstErrors".

- "emmeans": Estimated marginal means for the exposure calculated using
  [emmeans](https://rvlenth.github.io/emmeans/reference/emmeans.html)
  with proportional weights for factor variables that are averaged over.
  Output will have additional columns "emmean", "emmean_se" and
  "emmean_quasi_se".

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
linear_model(ckbtools_participant_data,
             exposure = "age_at_study_date",
             outcome = "sbp_mean",
             adjust_for = "sex",
             return = "coefs")
#> # A tibble: 3 × 7
#>   term              estimate std.error statistic  p.value     n formula         
#>   <chr>                <dbl>     <dbl>     <dbl>    <dbl> <int> <chr>           
#> 1 age_at_study_date    0.125   0.00945      13.2 1.42e-39 25000 sbp_mean ~ age_…
#> 2 (Intercept)        102.      0.516       198.  0           NA NA              
#> 3 sexFemale           -2.74    0.205       -13.4 7.53e-41    NA NA              

# factor exposure, calculate estimated marginal means
linear_model(ckbtools_participant_data,
             exposure = "smoking",
             outcome = "sbp_mean",
             adjust_for = c("age_at_study_date", "sex"),
             additional = "emmeans",
             return = "coefs")
#> # A tibble: 7 × 11
#>   term   estimate std.error statistic   p.value emmean emmean_se emmean_quasi_se
#>   <chr>     <dbl>     <dbl>     <dbl>     <dbl>  <dbl>     <dbl>           <dbl>
#> 1 smoki…   NA      NA          NA     NA          108.     0.157           0.197
#> 2 smoki…   -0.673   0.352      -1.91   5.61e- 2   107.     0.305           0.296
#> 3 smoki…   -0.297   0.352      -0.845  3.98e- 1   108.     0.304           0.295
#> 4 smoki…   -1.98    0.299      -6.64   3.15e-11   106.     0.224           0.220
#> 5 (Inte…  103.      0.555     186.     0           NA     NA              NA    
#> 6 age_a…    0.125   0.00944    13.2    9.49e-40    NA     NA              NA    
#> 7 sexFe…   -3.66    0.256     -14.3    2.55e-46    NA     NA              NA    
#> # ℹ 3 more variables: n <int>, n_group <dbl>, formula <chr>

# return the model object
linear_model(ckbtools_participant_data,
             exposure = "bmi_grp",
             outcome = "sbp_mean",
             adjust_for = c("age_at_study_date", "sex"),
             return = "model")
#> 
#> Call:
#> stats::lm(formula = model_formula, data = data)
#> 
#> Coefficients:
#>       (Intercept)     bmi_grp18.5-25       bmi_grp25-30         bmi_grp30+  
#>           95.9089             4.7895             9.7153            15.7819  
#> age_at_study_date          sexFemale  
#>            0.1386            -2.9989  
#> 
```
