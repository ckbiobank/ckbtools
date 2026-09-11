# Variance of the log risk

## Introduction

It is often necessary to compare estimates for a categorical variable
with more than two categories. In order to make these comparisons, it is
necessary to select an appropriate reference category; this can often be
challenging as, in some cases, the choice of reference category is not
obvious and can be rather arbitrary.

When there is a categorical variable with more than two levels in a
regression model, the estimates and corresponding standard errors only
allow a direct comparison between each specified category relative to
the category selected as the reference group. To address this
limitation, *‘floating’ standard errors*, or standard errors based on
*quasi-variances*, may be used. These are also called *floating absolute
risks*. This method allows for the estimation of standard errors across
all levels of the categorical variable, including the reference group,
enabling comparisons between any pair of groups.

In the context of floating absolute risks, various algorithms are
available which use slightly different methods for estimation of the
variance-covariance matrix. The variance-covariance matrix captures the
variability and relationship between estimates across levels of a
categorical variable. This enables consistent standard error estimation
across all groups, thereby facilitating accurate comparisons between any
pair of groups by accounting for the relative differences in their
standard errors (Easton et al, 1991; Plummer, 2004; Firth and de
Menezes, 2004).

In R there are two functions in two packages that implement this:
[`Epi::float`](https://rdrr.io/pkg/Epi/man/float.html) uses the method
of Plummer (2004) and
[`qvcalc::qvcalc`](https://davidfirth.github.io/qvcalc/reference/qvcalc.html)
the method of Firth and de Menezes (2004).

Both functions are re-exported in the `ckbtools` package, so they can be
used after calling
[`library(ckbtools)`](https://ckbiobank.github.io/ckbtools).

## `float`

Usage:

``` r

float(object)
```

where `object` is any type of regression model.

### Example

Cox regression of `ihd` on categories of `bmi_calc`.

``` r

library(survival)
ckbtools_participant_data$bmi_grp <- cut(ckbtools_participant_data$bmi_calc, c(0, 18.5, 25, 30, 100), right = FALSE)
ckbtools_participant_data$ihd_time <- as.numeric(ckbtools_participant_data$ihd_date - ckbtools_participant_data$study_date) / 365.25

cox_model <- coxph(Surv(ihd_time, ihd) ~ bmi_grp +
                     age_at_study_date + is_female + as.factor(region_code),
                   data = ckbtools_participant_data)
fl <- float(cox_model)
data.frame(coef = fl$coef, var = fl$var)
```

``` outp
#                coef          var
# [0,18.5)  0.0000000 0.0074317718
# [18.5,25) 0.4370379 0.0007175254
# [25,30)   0.6873568 0.0015079493
# [30,100)  1.2335243 0.0135786368
```

## `qvcalc`

Usage:

``` r

qvcalc(object)
```

where `object` is any type of regression model.

### Example

Cox regression of `ihd` on categories of `bmi_calc`.

``` r

qvcalc(cox_model, factorname = "bmi_grp")
```

You can also include the `coef.indices` argument to determine which rows
and columns of the variance-covariance matrix to use (i.e. which factors
to consider in the calculation). Add a zero to include the reference
level. This will enable the estimation of the reference group’s
‘floating’ standard error, which allows the calculation of confidence
intervals. This is particularly useful for creating display items such
as shape plots of hazard ratios/risk ratios etc. across all categories,
as it allows for the inclusion of the reference group estimate and
confidence intervals.

``` r

qvcalc(cox_model, coef.indices = c(0, 1, 2, 3))
```

## References

[Plummer, M. Improved estimates of floating absolute risk. Statistics in
Medicine 23, 93–104 (2004).](https://doi.org/10.1002/sim.1485) (Not open
access)

[Firth, D. & De Menezes, R. X. Quasi‐variances. Biometrika 91, 65–80
(2004)](https://doi.org/10.1093/biomet/91.1.65)

[Firth, D. Overcoming the Reference Category Problem in the Presentation
of Statistical Models. Sociological Methodology 33, 1–18
(2003).](https://doi.org/10.1111/j.0081-1750.2003.t01-1-00125.x) (Not
open access)

[Firth, D. Quasi-variances in Xlisp-Stat and on the web. Journal of
Statistical Software 5, 1–13
(2000).](https://doi.org/10.18637/jss.v005.i04)

[Easton, D. F., Peto, J. & Babiker, A. G. Floating absolute risk: An
alternative to relative risk in survival and case-control analysis
avoiding an arbitrary reference group. Statistics in Medicine 10,
1025–1035 (1991).](https://doi.org/10.1002/sim.4780100703) (Not open
access)
