# Age-at-risk expansion

## Introduction

Analysis of time-to-event data (survival analysis) is often used
estimate associations between exposures and disease outcomes. The Cox
proportional hazards model is widely used in such analyses. This model
has the key “proportional hazards” assumption that effect of covariates
(hazard ratios) are constant over time.

Age is frequently a confounder of the association between exposure and
disease, or a covariate of interest itself. When time-on-study is used
as the time scale for a Cox model, age at baseline can be adjusted for
but this assumes a constant effect of baseline age over time and does
not account for the effect of ‘acquired’ age during follow-up.

Lexis expansion involves splitting each participant’s follow-up time
into intervals. Lexis expansion by age-at-risk will create observations
for each subject for each age band (e.g. 5-year intervals) they pass
through during follow-up. For example, using 5-year age bands between 30
years and 80 years, a participant who enters the study at age 42 years
and leaves age 53 years will have three rows of data: from baseline
until their 45th birthday; from their 45th birthday to 50th birthday;
and from their 50th birthday to the end of follow-up.

The newly created age-at-risk bands can be used as a stratification
variable. When the effect of age-at-risk is not the primary focus,
baseline hazard can be stratified by age-at-risk. Each age-at-risk band
then has its own baseline hazard. The effect of all other covariates is
assumed to be constant across age-at-risk bands and throughout
follow-up.

It may of interest to estimate an effect for a covariate in each
age-at-risk band. This can be achieved by including an interaction
between age-at-risk band and the covariate. Further, an additional
nuisance term nested within age-at-risk band can be included to avoid
assuming that proportional associations are similar across different,
wider, bands of age. See [Multi-age group](#multi-age-group).

The following examples demonstrates how to implement such models in R,
using the
[`expand_age_at_risk()`](https://ckbiobank.github.io/ckbtools/reference/expand_age_at_risk.md)
function to apply Lexis expansion for age-at-risk.

## Examples

### Expand by age-at-risk

The data frame contains one row per participant. For example:

``` r

ckbtools_participant_data[ckbtools_participant_data$csid == "1000005",
                          c("csid", "study_date", "dob_anon", "ihd", "ihd_date")]
```

``` outp
#      csid study_date   dob_anon ihd   ihd_date
# 5 1000005 2007-08-08 1942-06-15   0 2017-12-31
```

To expand our data frame we specify the age-at-risk-bands, and identify
the columns that contain participants’ unique ID (default: “csid”), date
of birth (“dob_anon”), date of study entry (“study_date”), date of study
exit (“endpoint_date”) and status at exit (“endpoint”).

The age-at-risk strata are chosen by specifying a vector of cut points.
The first will be the minimum age-at-risk to enter the analysis, and the
last will be the age at which participants leave the analysis. For
example, `seq(from = 40, to = 80, by = 5)` will use 5-year age-at-risk
strata bands from 40 to 80 years.

``` r

# Expand data frame by age-at-risk
participant_data_long <- expand_age_at_risk(data          = ckbtools_participant_data,
                                            col.exit_date = "ihd_date",
                                            col.status    = "ihd",
                                            ages          = seq(from = 40, to = 80, by = 5),
                                            col.keep      = "bmi_grp")
```

The expanded data contains multiple rows per participant: one for each
age-at-risk strata during which they have some follow-up time. For
example:

``` r

participant_data_long[participant_data_long$csid == "1000005",
                      c("csid", "agegrp",
                        "start_int", "end_int",
                        "status", "t_start_days", "t_end_days")]
```

``` outp
# # A tibble: 3 × 7
#      csid agegrp start_int  end_int    status t_start_days t_end_days
#     <dbl> <chr>  <date>     <date>      <dbl>        <dbl>      <dbl>
# 1 1000005 65-70  2007-08-08 2012-06-14      0            0      1773.
# 2 1000005 70-75  2012-06-15 2017-06-14      0         1773      3599.
# 3 1000005 75-80  2017-06-15 2017-12-31      0         3599      3799.
```

The expanded data frame contains the participant ID columns used to
expand the data (set by arguments `col.id`, `col.dob`, `col.entry_date`,
`col.exit_date`, `col.status`) and new columns:

- **agegrp_start**: Lower bound for age-at-risk strata
- **agegrp_end**: Upper bound for age-at-risk strata
- **agegrp**: Age-at-risk strata label
- **start_int**: Start date for time interval
- **end_int**: End date for time interval
- **status**: Participant status at end of age-at-risk interval
- **t_start_days**: Start time for interval (days)
- **t_end_days**: End time for interval (days)
- **time_in**: Start time for interval (years, calculated as
  days/365.25)
- **time_out**: End time for interval (years, calculated as days/365.25)

### Cox proportional hazards model stratified by age-at-risk

The expanded data frame contains column `agegrp` which can be used to
stratify the Cox proportional hazards model by including a
`strata(agegrp)` in the model formula and using the time columns to
specify the starting and ending time for intervals.

``` r

library(survival)
fit <- coxph(Surv(time_in, time_out, status) ~ bmi_grp + strata(agegrp),
             data = participant_data_long)
broom::tidy(fit)
```

``` outp
# # A tibble: 3 × 5
#   term           estimate std.error statistic  p.value
#   <chr>             <dbl>     <dbl>     <dbl>    <dbl>
# 1 bmi_grp18.5-25    0.458    0.0952      4.81 1.53e- 6
# 2 bmi_grp25-30      0.728    0.0992      7.33 2.24e-13
# 3 bmi_grp30+        1.20     0.153       7.84 4.68e-15
```

``` r

print(fit$xlevels$`strata(agegrp)`)
```

``` outp
# [1] "40-45" "45-50" "50-55" "55-60" "60-65" "65-70" "70-75" "75-80"
```

### Multi-age group

The same expanded dataset can be used for other models, such as those
estimating HRs with an interaction between the risk factor and
age-at-risk. For example, an interaction with 10-year age-at-risk band
can be included while allowing hazards to differ between the first and
second half of that band. By including an additional nuisance term
nested within age band, the HR in each 10-year band is estimated as the
geometric mean of the HRs in the first and second half of that decade,
avoiding assumptions that proportional associations are similar across
different decades of age. See [Lacey et al.,
2018](https://doi.org/10.1016/S2214-109X(18)30217-1) for an example with
major vascular disease and systolic blood pressure.

To fit such a model, first create a new column from `agegrp_start`
called `age_band` that is the 10-year age-at-risk band. Then create a
numeric column `age_band_midpoint` that is the midpoint of each decade.
Also create a column `agegrp_midpoint` that is the mid-point of each
5-year age-at-risk group. The nuisance term `nuis` is the different
between the midpoint of the 10-year age-at-risk band and the 5-year
age-at-risk group. Finally, create a column that is the interaction
between risk factor group and 10-year age-at-risk band (here,
`bmi_grp_x_age_band` for BMI group).

``` r

participant_data_long <- participant_data_long %>%
  dplyr::mutate(age_band = cut(agegrp_start,
                               breaks = c(40, 50, 60, 70, 80),
                               labels = c("40-49", "50-59", "60-69", "70-79"),
                               right = FALSE)) %>%
  dplyr::mutate(age_band_midpoint = dplyr::case_match(age_band,
                                                      "40-49" ~ 45,
                                                      "50-59" ~ 55,
                                                      "60-69" ~ 65,
                                                      "70-79" ~ 75),
                agegrp_midpoint = (agegrp_start + agegrp_end) / 2,
                nuis = age_band_midpoint - agegrp_midpoint,
                bmi_grp_x_age_band = interaction(bmi_grp, age_band, sep = "_"))
```

``` outp
# Warning: There was 1 warning in `dplyr::mutate()`.
# ℹ In argument: `age_band_midpoint = dplyr::case_match(...)`.
# Caused by warning:
# ! `case_match()` was deprecated in dplyr 1.2.0.
# ℹ Please use `recode_values()` instead.
```

In the Cox model formula, include the BMI group by 10-year age-at-risk
band interaction (`bmi_grp_x_age_band`) and the nuisance term nested in
10-year age-at-risk bands (`nuis:age_band`).

``` r

fit <- coxph(Surv(time_in, time_out, status) ~ bmi_grp_x_age_band + nuis:age_band,
             data = participant_data_long)
```

A log HR is estimated for each combination of BMI group and 10-year
age-at-risk band. The reference group is BMI \<18.5 and age-at-risk
40-49, and so this group has a hazard ratio of one.

``` r

broom::tidy(fit)
```

``` outp
# # A tibble: 19 × 5
#    term                            estimate std.error statistic  p.value
#    <chr>                              <dbl>     <dbl>     <dbl>    <dbl>
#  1 bmi_grp_x_age_band18.5-25_40-49   1.34      0.586       2.29 2.23e- 2
#  2 bmi_grp_x_age_band25-30_40-49     1.72      0.591       2.91 3.67e- 3
#  3 bmi_grp_x_age_band30+_40-49       2.11      0.707       2.98 2.88e- 3
#  4 bmi_grp_x_age_band<18.5_50-59     2.06      0.604       3.41 6.49e- 4
#  5 bmi_grp_x_age_band18.5-25_50-59   2.25      0.581       3.87 1.08e- 4
#  6 bmi_grp_x_age_band25-30_50-59     2.51      0.583       4.31 1.60e- 5
#  7 bmi_grp_x_age_band30+_50-59       2.79      0.637       4.39 1.14e- 5
#  8 bmi_grp_x_age_band<18.5_60-69     2.07      0.607       3.40 6.65e- 4
#  9 bmi_grp_x_age_band18.5-25_60-69   2.89      0.580       4.98 6.42e- 7
# 10 bmi_grp_x_age_band25-30_60-69     3.03      0.582       5.20 1.96e- 7
# 11 bmi_grp_x_age_band30+_60-69       3.68      0.613       6.01 1.89e- 9
# 12 bmi_grp_x_age_band<18.5_70-79     3.27      0.593       5.51 3.61e- 8
# 13 bmi_grp_x_age_band18.5-25_70-79   3.54      0.580       6.11 1.01e- 9
# 14 bmi_grp_x_age_band25-30_70-79     3.91      0.582       6.72 1.85e-11
# 15 bmi_grp_x_age_band30+_70-79       4.36      0.616       7.08 1.42e-12
# 16 nuis:age_band40-49               -0.0374    0.0318     -1.18 2.40e- 1
# 17 nuis:age_band50-59               -0.0622    0.0176     -3.53 4.18e- 4
# 18 nuis:age_band60-69               -0.0827    0.0153     -5.42 6.06e- 8
# 19 nuis:age_band70-79               -0.0719    0.0151     -4.75 1.99e- 6
```

## Note

The end times for each age-at-risk strata (columns t_end_days and
time_out) are always x.95. This helps avoid ambiguity about risk sets,
but the time between the start and end times does not equal a whole
number of days. So these columns should not be used to calculate rates
or time at risk.

## Using cox_model() and multiple_models()

The Cox model can also be fit using the
[`cox_model()`](https://ckbiobank.github.io/ckbtools/reference/cox_model.md)
function in this package with the expanded data frame

``` r

cox_model(
  participant_data_long,
  exposure   = "bmi_grp",
  outcome    = "status",
  time_in    = "time_in",
  time_out   = "time_out",
  adjust_for = "strata(agegrp)",
  return     = "coefs"
)
```

``` outp
# # A tibble: 4 × 10
#   term  estimate std.error statistic   p.value     n nevent n_group nevent_group
#   <chr>    <dbl>     <dbl>     <dbl>     <dbl> <int>  <int>   <dbl>        <dbl>
# 1 bmi_…   NA       NA          NA    NA        71887   2098    6199          121
# 2 bmi_…    0.458    0.0952      4.81  1.53e- 6    NA     NA   45902         1269
# 3 bmi_…    0.728    0.0992      7.33  2.24e-13    NA     NA   18584          642
# 4 bmi_…    1.20     0.153       7.84  4.68e-15    NA     NA    1202           66
# # ℹ 1 more variable: formula <chr>
```

If using the
[`multiple_models()`](https://ckbiobank.github.io/ckbtools/reference/multiple_models.md)
function, age-at-risk expansion and stratification can be used via the
`age_at_risk` argument. A list of arguments for the
[`expand_age_at_risk()`](https://ckbiobank.github.io/ckbtools/reference/expand_age_at_risk.md)
function is supplied as a named list.

Note the it is the original non-expanded data set that must be used, and
that `"strata(agegrp)"` is automatically added to the model formula.

``` r

results <- multiple_models(
  ckbtools_participant_data,
  model       = "cox",
  exposures   = "bmi_grp",
  outcomes    = "ihd",
  age_at_risk = list(ages = seq(from = 40, to = 80, by = 5)),
  return      = "coefs",
  verbose     = "quiet"
)

results[, c("term", "estimate", "std.error", "n_group", "nevent_group")]
```

``` outp
# # A tibble: 4 × 5
#   term           estimate std.error n_group nevent_group
#   <chr>             <dbl>     <dbl>   <dbl>        <dbl>
# 1 bmi_grp<18.5     NA       NA         6199          121
# 2 bmi_grp18.5-25    0.458    0.0952   45902         1269
# 3 bmi_grp25-30      0.728    0.0992   18584          642
# 4 bmi_grp30+        1.20     0.153     1202           66
```

If there is only one outcome, then
[`multiple_models()`](https://ckbiobank.github.io/ckbtools/reference/multiple_models.md)
can be used with the expanded data frame, but the arguments `outcomes`,
`time_in` and `time_out` must be set appropriately.

``` r

results <- multiple_models(
  participant_data_long,
  model      = "cox",
  exposures  = "bmi_grp",
  outcomes   = "status",
  time_in    = "time_in",
  time_out   = "time_out",
  adjust_for = list(none = "strata(agegrp)"),
  return     = "coefs",
  verbose     = "quiet"
)

results[, c("term", "estimate", "std.error", "n_group", "nevent_group")]
```

``` outp
# # A tibble: 4 × 5
#   term           estimate std.error n_group nevent_group
#   <chr>             <dbl>     <dbl>   <dbl>        <dbl>
# 1 bmi_grp<18.5     NA       NA         6199          121
# 2 bmi_grp18.5-25    0.458    0.0952   45902         1269
# 3 bmi_grp25-30      0.728    0.0992   18584          642
# 4 bmi_grp30+        1.20     0.153     1202           66
```

## References

[San H, Lewington S. Lexis expansion-age-at-risk adjustment for survival
analysis. InBrussels: Statistics and Pharmacokinetics Pharmaceutical
Users Software Exchange conference 2013 Oct
13.](https://www.lexjansen.com/phuse/2013/sp/SP09.pdf)

[Moolgavkar, S. H., Chang, E. T., Watson, H. N. & Lau, E. C. An
Assessment of the Cox Proportional Hazards Regression Model for
Epidemiologic Studies. Risk Analysis 38, 777–794
(2018).](https://doi.org/10.1111/risa.12865)

[Vyas, M. V., Fang, J., Kapral, M. K. & Austin, P. C. Choice of
time-scale in time-to-event analysis: evaluating age-dependent
associations. Annals of Epidemiology 62, 69–76
(2021).](https://doi.org/10.1016/j.annepidem.2021.06.006) (Not open
access)

[Hurley, M. A. A reference relative time-scale as an alternative to
chronological age for cohorts with long follow-up. Emerg Themes
Epidemiol 12, 18 (2015).](https://doi.org/10.1186/s12982-015-0043-6)

[Canchola AJ, Stewart SL, Bernstein L, West DW, Ross RK, Deapen D,
Pinder R, Reynolds P, Wright W, Anton-Culver H, Peel D. Cox regression
using different time-scales. Western Users of SAS Software. San
Francisco, California. 2003 Nov
5.](https://www.lexjansen.com/wuss/2003/DataAnalysis/i-cox_time_scales.pdf)

[Lacey, B. et al. Age-specific association between blood pressure and
vascular and non-vascular chronic diseases in 0·5 million adults in
China: a prospective cohort study. The Lancet Global Health 6, e641–e649
(2018).](https://doi.org/10.1016/S2214-109X(18)30217-1)
