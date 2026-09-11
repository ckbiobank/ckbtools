# Standardised incidence rates

## Introduction

The
[`incidence_rates()`](https://ckbiobank.github.io/ckbtools/reference/incidence_rates.md)
function calculates standardised incidence rates, including
standardising by age-at-risk.

The function uses
[`Epi::Lexis()`](https://rdrr.io/pkg/Epi/man/Lexis.html) and
[`Epi::splitLexis()`](https://rdrr.io/pkg/Epi/man/splitLexis.html) to
expand each participant row into separate rows for each age-at-risk
group (defined by the `ages` argument). A new column called `agegrp` is
created which can be used to calculate separate rates for each
age-at-risk group or to standardise the rates by age group.

By default, the function uses the whole data set to create the weights.
Or the `weights` argument can be used to provide weights from a
reference/standard population.

The function returns a list of two data frames:

**incidence_rates_in_all_strata** containing a row for each stratum
formed by the grouping variable(s) (specified by the `group` argument)
and the standardisation variable(s) (specified by the `adjust_for`
argument), with results columns:

- **std_pyar**: Person-years at risk in the reference/standard
  population
- **pyar**: Person-years at risk in this stratum
- **nevents**: Number of events in this stratum
- **incidencerate**: Crude incidence rate in this stratum
- **weight**: Weight for this stratum

**standardised_incidence_rates** containing a row for each group formed
by the grouping variable(s), with results columns:

- **crude_nevents**: Number of events
- **crude_pyar**: Person-years at risk
- **crude_ir**: Crude incidence rate
- **std_ir**: Standardised incidence rate
- **std_ir_se**: Standard error of standardised incidence rate
- **std_ir_lci**: Lower limit of 95% confidence interval for
  standardised incidence rate
- **std_ir_uci**: Upper limit of 95% confidence interval for
  standardised incidence rate

The incidence incidence rates in this data frame are multiplied by a
scale factor set by the `per_pyar` argument (default: 1000).

## Number of events

The function gives a warning when any stratum contains fewer than 20
events, and when any stratum contains zero events. Results may be
unreliable in these cases. See [Confidence
intervals](#confidence-intervals) for available methods for calculating
confidence intervals which may be more suitable when there are few
events.

## Examples

#### Prepare example data

``` r

# Create BMI groups column
ckbtools_participant_data$bmi_grp <- cut(ckbtools_participant_data$bmi_calc,
                                         c(0, 18.5, 25, 30, 100),
                                         labels = c("<18.5", "18.5-25", "25-30", "30+"),
                                         right = FALSE)
```

#### Identify relevant columns in the data set

To calculate incidence rates identify the columns that contain
participants’ date of birth (argument: `col.dob`, default: “dob_anon”),
date of study entry (`col.entry_date`, “study_date”), date of study exit
(`col.exit_date`, “endpoint_date”) and status at exit (`col.status`,
“endpoint”).

#### Overall incidence rate

Not specifying the `group` and `adjust_for` arguments will calculate
overall incidence rate. The standardised incidence rate will equal the
crude incidence rate.

``` r

incidence_rates(ckbtools_participant_data,
                col.exit_date  = "ihd_date",
                col.status     = "ihd")
```

``` outp
# $incidence_rates_in_all_strata
# # A tibble: 1 × 6
#     all std_pyar    pyar nevents incidencerate weight
#   <dbl>    <dbl>   <dbl>   <dbl>         <dbl>  <dbl>
# 1     1  264327. 264327.    2271       0.00859      1
# 
# $standardised_incidence_rates
# # A tibble: 1 × 7
#   crude_nevents crude_pyar crude_ir std_ir std_ir_se std_ir_lci std_ir_uci
#           <dbl>      <dbl>    <dbl>  <dbl>     <dbl>      <dbl>      <dbl>
# 1          2271    264327.     8.59   8.59     0.180       8.24       8.94
```

#### Crude incidence rates

Not specifying the `adjust_for` argument will calculate the crude
incidence rate for each group. The standardised incidence rate for each
group will equal the crude incidence rate for that group. For example,
by BMI group:

``` r

incidence_rates(ckbtools_participant_data,
                col.exit_date = "ihd_date",
                col.status    = "ihd",
                group         = "bmi_grp")
```

``` outp
# $incidence_rates_in_all_strata
# # A tibble: 4 × 7
#     all std_pyar bmi_grp    pyar nevents incidencerate weight
#   <dbl>    <dbl> <fct>     <dbl>   <dbl>         <dbl>  <dbl>
# 1     1  264327. <18.5    22826.     136       0.00596      1
# 2     1  264327. 18.5-25 168758.    1378       0.00817      1
# 3     1  264327. 25-30    68371.     683       0.00999      1
# 4     1  264327. 30+       4372.      74       0.0169       1
# 
# $standardised_incidence_rates
# # A tibble: 4 × 8
#   bmi_grp crude_nevents crude_pyar crude_ir std_ir std_ir_se std_ir_lci
#   <fct>           <dbl>      <dbl>    <dbl>  <dbl>     <dbl>      <dbl>
# 1 <18.5             136     22826.     5.96   5.96     0.511       4.96
# 2 18.5-25          1378    168758.     8.17   8.17     0.220       7.73
# 3 25-30             683     68371.     9.99   9.99     0.382       9.24
# 4 30+                74      4372.    16.9   16.9      1.97       13.1 
# # ℹ 1 more variable: std_ir_uci <dbl>
```

To calculate incidence rates by age-at-risk groups (30-49, 50-59, 60-69,
70-79 years), set `group = "agegrp"`. Any follow-up time outside these
age groups is ignored.

``` r

incidence_rates(ckbtools_participant_data,
                col.exit_date = "ihd_date",
                col.status    = "ihd",
                ages          = c(30, 50, 60, 70, 80),
                group         = "agegrp")
```

``` outp
# $incidence_rates_in_all_strata
# # A tibble: 4 × 7
#     all std_pyar agegrp    pyar nevents incidencerate weight
#   <dbl>    <dbl> <fct>    <dbl>   <dbl>         <dbl>  <dbl>
# 1     1  260331. [30,50) 78644.     184       0.00234      1
# 2     1  260331. [50,60) 85910.     530       0.00617      1
# 3     1  260331. [60,70) 64178.     688       0.0107       1
# 4     1  260331. [70,80) 31599.     708       0.0224       1
# 
# $standardised_incidence_rates
# # A tibble: 4 × 8
#   agegrp  crude_nevents crude_pyar crude_ir std_ir std_ir_se std_ir_lci
#   <fct>           <dbl>      <dbl>    <dbl>  <dbl>     <dbl>      <dbl>
# 1 [30,50)           184     78644.     2.34   2.34     0.172       2.00
# 2 [50,60)           530     85910.     6.17   6.17     0.268       5.64
# 3 [60,70)           688     64178.    10.7   10.7      0.409       9.92
# 4 [70,80)           708     31599.    22.4   22.4      0.842      20.8 
# # ℹ 1 more variable: std_ir_uci <dbl>
```

#### Standardised by sex and age

To calculate incidence rates by BMI group, standardised by sex (column
“is_female”) and 10 year age groups (and excluding follow-up time where
participants are younger than 40 or older than 79):

##### Using the Normal approximation for confidence intervals

``` r

incidence_rates(ckbtools_participant_data,
                col.exit_date  = "ihd_date",
                col.status     = "ihd",
                ages           = seq(40, 80, 10),
                group          = "bmi_grp",
                adjust_for     = c("is_female", "agegrp"),
                method         = "norm")
```

``` outp
# Warning: is_female is numeric and has 2 levels, so treating as a group variable.
# ℹ Change is_female to a factor to suppress this warning.
```

``` outp
# Warning: There are strata with zero events.
```

``` outp
# $incidence_rates_in_all_strata
# # A tibble: 32 × 8
#    is_female agegrp  std_pyar bmi_grp   pyar nevents incidencerate weight
#        <int> <fct>      <dbl> <fct>    <dbl>   <dbl>         <dbl>  <dbl>
#  1         0 [40,50)   28706. <18.5    2304.       0       0       0.115 
#  2         0 [50,60)   36260. <18.5    3272.      14       0.00428 0.145 
#  3         0 [60,70)   26611. <18.5    2759.       9       0.00326 0.107 
#  4         0 [70,80)   12913. <18.5    1700.      27       0.0159  0.0517
#  5         1 [40,50)   39156. <18.5    2628.       3       0.00114 0.157 
#  6         1 [50,60)   49647. <18.5    3681.      19       0.00516 0.199 
#  7         1 [60,70)   37568. <18.5    3367.      20       0.00594 0.151 
#  8         1 [70,80)   18686. <18.5    1889.      29       0.0154  0.0749
#  9         0 [40,50)   28706. 18.5-25 18612.      42       0.00226 0.115 
# 10         0 [50,60)   36260. 18.5-25 23328.     134       0.00574 0.145 
# # ℹ 22 more rows
# 
# $standardised_incidence_rates
# # A tibble: 4 × 8
#   bmi_grp crude_nevents crude_pyar crude_ir std_ir std_ir_se std_ir_lci
#   <fct>           <dbl>      <dbl>    <dbl>  <dbl>     <dbl>      <dbl>
# 1 <18.5             121     21600.     5.60   5.04     0.467       4.13
# 2 18.5-25          1268    159496.     7.95   7.93     0.223       7.50
# 3 25-30             642     64328.     9.98  10.3      0.410       9.54
# 4 30+                66      4123.    16.0   16.6      2.05       12.6 
# # ℹ 1 more variable: std_ir_uci <dbl>
```

##### Using the gamma method for confidence intervals

``` r

incidence_rates(ckbtools_participant_data,
                col.exit_date  = "ihd_date",
                col.status     = "ihd",
                ages           = seq(40, 80, 10),
                group          = "bmi_grp",
                adjust_for     = c("is_female", "agegrp"),
                method         = "gamma")
```

``` outp
# Warning: is_female is numeric and has 2 levels, so treating as a group variable.
# ℹ Change is_female to a factor to suppress this warning.
```

``` outp
# Warning: There are strata with zero events.
```

``` outp
# $incidence_rates_in_all_strata
# # A tibble: 32 × 8
#    is_female agegrp  std_pyar bmi_grp   pyar nevents incidencerate weight
#        <int> <fct>      <dbl> <fct>    <dbl>   <dbl>         <dbl>  <dbl>
#  1         0 [40,50)   28706. <18.5    2304.       0       0       0.115 
#  2         0 [50,60)   36260. <18.5    3272.      14       0.00428 0.145 
#  3         0 [60,70)   26611. <18.5    2759.       9       0.00326 0.107 
#  4         0 [70,80)   12913. <18.5    1700.      27       0.0159  0.0517
#  5         1 [40,50)   39156. <18.5    2628.       3       0.00114 0.157 
#  6         1 [50,60)   49647. <18.5    3681.      19       0.00516 0.199 
#  7         1 [60,70)   37568. <18.5    3367.      20       0.00594 0.151 
#  8         1 [70,80)   18686. <18.5    1889.      29       0.0154  0.0749
#  9         0 [40,50)   28706. 18.5-25 18612.      42       0.00226 0.115 
# 10         0 [50,60)   36260. 18.5-25 23328.     134       0.00574 0.145 
# # ℹ 22 more rows
# 
# $standardised_incidence_rates
# # A tibble: 4 × 7
#   bmi_grp crude_nevents crude_pyar crude_ir std_ir std_ir_lci std_ir_uci
#   <fct>           <dbl>      <dbl>    <dbl>  <dbl>      <dbl>      <dbl>
# 1 <18.5             121     21600.     5.60   5.04       4.17       6.06
# 2 18.5-25          1268    159496.     7.95   7.93       7.50       8.38
# 3 25-30             642     64328.     9.98  10.3        9.56      11.2 
# 4 30+                66      4123.    16.0   16.6       12.8       21.2
```

Note how the standardised incidence rates are exactly the same, but the
confidence interval limits are slightly different.

#### External weights

To calculate incidence rates by BMI group, standardised by 10 year age
groups with external weights:

``` r

std_pop_w <- data.frame(agegrp = c("[40,50)",
                                   "[50,60)",
                                   "[60,70)",
                                   "[70,80)"),
                        std_pyar = c(50, 30, 20, 10))

incidence_rates(ckbtools_participant_data,
                col.exit_date  = "ihd_date",
                col.status     = "ihd",
                ages           = seq(40, 80, 10),
                group          = "bmi_grp",
                adjust_for     = "agegrp",
                weights        = std_pop_w,
                method         = "gamma")
```

``` outp
# Warning: There are strata with fewer than 20 events.
```

``` outp
# $incidence_rates_in_all_strata
# # A tibble: 16 × 7
#    agegrp  std_pyar bmi_grp   pyar nevents incidencerate weight
#    <chr>      <dbl> <fct>    <dbl>   <dbl>         <dbl>  <dbl>
#  1 [40,50)       50 <18.5    4933.       3      0.000608 0.455 
#  2 [50,60)       30 <18.5    6953.      33      0.00475  0.273 
#  3 [60,70)       20 <18.5    6127.      29      0.00473  0.182 
#  4 [70,80)       10 <18.5    3588.      56      0.0156   0.0909
#  5 [40,50)       50 18.5-25 43114.      99      0.00230  0.455 
#  6 [50,60)       30 18.5-25 54915.     314      0.00572  0.273 
#  7 [60,70)       20 18.5-25 41260.     441      0.0107   0.182 
#  8 [70,80)       10 18.5-25 20206.     414      0.0205   0.0909
#  9 [40,50)       50 25-30   18618.      63      0.00338  0.455 
# 10 [50,60)       30 25-30   22617.     169      0.00747  0.273 
# 11 [60,70)       20 25-30   15774.     194      0.0123   0.182 
# 12 [70,80)       10 25-30    7318.     216      0.0295   0.0909
# 13 [40,50)       50 30+      1197.       6      0.00501  0.455 
# 14 [50,60)       30 30+      1423.      14      0.00984  0.273 
# 15 [60,70)       20 30+      1017.      24      0.0236   0.182 
# 16 [70,80)       10 30+       486.      22      0.0453   0.0909
# 
# $standardised_incidence_rates
# # A tibble: 4 × 7
#   bmi_grp crude_nevents crude_pyar crude_ir std_ir std_ir_lci std_ir_uci
#   <fct>           <dbl>      <dbl>    <dbl>  <dbl>      <dbl>      <dbl>
# 1 <18.5             121     21600.     5.60   3.85       3.16       4.73
# 2 18.5-25          1268    159496.     7.95   6.41       6.04       6.80
# 3 25-30             642     64328.     9.98   8.50       7.83       9.22
# 4 30+                66      4123.    16.0   13.4       10.2       17.4
```

Note that there must be columns in the `weights` data frame that match
the `adjust_for` variables, and one row for each combination of levels
of those variables, and a column named `std_pyar` with the person-years
at risk for that stratum in the reference/standard population.

#### Multiple groups

Multiple columns can be specified in `group` to calculate incidence
rates for all combinations. For example, to calculate incidence rates by
sex and 10 year age groups (between 30 and 80 years):

``` r

incidence_rates(ckbtools_participant_data,
                col.dob        = "dob_anon",
                col.entry_date = "study_date",
                col.exit_date  = "ihd_date",
                col.status     = "ihd",
                ages           = seq(30, 80, 10),
                group          = c("is_female", "agegrp"))
```

``` outp
# Warning: is_female is numeric and has 2 levels, so treating as a group variable.
# ℹ Change is_female to a factor to suppress this warning.
```

``` outp
# Warning: There are strata with fewer than 20 events.
```

``` outp
# $incidence_rates_in_all_strata
# # A tibble: 10 × 8
#      all std_pyar is_female agegrp    pyar nevents incidencerate weight
#    <dbl>    <dbl>     <int> <fct>    <dbl>   <dbl>         <dbl>  <dbl>
#  1     1  260331.         0 [30,40)  4596.      10      0.00218       1
#  2     1  260331.         0 [40,50) 28706.      62      0.00216       1
#  3     1  260331.         0 [50,60) 36260.     213      0.00587       1
#  4     1  260331.         0 [60,70) 26611.     278      0.0104        1
#  5     1  260331.         0 [70,80) 12913.     290      0.0225        1
#  6     1  260331.         1 [30,40)  6188.       3      0.000485      1
#  7     1  260331.         1 [40,50) 39156.     109      0.00278       1
#  8     1  260331.         1 [50,60) 49647.     317      0.00639       1
#  9     1  260331.         1 [60,70) 37568.     410      0.0109        1
# 10     1  260331.         1 [70,80) 18686.     418      0.0224        1
# 
# $standardised_incidence_rates
# # A tibble: 10 × 9
#    is_female agegrp  crude_nevents crude_pyar crude_ir std_ir std_ir_se
#        <int> <fct>           <dbl>      <dbl>    <dbl>  <dbl>     <dbl>
#  1         0 [30,40)            10      4596.    2.18   2.18      0.688
#  2         0 [40,50)            62     28706.    2.16   2.16      0.274
#  3         0 [50,60)           213     36260.    5.87   5.87      0.402
#  4         0 [60,70)           278     26611.   10.4   10.4       0.627
#  5         0 [70,80)           290     12913.   22.5   22.5       1.32 
#  6         1 [30,40)             3      6188.    0.485  0.485     0.280
#  7         1 [40,50)           109     39156.    2.78   2.78      0.267
#  8         1 [50,60)           317     49647.    6.39   6.39      0.359
#  9         1 [60,70)           410     37568.   10.9   10.9       0.539
# 10         1 [70,80)           418     18686.   22.4   22.4       1.09 
# # ℹ 2 more variables: std_ir_lci <dbl>, std_ir_uci <dbl>
```

## Methods

If there have been $`n`$ events during $`t`$ person-years at risk, then
the incidence rate is given by
``` math
r = \frac{n}{t} \text{ .}
```
The standardised incidence rate across strata (indexed by $`i`$) is
calculated as
``` math
\textrm{SIR} = \frac{\sum_i{r_i \times s_i}}{\sum_i{s_i}}
```
where $`r_i = n_i / t_i`$ and $`s_i`$ are the person-years at risk for
each stratum in a standard population.

By default, this function will calculate the $`s_i`$ using all
participants in the produced data, and calculate standard incidence
rates for groups defined by `group`.

The function returns the weight for each stratum calculated as
$`s_i/ \sum_i{s_i}`$.

### Confidence intervals

By treating the number of events in each stratum as a Poisson variable,
the variance of the SIR can be calculated by
``` math
\textrm{var}(\textrm{SIR}) = \sum_i \left( \left(\frac{w_i}{\sum w_i}\right)^2 \times \frac{r_i}{t_i}\right) \text{ ,}
```
then the standard error is given by $`\sqrt{\textrm{var}}`$. Using this
standard error to construct confidence intervals by a normal
approximation may not be suitable when there are few events or rates are
close to zero, or there are a small number of events within strata.

Other methods for calculating confidence intervals including the gamma
method of Fay and Feuer (1997). This is available in the
`wspoissonTest()` function in the `asht` package, and can be used in
[`incidence_rates()`](https://ckbiobank.github.io/ckbtools/reference/incidence_rates.md)
by setting the argument `method = "gamma"`.

Using a normal approximation method is not suitable when there are fewer
than 100 events. Methods similar to the gamma method of Fay and Feuer
(1997) have been shown to be suitable when there are 10 or more events
(Morris et al. 2018).

### Date calculations using Epi package

Dates are converted to a numerical value, giving the calendar year as a
fractional number, assuming that years are all 365.25 days long, so
inaccuracies may arise see \[Epi::cal.yr\] for more details

## References

[Armitage, Berry, and Matthews. *Statistical Methods in Medical
Research*. 2002. (Chapter: Statistical methods in epidemiology. Section:
19.3 Rates and standardization. Pages 659 -
633.)](https://www.wiley.com/en-ae/Statistical+Methods+in+Medical+Research%2C+4th+Edition-p-9780632052578)
(Not open access)

[Fay, M. P. & Feuer, E. J. Confidence Intervals for Directly
Standardized Rates: A Method Based on the Gamma Distribution. Statistics
in Medicine 16, 791–801
(1997).](https://doi.org/10.1002/(SICI)1097-0258(19970415)16:7%3C791::AID-SIM500%3E3.0.CO;2-%23)
(Not open access)

[Morris, J. K., Tan, J., Fryers, P. & Bestwick, J. Evaluation of
stability of directly standardized rates for sparse data using
simulation methods. Popul Health Metrics 16, 19
(2018).](https://doi.org/10.1186/s12963-018-0177-1)

[Ng, H. K. T., Filardo, G. & Zheng, G. Confidence interval estimating
procedures for standardized incidence rates. Computational Statistics &
Data Analysis 52, 3501–3516
(2008).](https://doi.org/10.1016/j.csda.2007.11.004) (Not open access)
