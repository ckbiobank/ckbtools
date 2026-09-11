# Trend and heterogeneity tests

## Trend test

The [`trend()`](https://ckbiobank.github.io/ckbtools/reference/trend.md)
function carries out a chi-squared test for trend (on 1 degree of
freedom) on a vector of independent estimates. The function arguments
are a vector `beta` of independent estimates and a vector `se` of their
standard errors. The function returns the test statistic and p-value in
a named list.

### Example

These data are taken from Figure 3 of the paper ‘Fresh fruit consumption
and major cardiovascular disease in the China Kadoorie Biobank’ by Du,
H., Li, L., Bennett, D. et al. N Engl J Med 2016; 374: 1332–1343.

``` r

b <- c(-0.59784, -0.28768, -0.4943)
s <- c(0.127929, 0.121515, 0.094959)

trend(b, s)
```

``` outp
# $`test statistic`
# [1] 0.1386369
# 
# $p
# [1] 0.7096399
```

## Heterogeneity test

The
[`heterogeneity()`](https://ckbiobank.github.io/ckbtools/reference/heterogeneity.md)
function carries out a a chi-squared test (on $`k - 1`$ degrees of
freedom) for heterogeneity on a vector of $`k`$ independent estimates.
The function arguments are a vector `beta` of independent estimates and
a vector `se` of their standard errors. The function returns the test
statistic, degrees of freedom and p-value in a named list.

### Example

These data are taken from Figure 3 of Du et al.

``` r

beta <- c(-0.47804, -0.46204, -0.40048)
se <- c(0.090005, 0.176823, 0.106565)

heterogeneity(beta, se)
```

``` outp
# $`test statistic`
# [1] 0.3165459
# 
# $`degrees of freedom`
# [1] 2
# 
# $p
# [1] 0.8536168
```

## Methods

### Trend

Suppose that we want to check whether estimates change progressively
from one stratum (subgroup) to the next then we could use a chi-squared
test for trend (on 1 degree of freedom).

The test statistic for trend given a vector of $`k`$ independent
estimates $`\hat{\beta}_i`$ and their standard errors $`\sigma_i`$ is

``` math
 \frac{\sum_{i = 1}^{k} \left(w_i \hat{\beta}_i (i - A) \right)^2}{\sum_{i = 1}^{k}{w_i (i - A)^2}}
```
where

``` math
A = \frac{\sum_{j = 1}^{k}{j w_j}}{\sum_{j = 1}^{k}{w_j}}
```

and $`w_i = 1/\sigma^2`$. This test statistic has a chi-squared
distribution with 1 degree of freedom under the null hypothesis of no
linear trend.

### Heterogeneity

Suppose that information on the estimates for different strata
(subgroups) are to be assessed in order to see if they differ between
strata (i.e. effect modification) then a chi-squared test (on $`k - 1`$
degrees of freedom) for heterogeneity between the estimates for the
different strata can used.

The test statistic for heterogeneity given a vector of $`k`$ independent
estimates $`\hat{\beta}_i`$ and their standard errors $`\sigma_i`$ is

``` math
 \left(\sum_{i=1}^{k} w_i \hat{\beta}_i^2 \right) - \frac{\left(\sum_{i=1}^{k} w_i \hat{\beta}_i \right)^2}{\sum_{i=1}^{k} w_i}, 
```

where $`w_i = 1/\sigma^2`$, which has a chi-squared distribution with
$`k-1`$ degrees of freedom under the null hypothesis.

The test statistic can also be written as:

``` math
\sum_{i=1}^{k} \left[ w_i \left( \hat{\beta}_i - \frac{\sum_{j=1}^{k} w_j\hat{\beta}_j}{\sum_{j=1}^{k} w_j} \right)^2 \right]
```

If there are only two strata (subgroups) then the tests for trend and
heterogeneity are identical.

### Notes

Both these formulae can be expressed in terms of $`(o-e)`$ and $`v`$
using
``` math
\hat{\beta}_i = (o_i - e_i) / v_i
```
and
``` math
\mathrm{se}(\hat{\beta}_i)^2 = 1 / w_i = 1 / v_i
```

First calculate the logrank statistic $`(o-e)`$ and its variance $`v`$
in each separate stratum (subgroup), and their sums:
``` math
O - E = \sum_{i=1}^{k} (o_i - e_i)
```
and
``` math
V = \sum_{i=1}^{k} v_i
```

#### Trend

Define

``` math
m_i = \sum_{i = 1}^{k} i v_i/V
```

and

``` math
T = \sum_{i = 1}^{k} (i-m_i)(o_i - e_i)
```

The variance of $`T`$, $`\mathrm{var}(T)`$, is then

``` math
\mathrm{var}(T) = \sum_{i=1}^{k} (i - m)^2v
```

and the chi-squared test statistic (1 degree of freedom) for trend is
``` math
\frac{T^2}{\mathrm{var}(T)}
```

#### Heterogeneity

The chi-squared test statistic ($`k - 1`$ degrees of freedom) for
heterogeneity is

``` math
\sum_{i=1}^k \left( \frac{(o_i-e_i)^2}{v_i} \right) - \frac{(O-E)^2}{V}
```

## References

[Early Breast Cancer Trialists’ Collaborative Group, Section 5 in
“Treatment of Early Breast Cancer. Volume 1. Worldwide Evidence
1985-1990”](https://www.ctsu.ox.ac.uk/research/the-early-breast-cancer-trialists-collaborative-group-ebctcg/original-methods-for-ebctcg-meta-analyses/section-5-statistical-methods)

[Effects of chemotherapy and hormonal therapy for early breast cancer on
recurrence and 15-year survival: an overview of the randomised trials.
The Lancet 365, 1687–1717
(2005).](https://doi.org/10.1016/S0140-6736(05)66544-0) (Not open
access)

[Armitage, P. Tests for Linear Trends in Proportions and Frequencies.
Biometrics 11, 375–386 (1955).](https://doi.org/10.2307/3001775) (Not
open access)

[Du, H. et al. Fresh Fruit Consumption and Major Cardiovascular Disease
in China. New England Journal of Medicine 374, 1332–1343
(2016).](https://doi.org/10.1056/NEJMoa1501451)
