# Pool comparative estimates by inverse-variance meta-analysis

Combines a drug's comparative estimates by weighting each by its
precision, the standard fixed-effect or random-effects
(DerSimonian–Laird) meta-analysis.

## Usage

``` r
pool_meta(
  de,
  method = c("fixed", "random"),
  balance = NULL,
  max_smd = 0.1,
  min_balanced_frac = 0.1
)
```

## Arguments

- de:

  A `directeffect_network`, or a data frame with `target`, `comparator`,
  `estimate` and `std_error`.

- method:

  `"fixed"` for a common-effect inverse-variance pool, or `"random"` for
  DerSimonian–Laird, which adds a between-comparison variance `tau^2`
  estimated by the method of moments.

- balance, max_smd, min_balanced_frac:

  Optional balance filter, applied exactly as in
  [`pool_bootstrap()`](https://ablack3.github.io/directeffect/reference/pool_bootstrap.md),
  so the two can be compared on identical subsets.

## Value

A data frame with one row per drug appearing as a target: `drug`,
`n_trials`, `n_balanced`, `frac_balanced`, `pooled_hr`, `ci_lower`,
`ci_upper`, `log_hr`, `se_log_hr`, `p_value` (two-sided), `tau2`,
`i_squared`, `q_p_value` (Cochran's Q test of homogeneity), `method`,
`reported`.

## Details

This is the efficient version of what
[`pool_bootstrap()`](https://ablack3.github.io/directeffect/reference/pool_bootstrap.md)
does crudely. Both average a drug's estimates as they stand;
meta-analysis weights each by `1/se^2` instead of equally, and its
interval comes from the pooled variance rather than from resampling
trials.

**The two share an estimand, and therefore share a bias.** Both answer
"how did this drug fare against the comparators it happened to be tested
against?" Neither can distinguish a protective drug from one that was
simply compared against harmful drugs, because a comparator's own effect
never enters the calculation. Weighting fixes efficiency; it does not
fix that. Only solving the network
([`fit_surface()`](https://ablack3.github.io/directeffect/reference/fit_surface.md))
estimates the comparators' positions rather than assuming them to be
zero.

## See also

[`pool_bootstrap()`](https://ablack3.github.io/directeffect/reference/pool_bootstrap.md),
[`fit_surface()`](https://ablack3.github.io/directeffect/reference/fit_surface.md),
[`compare_pooling()`](https://ablack3.github.io/directeffect/reference/compare_pooling.md).

## Examples

``` r
de <- direct_effect_network(example_comparisons, effect_measure = "HR")
pool_meta(de, method = "random")
#>           drug n_trials n_balanced frac_balanced pooled_hr  ci_lower  ci_upper
#> 1 atorvastatin        4          4             1 0.9112081 0.8226495 1.0093000
#> 2   lovastatin        1          1             1 1.0140985 0.8336013 1.2336781
#> 3  pravastatin        2          2             1 0.9088928 0.7677309 1.0760101
#> 4 rosuvastatin        3          3             1 0.8907608 0.8115541 0.9776980
#> 5  simvastatin        2          2             1 0.8916816 0.8001586 0.9936731
#>        log_hr  se_log_hr    p_value        tau2 i_squared  q_p_value method
#> 1 -0.09298400 0.05216377 0.07466124 0.005648022 0.5345751 0.09182779 random
#> 2  0.01400000 0.10000000 0.88865999 0.000000000        NA         NA random
#> 3 -0.09552809 0.08611628 0.26730427 0.005442000 0.3435172 0.21712543 random
#> 4 -0.11567934 0.04751271 0.01490401 0.000000000 0.0000000 0.60939432 random
#> 5 -0.11464615 0.05525466 0.03799871 0.000000000 0.0000000 0.50505061 random
#>   reported
#> 1     TRUE
#> 2     TRUE
#> 3     TRUE
#> 4     TRUE
#> 5     TRUE
```
