# Measure how well a fit recovers simulation truth

Compares a fit against the known true effects of the simulation that
generated its data. Works from the fit contract only, so it treats every
engine identically. For a surface fit the truth is re-centred at the
fit's arbitrary reference drug before comparison, and the reference row
itself (estimated as exactly 0 by construction) is excluded from all
four metrics — bias, RMSE, coverage, and rank correlation — so its
phantom exact (0, 0) pair cannot flatter any of them. For an anchored
fit (`reference = "placebo"`) the truth is already on the placebo = 0
scale, so it is compared directly and every drug contributes.

## Usage

``` r
validate_recovery(fit, simulation)
```

## Arguments

- fit:

  A `directeffect_fit` from
  [`fit_surface()`](https://ablack3.github.io/directeffect/reference/fit_surface.md),
  fitted to `simulation$network`.

- simulation:

  The result of
  [`simulate_direct_effect_network()`](https://ablack3.github.io/directeffect/reference/simulate_direct_effect_network.md).

## Value

A list with `bias` (mean error), `rmse`, `coverage` (proportion of
intervals containing the true value), `rank_correlation` (Spearman
correlation of estimated and true effect ordering), `n_drugs` (drugs
contributing to the summaries), and `reference`.

## Examples

``` r
simulation <- simulate_direct_effect_network(
  n_drugs = 5, n_comparisons = 12, n_anchors = 0,
  heterogeneity = 0, seed = 1
)
if (requireNamespace("netmeta", quietly = TRUE)) {
  fit <- fit_surface(simulation$network, engine = "netmeta")
  validate_recovery(fit, simulation)
}
#> $bias
#> [1] 0.02129941
#> 
#> $rmse
#> [1] 0.04496959
#> 
#> $coverage
#> [1] 1
#> 
#> $rank_correlation
#> [1] 1
#> 
#> $n_drugs
#> [1] 4
#> 
#> $reference
#> [1] "drug_01"
#> 
```
