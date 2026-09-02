# Compare the two engines' reconstructions of the same surface

The package's central v0.1 claim is that netmeta and Stan reconstruct
the same direct-effect surface independently. `compare_engines()` makes
that a routine check: it lines the two fits up per drug and reports the
difference and standardized difference between the engines' estimates.
The fits may be passed in either order; they must cover the same drugs
and use the same reference so the comparison is apples-to-apples (refit
with the same `reference`, or anchor both, otherwise).

## Usage

``` r
compare_engines(fit_a, fit_b)

plot_engine_comparison(fit_a, fit_b)
```

## Arguments

- fit_a, fit_b:

  Two `directeffect_fit` objects for the same network — one from the
  `"netmeta"` engine and one from `"stan"`, in either order.

## Value

A data frame with one row per drug: `drug`, `netmeta` (the frequentist
estimate), `stan_mean` (the posterior mean), `difference`
(`stan_mean - netmeta`), and `standardized_difference` (difference over
the combined standard error; `NA` for the reference drug, whose estimate
is exact 0 in both engines). `standardized_difference` is a yardstick
for judging whether a difference is large relative to the estimates'
uncertainty, not a test statistic: the two fits share the same data, so
it has no reference distribution.

`plot_engine_comparison()` returns a ggplot of the netmeta estimates
against the Stan posterior means with the identity line y = x; points on
the line mean the engines agree.

## Examples

``` r
comparisons <- data.frame(
  study_id   = c("S1", "S2", "S3"),
  target     = c("A", "A", "B"),
  comparator = c("B", "C", "C"),
  estimate   = c(0.0, 0.4, 0.4),
  std_error  = c(0.05, 0.05, 0.05)
)
de <- direct_effect_network(comparisons, effect_measure = "HR")
if (requireNamespace("netmeta", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE) && interactive()) {
  comparison <- compare_engines(
    fit_surface(de, engine = "netmeta"),
    fit_surface(de, engine = "stan")
  )
  comparison
}
```
