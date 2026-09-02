# Cycle consistency: do the comparisons agree around closed loops?

Sums the observed effects around each cycle in a cycle basis of the
evidence network. Under perfect consistency every cycle sum is 0; a
cycle whose standardized sum is large contains comparisons that cannot
all be right. Repeat comparisons of the same drug pair are pooled by
precision first, and the basis is built from a spanning tree — every
cycle of the network is a combination of basis cycles, so no cycle
enumeration is ever needed.

## Usage

``` r
cycle_consistency(fit)

plot_cycle_consistency(fit)
```

## Arguments

- fit:

  A `directeffect_fit` from
  [`fit_surface()`](https://ablack3.github.io/directeffect/reference/fit_surface.md).

## Value

A data frame with one row per basis cycle: `cycle` (the drug walk, e.g.
`"A - B - C - A"`), `n_edges`, `inconsistency` (the signed cycle sum on
the log scale), `std_error`, and `z`. Zero rows when the network has no
cycles.

`plot_cycle_consistency()` returns a ggplot showing each basis cycle's
standardized inconsistency with reference lines at z = -1.96, 0, and
1.96.

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
if (requireNamespace("netmeta", quietly = TRUE)) {
  fit <- fit_surface(de, engine = "netmeta")
  cycle_consistency(fit)
}
#>           cycle n_edges inconsistency  std_error z
#> 1 B - A - C - B       3             0 0.08660254 0
```
