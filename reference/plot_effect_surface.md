# Plot the one-dimensional direct-effect surface

Renders the fitted surface for a single outcome as it really is —
one-dimensional: each drug's position from harmful to beneficial with
its uncertainty interval, ordered by estimate. No 2-D layout pretends to
represent causal distance. For an anchored fit the placebo = 0 sea-level
line is drawn and effects read as absolute; for a surface fit the plot
is labeled as relative to the arbitrary reference.

## Usage

``` r
plot_effect_surface(fit, scale = c("log", "natural"))
```

## Arguments

- fit:

  A `directeffect_fit` from
  [`fit_surface()`](https://ablack3.github.io/directeffect/reference/fit_surface.md)
  or
  [`anchor_surface()`](https://ablack3.github.io/directeffect/reference/anchor_surface.md);
  both engines plot identically.

- scale:

  `"log"` draws effects on the log scale (0 = reference / sea level);
  `"natural"` back-transforms to the natural scale of a multiplicative
  measure on a logarithmic axis (1 = reference / sea level; log 0
  corresponds to e.g. HR = 1).

## Value

A ggplot object.

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
    requireNamespace("ggplot2", quietly = TRUE)) {
  plot_effect_surface(fit_surface(de, engine = "netmeta"))
}
```
