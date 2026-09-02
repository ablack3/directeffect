# Plot the surface positioned against sea level

The signature figure of the package's conceptual decomposition: the
relative structure inferred from comparative evidence, positioned
against the placebo = 0 sea-level line by the absolute anchors. Drugs
spread along the horizontal axis in order of effect; their heights are
the absolute direct effects. Rendering stays one-dimensional in the
effect — the horizontal ordering is presentation, not geometry. Only an
anchored fit can be drawn: a relative surface has no sea level.

## Usage

``` r
plot_sea_level(fit, scale = c("log", "natural"))
```

## Arguments

- fit:

  An anchored `directeffect_fit` from
  [`anchor_surface()`](https://ablack3.github.io/directeffect/reference/anchor_surface.md).

- scale:

  As in
  [`plot_effect_surface()`](https://ablack3.github.io/directeffect/reference/plot_effect_surface.md).

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
anchors <- data.frame(
  study_id  = "RCT1",
  drug      = "C",
  reference = "placebo",
  estimate  = 0.3,
  std_error = 0.04
)
de <- direct_effect_network(comparisons, anchors = anchors,
                            effect_measure = "HR")
if (requireNamespace("netmeta", quietly = TRUE) &&
    requireNamespace("ggplot2", quietly = TRUE)) {
  plot_sea_level(anchor_surface(fit_surface(de, engine = "netmeta")))
}
```
