# Debias a fitted surface using a network-level bias estimate

Per-edge calibration (e.g. empirical calibration against a
negative-control panel) corrects each comparison in isolation. When a
drug's own comparisons are consistently biased in the same direction
across many edges, inverse- variance pooling in
[`fit_surface()`](https://ablack3.github.io/directeffect/reference/fit_surface.md)
treats that consistency as *precision* rather than *shared bias* — it
concentrates the bias into a spuriously narrow surface position instead
of cancelling it out.

## Usage

``` r
debias_surface(fit, bias_fit, independent = TRUE)
```

## Arguments

- fit:

  A `directeffect_fit` for the outcome of interest, from
  [`fit_surface()`](https://ablack3.github.io/directeffect/reference/fit_surface.md).

- bias_fit:

  A second `directeffect_fit`, built by calling
  [`direct_effect_network()`](https://ablack3.github.io/directeffect/reference/direct_effect_network.md)
  and
  [`fit_surface()`](https://ablack3.github.io/directeffect/reference/fit_surface.md)
  on the same comparisons (same drugs, same edges) but with
  `estimate`/`std_error` replaced by each edge's own bias signal (e.g.
  the mean and standard error of the mean of its negative-control
  panel's log effect). Every drug's fitted position here is an estimate
  of its systematic bias, not a real effect.

- independent:

  Are `fit` and `bias_fit`'s sampling errors independent? Default
  `TRUE`, in which case `var(debiased) = var(fit) + var(bias_fit)`. This
  likely *understates* the true variance somewhat, since both fits
  typically come from the same underlying cohorts and design and so
  share correlated error sources that this simple sum does not capture —
  there is currently no general way to estimate the covariance between
  two independently-fit surfaces, so that understatement is a known
  limitation, not something this function corrects for.
  `independent = FALSE` is not yet implemented and errors.

## Value

A `directeffect_fit` restricted to the drugs `fit` and `bias_fit` have
in common (a warning names any dropped), with `estimate` equal to
`fit`'s estimate minus `bias_fit`'s estimate per drug and `covariance`
equal to the (assumed-independent) sum. `$bias_fit` on the result
carries the bias surface used, for provenance.

## Details

`debias_surface()` corrects for this by taking a second surface fit on
the *same* comparison network, but with the outcome of interest replaced
by something whose true effect is known to be null for every drug —
typically the mean of each edge's own negative-control panel. Because
every drug's true position on that second surface is 0, any nonzero
fitted position is a network-propagated estimate of that drug's
systematic bias: connectivity- aware, unlike a raw per-drug average of
edge-level bias, which ignores how bias can concentrate through well-
versus poorly-connected nodes. Subtracting it from the surface of
interest is a network-level analogue of empirical calibration —
corrected once, for the whole surface, rather than edge by edge.

## Examples

``` r
comparisons <- data.frame(
  study_id   = c("S1", "S2", "S3"),
  target     = c("A", "A", "B"),
  comparator = c("B", "C", "C"),
  estimate   = c(0.5, 0.9, 0.4),
  std_error  = c(0.05, 0.05, 0.05)
)
bias <- data.frame(
  study_id   = c("S1", "S2", "S3"),
  target     = c("A", "A", "B"),
  comparator = c("B", "C", "C"),
  estimate   = c(0.0, 0.2, 0.2),
  std_error  = c(0.05, 0.05, 0.05)
)
if (requireNamespace("netmeta", quietly = TRUE)) {
  fit <- fit_surface(direct_effect_network(comparisons), engine = "netmeta")
  bias_fit <- fit_surface(direct_effect_network(bias), engine = "netmeta")
  debiased <- debias_surface(fit, bias_fit)
  debiased$effects
}
#>   drug estimate  std_error      lower      upper scale reference  engine
#> 1    A      0.0 0.00000000  0.0000000  0.0000000   log         A netmeta
#> 2    B     -0.5 0.05773503 -0.6131586 -0.3868414   log         A netmeta
#> 3    C     -0.7 0.05773503 -0.8131586 -0.5868414   log         A netmeta
```
