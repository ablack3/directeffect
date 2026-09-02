# Position a fitted surface against placebo (sea level)

Applies absolute placebo anchors to a fitted relative surface as a
distinct second step, turning relative positions into absolute direct
effects against `placebo = 0`. Anchors retain their uncertainty: they
pull the surface to the right height without pinning any drug exactly,
and their standard errors propagate into every absolute interval.

## Usage

``` r
anchor_surface(fit, anchors = NULL, ...)
```

## Arguments

- fit:

  A `directeffect_fit` from
  [`fit_surface()`](https://ablack3.github.io/directeffect/reference/fit_surface.md).

- anchors:

  Absolute estimates to anchor with (same schema as in
  [`direct_effect_network()`](https://ablack3.github.io/directeffect/reference/direct_effect_network.md)).
  Defaults to the anchors already attached to the fit's network.

- ...:

  Passed to
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
  when re-fitting with the anchored Stan model (`chains`, `iter`,
  `seed`, ...). Unused by the frequentist path. As in
  [`fit_surface()`](https://ablack3.github.io/directeffect/reference/fit_surface.md),
  when no `seed` is supplied the sampler is seeded from the wall clock
  and process id and the caller's global random-number state is left
  untouched.

## Value

A `directeffect_fit` whose effects are absolute: `reference` is
`"placebo"` and each estimate is the drug's direct effect versus placebo
on the log scale.

## Details

The Bayesian path refits the network with the anchored Stan model, in
which no arbitrary identification constraint exists — the anchors
determine the absolute location. The frequentist path estimates the
single location offset that best reconciles the fitted surface with the
anchors by generalized least squares over the anchors' proposals, whose
covariance combines the anchor variances with the surface covariance of
the anchored drugs. Every absolute variance is
`Var(theta_d) + Var(offset) + 2 Cov(theta_d, offset)`, computed from the
fit's full surface covariance, so the absolute standard errors equal the
actual sampling variability of the estimator. With a single anchor the
anchored drug's surface contribution cancels exactly: its absolute
standard error is its anchor's standard error, and the whole fit equals
the joint generalized-least-squares solution of comparisons plus anchor
rows — the same answer the anchored Stan model gives. With several
disagreeing anchors the offset and its variance still equal that joint
solution at the surface's reference drug, but the joint (Bayesian) fit
also lets the anchors update the drugs' relative positions, which a
location shift by design does not; the frequentist path then reports the
honest — slightly larger — sampling variance of its own estimator.

A fit whose network component has no anchor cannot be positioned:
`anchor_surface()` refuses rather than silently picking a sea level.

## See also

[directeffect_formats](https://ablack3.github.io/directeffect/reference/directeffect_formats.md)
for the explicit anchors input schema and the effects table schema.

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
if (requireNamespace("netmeta", quietly = TRUE)) {
  surface <- fit_surface(de, engine = "netmeta")
  absolute <- anchor_surface(surface)
  absolute$effects
}
#>   drug estimate  std_error     lower     upper scale reference  engine
#> 1    A      0.7 0.05715476 0.5879787 0.8120213   log   placebo netmeta
#> 2    B      0.7 0.05715476 0.5879787 0.8120213   log   placebo netmeta
#> 3    C      0.3 0.04000000 0.2216014 0.3783986   log   placebo netmeta
```
