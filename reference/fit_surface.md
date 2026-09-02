# Fit the relative direct-effect surface

Estimates the relative positions of all drugs in a connected
direct-effect network from the comparative estimates alone. Anchors are
deliberately ignored at this stage: the surface answers "are my
comparative estimates internally coherent, and what relative structure
do they imply?" — not where that structure sits in absolute terms.
Identification uses an arbitrary reference drug fixed at 0; use
[`anchor_surface()`](https://ablack3.github.io/directeffect/reference/anchor_surface.md)
to position the surface absolutely.

## Usage

``` r
fit_surface(de, engine = c("netmeta", "stan"), reference = NULL, ...)
```

## Arguments

- de:

  A `directeffect_network` created by
  [`direct_effect_network()`](https://ablack3.github.io/directeffect/reference/direct_effect_network.md).
  Must be a single connected component; fit components separately (see
  [`check_connectivity()`](https://ablack3.github.io/directeffect/reference/check_connectivity.md))
  otherwise.

- engine:

  `"netmeta"` (frequentist network meta-analysis) or `"stan"`
  (Bayesian). For networks with one comparison per study, both engines
  fit the identical common-effect likelihood
  `y_k ~ N(theta_target - theta_comparator, se_k^2)` and return the same
  fit contract. Multi-arm evidence (several comparisons sharing a
  `study_id`) is supported by the `"netmeta"` engine only in this
  version: netmeta models the correlation between contrasts sharing an
  arm, while the Stan engine refuses such networks rather than mis-treat
  the rows as independent.

- reference:

  Drug fixed at 0 for identification. Defaults to the first treatment
  alphabetically. The choice is arbitrary and does not affect any
  estimated difference between drugs.

- ...:

  Engine-specific options. The Stan engine accepts `chains` (default 4),
  `iter` (default 2000), `seed`, `refresh` (default 0), and any further
  argument to
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).
  When no `seed` is supplied the sampler is seeded from the wall clock
  and process id; in either case the caller's global random-number state
  is left untouched, so fitting never disturbs the session's
  reproducibility. The netmeta engine accepts none.

## Value

An object of class `directeffect_fit` with components `effects` (tidy
per-drug table: `drug`, `estimate`, `std_error`, `lower`, `upper`,
`scale`, `reference`, `engine`), `covariance` (the full covariance of
the estimated effects; see
[directeffect_formats](https://ablack3.github.io/directeffect/reference/directeffect_formats.md)),
`comparisons`, `anchors`, `heterogeneity`, `diagnostics`, `engine`,
`engine_fit` (the raw engine object — the only place engine internals
appear), and `network`.

## See also

[directeffect_formats](https://ablack3.github.io/directeffect/reference/directeffect_formats.md)
for the explicit schema of every column in the effects table and the
other fit components.

## Examples

``` r
comparisons <- data.frame(
  study_id   = c("S1", "S2", "S3"),
  target     = c("A", "A", "B"),
  comparator = c("B", "C", "C"),
  estimate   = c(log(1.02), log(1.34), log(1.29)),
  std_error  = c(0.07, 0.09, 0.08)
)
de <- direct_effect_network(comparisons, effect_measure = "HR")
if (requireNamespace("netmeta", quietly = TRUE)) {
  fit <- fit_surface(de, engine = "netmeta")
  fit$effects
}
#>   drug    estimate  std_error      lower       upper scale reference  engine
#> 1    A  0.00000000 0.00000000  0.0000000  0.00000000   log         A netmeta
#> 2    B -0.02440579 0.06051753 -0.1430180  0.09420638   log         A netmeta
#> 3    C -0.28506030 0.06868800 -0.4196863 -0.15043430   log         A netmeta
```
