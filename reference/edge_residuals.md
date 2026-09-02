# Edge residuals: does each comparison agree with the fitted surface?

For every observed comparison, reports the observed effect, the effect
predicted from the fitted surface (`theta_target - theta_comparator`),
their difference, that difference divided by the residual's actual
standard deviation, and the comparison's leverage. Large standardized
residuals mark comparisons that conflict with the surface implied by the
rest of the network. Consumes only the fit contract, so every engine is
treated identically. Rows align with `fit$comparisons`.

## Usage

``` r
edge_residuals(fit)
```

## Arguments

- fit:

  A `directeffect_fit` from
  [`fit_surface()`](https://ablack3.github.io/directeffect/reference/fit_surface.md).

## Value

A data frame with columns `target`, `comparator`, `observed`,
`predicted`, `residual`, `standardized_residual` (`NA` on bridge
comparisons), and `leverage`, one row per comparison in
`fit$comparisons` order.

## Details

The standardization divides by the residual's actual standard deviation
`sqrt(se^2 - Var(predicted))`, not the raw comparison standard error:
the prediction is partly built from the observation itself, so the
residual varies less than the observation does. Standardized this way,
the residuals have unit variance under coherence. `leverage`
(`Var(predicted) / se^2`) says how much of the comparison's variance is
absorbed by its own prediction.

A bridge comparison — the only evidence connecting its endpoints, with
leverage 1 — determines its own prediction completely: nothing in the
network can corroborate or contradict it. Its standardized residual is
reported as `NA`, never as a reassuring 0: "uncheckable" is not
"perfectly consistent". Interpret large standardized residuals as
conflict only on comparisons with leverage clearly below 1; an `NA` row
is uncorroborated evidence whose correctness the network cannot assess.

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
  edge_residuals(fit)
}
#>   target comparator observed     predicted      residual standardized_residual
#> 1      A          B      0.0 -8.881784e-15  8.881784e-15          3.076740e-13
#> 2      A          C      0.4  4.000000e-01 -2.775558e-16         -9.614813e-15
#> 3      B          C      0.4  4.000000e-01 -9.159340e-15         -3.172888e-13
#>    leverage
#> 1 0.6666667
#> 2 0.6666667
#> 3 0.6666667
```
