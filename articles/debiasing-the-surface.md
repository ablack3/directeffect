# Debiasing the surface with negative controls

Comparative estimates can be confounded — by indication, by unmeasured
severity, by whatever a negative-control panel is meant to detect. The
standard fix is empirical calibration: for each comparison, fit a null
distribution from a panel of outcomes with no plausible causal
relationship to the exposures, and use it to widen or shift that one
comparison’s interval. That is a *per-edge* correction.
[`debias_surface()`](https://ablack3.github.io/directeffect/reference/debias_surface.md)
does something different, and the two are not interchangeable.

## What per-edge calibration misses

[`fit_surface()`](https://ablack3.github.io/directeffect/reference/fit_surface.md)
pools comparisons by inverse variance. If a drug’s own comparisons are
biased in the same direction across several edges — plausible whenever
the bias comes from something about the drug itself (its indication, its
typical patient mix) rather than from any one trial — inverse-variance
pooling cannot tell that apart from genuine precision. Several biased
edges agreeing with each other look, to the pooling step, exactly like
several unbiased edges agreeing with each other: the surface gets *more*
confident, not less, about the wrong number.

[`debias_surface()`](https://ablack3.github.io/directeffect/reference/debias_surface.md)
corrects this at the network level instead of the edge level. Fit the
surface twice on the same comparisons: once for the outcome of interest,
once with every edge’s `estimate`/`std_error` replaced by its own
negative-control panel’s summary (mean and SE of the panel’s log effect
— a quantity whose true value is 0 for every drug by construction).
Whatever nonzero position a drug gets on that second surface is its
network-propagated systematic bias. Subtract it.

## A worked example

Four drugs, `A` fixed as reference, plus a fifth (`E`) that matters for
a reason covered below. `C` is compared cleanly to `A` and `B`. `D` is
compared to all three of `A`, `B`, and `C` — and every one of those
three comparisons happens to share a **+0.3** bias, say from something
common to how `D`’s trials were run:

``` r

library(directeffect)

comparisons <- data.frame(
  study_id   = c("S1", "S2", "S3", "S4", "S5", "S6", "S7"),
  target     = c("A",  "A",  "B",  "A",  "B",  "C",  "D"),
  comparator = c("B",  "C",  "C",  "D",  "D",  "D",  "E"),
  estimate   = c(0.0, -0.2, -0.2,  0.2,  0.2,  0.4,  0.2),
  std_error  = rep(0.05, 7)
)

de  <- direct_effect_network(comparisons, effect_measure = "HR")
fit <- fit_surface(de, engine = "netmeta")
fit$effects[, c("drug", "estimate", "std_error", "lower", "upper")]
#>   drug      estimate  std_error       lower       upper
#> 1    A  0.000000e+00 0.00000000  0.00000000  0.00000000
#> 2    B  6.661338e-15 0.03535534 -0.06929519  0.06929519
#> 3    C  2.000000e-01 0.03535534  0.13070481  0.26929519
#> 4    D -2.000000e-01 0.03535534 -0.26929519 -0.13070481
#> 5    E -4.000000e-01 0.06123724 -0.52002279 -0.27997721
```

`D`’s three comparisons corroborate each other perfectly (they were
generated to), so the surface reports a *tighter* interval for `D` than
for `C`, which has only two clean comparisons — even though `D`’s point
estimate is the biased one. That is the exact pathology from the section
above, not a contrived edge case: nothing here looks wrong until you
know the truth.

## The bias network

A negative-control panel exists for every comparison except `D`–`E`
(suppose that trial’s panel wasn’t run, or didn’t clear a minimum
control count). Each panel’s mean log effect stands in for `estimate`;
`A`, `B`, `C`’s panels are centered at 0, `D`’s three are centered at
the injected 0.3:

``` r

bias <- data.frame(
  study_id   = c("S1", "S2", "S3", "S4", "S5", "S6"),
  target     = c("A",  "A",  "B",  "A",  "B",  "C"),
  comparator = c("B",  "C",  "C",  "D",  "D",  "D"),
  estimate   = c(0.0,  0.0,  0.0,  0.3,  0.3,  0.3),
  std_error  = rep(0.05, 6)
)

bias_de  <- direct_effect_network(bias, effect_measure = "HR")
bias_fit <- fit_surface(bias_de, engine = "netmeta")
bias_fit$effects[, c("drug", "estimate", "std_error")]
#>   drug      estimate  std_error
#> 1    A  0.000000e+00 0.00000000
#> 2    B  1.332268e-14 0.03535534
#> 3    C  6.661338e-15 0.03535534
#> 4    D -3.000000e-01 0.03535534
```

The bias surface recovers what was actually injected: `D` sits at `-0.3`
(the sign flips because `D` was consistently the *comparator* in its
biased edges above), everyone else at 0.

## Debiasing

``` r

debiased <- debias_surface(fit, bias_fit)
#> Warning: 1 of 5 drug(s) in `fit` have no bias estimate in `bias_fit` and are
#> dropped: E
debiased$effects[, c("drug", "estimate", "std_error", "lower", "upper")]
#>   drug      estimate std_error        lower     upper
#> 1    A  0.000000e+00      0.00  0.000000000 0.0000000
#> 2    B -6.661338e-15      0.05 -0.097998199 0.0979982
#> 3    C  2.000000e-01      0.05  0.102001801 0.2979982
#> 4    D  1.000000e-01      0.05  0.002001801 0.1979982
```

`D`’s estimate moves from the biased `-0.2` back to `0.1` — the value
the comparisons would have shown without the shared confound — and its
interval *widens* rather than staying artificially tight, because the
debiased variance correctly adds the bias surface’s own uncertainty back
in. The false precision from three agreeing-but-biased edges is gone.

## Shrinking the bias estimate

[`debias_surface()`](https://ablack3.github.io/directeffect/reference/debias_surface.md)
above subtracted `D`’s network-propagated bias estimate (`-0.3`) exactly
as fit. That estimate came from only three edges — precise enough here
because they were generated to agree perfectly, but in practice a drug
reached mostly through connectivity rather than direct negative-control
evidence can have a bias estimate that is itself noisy. Subtracting a
noisy estimate at face value can add noise instead of removing bias.

`shrink` pulls each drug’s bias estimate toward a common target before
subtracting, by DerSimonian–Laird empirical-Bayes shrinkage — the same
between-group variance estimator `pool_meta(method = "random")` uses,
but estimated across drugs’ bias estimates here rather than across one
drug’s trials:

``` r

shrunk <- debias_surface(fit, bias_fit, shrink = "zero")
#> Warning: 1 of 5 drug(s) in `fit` have no bias estimate in `bias_fit` and are
#> dropped: E
shrunk$shrinkage
#> $method
#> [1] "zero"
#> 
#> $tau2
#> [1] 0.02875
#> 
#> $target
#> [1] 0
#> 
#> $weight
#>         A         B         C         D 
#> 1.0000000 0.9583333 0.9583333 0.9583333
shrunk$effects[, c("drug", "estimate", "std_error")]
#>   drug      estimate  std_error
#> 1    A  0.000000e+00 0.00000000
#> 2    B -6.106227e-15 0.04947643
#> 3    C  2.000000e-01 0.04947643
#> 4    D  8.750000e-02 0.04947643
```

`A` — `bias_fit`’s own pinned reference, fixed at exactly 0 by
construction rather than measured — gets shrinkage weight 1 (trusted
exactly, never pulled) and is excluded from estimating `tau2`; dividing
by its zero variance would otherwise blow up the whole calculation.
Among the other three drugs, `B` and `C`’s bias estimates already agree
with 0 and `D`’s disagrees sharply, so the fitted between-drug variance
is large relative to each drug’s own sampling variance — weights come
out close to 1, and `D`’s debiased estimate above moves only slightly
from the unshrunk `debiased` answer, not back toward its still-biased
raw value.

The real judgment call is the target. `shrink = "zero"` assumes negative
controls should show no systematic bias for a typical drug.
`shrink = "mean"` shrinks toward the panel’s own precision-weighted
grand mean instead — appropriate if some shared baseline bias (a
database- or design-level effect common to every drug) is expected, and
only *deviations* from that baseline count as real drug-specific signal:

``` r

shrunk_mean <- debias_surface(fit, bias_fit, shrink = "mean")
#> Warning: 1 of 5 drug(s) in `fit` have no bias estimate in `bias_fit` and are
#> dropped: E
shrunk_mean$shrinkage$target
#> [1] -0.1
shrunk_mean$effects[, c("drug", "estimate", "std_error")]
#>   drug    estimate  std_error
#> 1    A 0.000000000 0.00000000
#> 2    B 0.004166667 0.04947643
#> 3    C 0.204166667 0.04947643
#> 4    D 0.091666667 0.04947643
```

The grand mean is pulled negative by `D`’s large bias, so `"mean"`
treats a small negative bias as the *expected* baseline for every drug,
not just `D`. `B` and `C`’s debiased estimates tick up slightly relative
to `"zero"` even though their own negative controls showed no bias at
all — some of `D`’s problem gets diluted across the whole panel. (`A`,
the pinned reference, is untouched either way: its shrink weight is
always 1.) That is the cost of assuming a shared baseline: if the
panel’s heterogeneity is really just one bad drug rather than a common
database effect, `"mean"` quietly launders part of `D`’s bias into
everyone else’s estimate. Neither target is a safe default; which one is
right depends on what actually generated the panel’s bias, not on which
fits better.

## Debiasing reaches only as far as the bias network does

[`debias_surface()`](https://ablack3.github.io/directeffect/reference/debias_surface.md)
warns and drops any drug present in `fit` but absent from `bias_fit` —
here, `E`:

``` r

collect_warnings <- function(expr) {
  messages <- character()
  withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  messages
}
collect_warnings(debias_surface(fit, bias_fit))
#> [1] "1 of 5 drug(s) in `fit` have no bias estimate in `bias_fit` and are dropped: E"
```

`E`’s own trial against `D` may be perfectly clean, but `E`’s *surface*
position is still entangled with `D`’s through the network fit — and
with no negative-control edge of its own, there is nothing to correct it
with. Dropping `E` rather than reporting it half-corrected is
deliberate: a network-propagated bias estimate is only as complete as
the bias network it was propagated through.

## What debiasing does not fix

The debiased variance is `var(fit) + var(bias_fit)`, which assumes the
two fits’ sampling errors are independent (`independent = FALSE` is not
yet implemented and errors). In practice the outcome-of-interest fit and
the negative-control fit usually come from the same cohorts, the same
design, the same database — correlated error sources this simple sum
does not capture. The documented consequence is that debiased intervals
are, if anything, somewhat too narrow, not too wide. Debiasing corrects
a *point estimate* bias with a network-aware method; it is not a
substitute for getting the negative-control panel itself right.

## Learn more

- [`?debias_surface`](https://ablack3.github.io/directeffect/reference/debias_surface.md)
  — the full argument and return-value reference.
- [`vignette("pooling-vs-surface")`](https://ablack3.github.io/directeffect/articles/pooling-vs-surface.md)
  — a different correction problem (wrong estimand, not confounding)
  that
  [`debias_surface()`](https://ablack3.github.io/directeffect/reference/debias_surface.md)
  does not address.
