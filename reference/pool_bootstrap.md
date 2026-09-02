# Pool comparative estimates by bootstrapped mean across trials

Averages a drug's comparative estimates directly, with a non-parametric
bootstrap over the estimates for the interval. This is the pooling used
in high-throughput target-trial emulation (Zang et al., *Nature
Communications* 2023;14:8180): per drug, take the sample mean of the
adjusted hazard ratios over the trials that balanced, bootstrap the
trials 1,000 times for a 95% interval, and test whether that mean is
below 1.

## Usage

``` r
pool_bootstrap(
  de,
  balance = NULL,
  max_smd = 0.1,
  min_balanced_frac = 0.1,
  R = 1000,
  scale = c("hr", "log"),
  seed = NULL
)
```

## Arguments

- de:

  A `directeffect_network` from
  [`direct_effect_network()`](https://ablack3.github.io/directeffect/reference/direct_effect_network.md),
  or a plain data frame of comparisons with `target`, `comparator` and
  `estimate`.

- balance:

  Optional data frame with `target`, `comparator` and `max_abs_smd`,
  used to keep only balanced trials. The source paper counts a trial as
  balanced when at most 2% of covariates exceed \|SMD\| 0.1; supply
  whatever balance summary the upstream pipeline produced and set
  `max_smd` to the threshold.

- max_smd:

  Balance threshold applied to `balance$max_abs_smd`. Trials above it
  are dropped. Ignored when `balance` is `NULL`.

- min_balanced_frac:

  A drug is reported only if at least this fraction of its trials
  survive the balance filter. The source paper uses 0.10. Ignored when
  `balance` is `NULL`.

- R:

  Bootstrap resamples. Default 1000, as in the source paper.

- scale:

  `"hr"` averages on the hazard-ratio (linear) scale, which is what the
  source paper does. `"log"` averages log hazard ratios and
  exponentiates, which is the scale the estimates are actually symmetric
  on and is usually the better choice; it is not what the paper did.

- seed:

  Optional integer seed. The caller's RNG state is restored on exit, so
  pooling never disturbs the session's reproducibility.

## Value

A data frame with one row per drug that appears as a target: `drug`,
`n_trials`, `n_balanced`, `frac_balanced`, `pooled_hr`, `ci_lower`,
`ci_upper`, `p_value` (bootstrap, one-sided against HR \< 1, as in the
source paper), `scale`, and `reported` (whether it cleared
`min_balanced_frac`).

## Details

It is offered here as a deliberate contrast to
[`fit_surface()`](https://ablack3.github.io/directeffect/reference/fit_surface.md),
not as an alternative implementation of it. The two answer different
questions:

- `pool_bootstrap()` asks **"how did this drug fare against the
  comparators it happened to be tested against?"** Each estimate enters
  as-is, so the answer inherits whatever the comparators' own effects
  were. A drug compared mostly against harmful drugs looks protective,
  and nothing in the method can tell that apart from the drug being
  protective.

- [`fit_surface()`](https://ablack3.github.io/directeffect/reference/fit_surface.md)
  asks **"where does this drug sit relative to every other drug?"**,
  solving the whole network at once so a comparator's own position is
  estimated rather than assumed to be neutral.

Where the two disagree, the disagreement is informative: it measures how
far a drug's comparator set sits from the network average;
[`compare_pooling()`](https://ablack3.github.io/directeffect/reference/compare_pooling.md)
reports that quantity per drug.

## References

Zang C, Zhang H, Xu J, et al. High-throughput target trial emulation for
Alzheimer's disease drug repurposing with real-world data. *Nature
Communications* 2023;14:8180.

## See also

[`fit_surface()`](https://ablack3.github.io/directeffect/reference/fit_surface.md)
for the network alternative, and
[`compare_pooling()`](https://ablack3.github.io/directeffect/reference/compare_pooling.md)
to run both on one network and quantify the gap.

## Examples

``` r
de <- direct_effect_network(example_comparisons, effect_measure = "HR")
pool_bootstrap(de, R = 200, seed = 1)
#>           drug n_trials n_balanced frac_balanced pooled_hr  ci_lower  ci_upper
#> 1 atorvastatin        4          4             1 0.9254846 0.8572401 1.0521182
#> 2   lovastatin        1          1             1 1.0140985 1.0140985 1.0140985
#> 3  pravastatin        2          2             1 0.8922181 0.8130196 0.9714165
#> 4 rosuvastatin        3          3             1 0.9055725 0.8478937 0.9512294
#> 5  simvastatin        2          2             1 0.8840177 0.8504412 0.9175942
#>   p_value scale reported
#> 1   0.075    hr     TRUE
#> 2   1.000    hr     TRUE
#> 3   0.000    hr     TRUE
#> 4   0.000    hr     TRUE
#> 5   0.000    hr     TRUE
```
