# Simulate a direct-effect network with known truth

Generates true direct effects, a connected comparison graph, standard
errors, observed comparative estimates, and placebo anchor estimates.
The returned object retains the truth so recovery can be measured rather
than assumed — see
[`validate_recovery()`](https://ablack3.github.io/directeffect/reference/validate_recovery.md).

## Usage

``` r
simulate_direct_effect_network(
  n_drugs = 20,
  n_comparisons = 100,
  n_anchors = 3,
  heterogeneity = 0.1,
  seed = 1,
  effect_sd = 0.5,
  se_range = c(0.05, 0.15),
  anchor_se_range = c(0.03, 0.1)
)
```

## Arguments

- n_drugs:

  Number of drugs.

- n_comparisons:

  Number of comparative estimates; at least `n_drugs - 1` so the network
  can be connected.

- n_anchors:

  Number of placebo anchors, each on a distinct drug (at most
  `n_drugs`).

- heterogeneity:

  Between-study standard deviation `tau` added to every comparison's
  sampling variance. 0 matches the v0.1 common-effect model exactly.

- seed:

  Seed for reproducibility. The global random-number state is restored
  afterwards.

- effect_sd:

  Standard deviation of the true direct effects around 0 on the log
  scale.

- se_range:

  Range the comparison standard errors are drawn from, uniformly.

- anchor_se_range:

  Range the anchor standard errors are drawn from, uniformly.

## Value

A list with components `network` (a ready-to-fit
`directeffect_network`), `comparisons`, `anchors`, `truth` (data frame
of `drug`, `theta`, with the placebo row at 0), `tau`, and `seed`.

## Details

The generative model matches the package's fitting model: observed
comparisons are `y ~ N(theta_target - theta_comparator, se^2 + tau^2)`
where `tau` is `heterogeneity` (0 gives exact model match), and anchors
are `a ~ N(theta_drug, se^2)` against `reference = "placebo"` at 0. The
comparison graph is connected by construction: a random spanning tree
first, then the remaining comparisons between random drug pairs (repeat
comparisons of the same pair are allowed, as in real evidence networks).

## Examples

``` r
simulation <- simulate_direct_effect_network(
  n_drugs = 5, n_comparisons = 10, n_anchors = 1,
  heterogeneity = 0, seed = 1
)
simulation$truth
#>      drug       theta
#> 1 placebo  0.00000000
#> 2 drug_01 -0.31322691
#> 3 drug_02  0.09182166
#> 4 drug_03 -0.41781431
#> 5 drug_04  0.79764040
#> 6 drug_05  0.16475389
```
