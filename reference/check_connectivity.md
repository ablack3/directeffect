# Report per-component identifiability of a direct-effect network

Comparative estimates identify only differences within a connected
component; a component with no absolute anchor supports relative effects
only, and the package never implies such a component can be positioned
absolutely. `check_connectivity()` reports, for each connected
component, the drug count, comparison count, anchor count, and whether
absolute effects are identifiable there.

## Usage

``` r
check_connectivity(de)
```

## Arguments

- de:

  A `directeffect_network` created by
  [`direct_effect_network()`](https://ablack3.github.io/directeffect/reference/direct_effect_network.md).

## Value

A data frame of class `directeffect_connectivity` with one row per
connected component and columns `component`, `n_drugs`, `n_comparisons`,
`n_anchors`, and `absolute_identifiable`. The `drugs` attribute holds
the drug names per component. Printing renders the identifiability
report.

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
check_connectivity(de)
#> Component 1:
#>   3 drugs
#>   3 comparisons
#>   0 absolute anchors
#>   relative effects identifiable only
```
