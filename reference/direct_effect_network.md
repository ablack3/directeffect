# Construct a direct-effect network

Builds a validated evidence network from a table of comparative
(drug-versus-drug) effect estimates, optionally alongside a table of
absolute (placebo-anchored) estimates. The returned object defines the
evidence network — it contains no fitted model. Estimates for
multiplicative effect measures (HR, RR, OR) must be supplied on the log
scale.

## Usage

``` r
direct_effect_network(comparisons, anchors = NULL, effect_measure = "HR")
```

## Arguments

- comparisons:

  A data frame with one row per estimated drug-versus-drug effect.
  Required columns: `study_id`, `target`, `comparator`, `estimate`,
  `std_error`. Additional metadata columns (e.g. `database`,
  `population`, `design`, `outcome`) are preserved.

- anchors:

  Optional data frame of absolute (placebo-controlled) estimates.
  Required columns: `study_id`, `drug`, `reference` (normally
  `"placebo"`), `estimate`, `std_error`. Additional columns are
  preserved.

- effect_measure:

  A single string naming the effect measure the estimates are on, e.g.
  `"HR"`, `"RR"`, or `"OR"`.

## Value

An object of class `directeffect_network` with components `comparisons`,
`anchors`, `treatments`, `graph` (an igraph graph), `components`
(integer component membership named by treatment), `effect_measure`, and
`metadata`.

## See also

[directeffect_formats](https://ablack3.github.io/directeffect/reference/directeffect_formats.md)
for the explicit input and output format reference, including validation
rules.

## Examples

``` r
comparisons <- data.frame(
  study_id   = c("S1", "S2", "S3"),
  target     = c("A", "A", "B"),
  comparator = c("B", "C", "C"),
  estimate   = c(log(1.02), log(1.34), log(1.29)),
  std_error  = c(0.07, 0.09, 0.08)
)
anchors <- data.frame(
  study_id  = "RCT1",
  drug      = "C",
  reference = "placebo",
  estimate  = log(1.20),
  std_error = 0.04
)
de <- direct_effect_network(comparisons, anchors = anchors,
                            effect_measure = "HR")
de
#> <directeffect_network>
#>   Effect measure: HR (log scale)
#>   Drugs:          3
#>   Comparisons:    3
#>   Anchors:        1
#>   Components:     1
#> Use check_connectivity() for the per-component identifiability report.
```
