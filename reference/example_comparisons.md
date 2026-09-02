# Example head-to-head comparisons: a simulated statin network

A realistic but **entirely simulated** evidence base modeled on the
statin trial literature: twelve head-to-head estimates among six statins
(atorvastatin, fluvastatin, lovastatin, pravastatin, rosuvastatin,
simvastatin) for a major-cardiovascular-events outcome, on the log
hazard ratio scale. The drug names are real; every number is fake —
drawn from the invented truth in
[example_truth](https://ablack3.github.io/directeffect/reference/example_truth.md)
— and must not be interpreted clinically. Together with
[example_anchors](https://ablack3.github.io/directeffect/reference/example_anchors.md)
it powers the package vignettes. Generated deterministically by
`data-raw/example-network.R` with mild between-study heterogeneity (tau
= 0.02); the test suite verifies the shipped data matches that generator
and that the workflow recovers
[example_truth](https://ablack3.github.io/directeffect/reference/example_truth.md).

## Usage

``` r
example_comparisons
```

## Format

A data frame with 12 rows and 5 columns, in the comparisons input format
documented in
[directeffect_formats](https://ablack3.github.io/directeffect/reference/directeffect_formats.md):

- study_id:

  Simulated trial identifier.

- target, comparator:

  The two statins compared.

- estimate:

  Simulated log hazard ratio, target vs comparator.

- std_error:

  Its standard error.

## Source

Simulated; see `data-raw/example-network.R`. The generating truth ships
as
[example_truth](https://ablack3.github.io/directeffect/reference/example_truth.md).

## See also

[directeffect_formats](https://ablack3.github.io/directeffect/reference/directeffect_formats.md)
for the full format reference.
