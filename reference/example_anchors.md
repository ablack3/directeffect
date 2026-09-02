# Example placebo anchors: simulated landmark trials

Three simulated placebo-controlled estimates for statins in
[example_comparisons](https://ablack3.github.io/directeffect/reference/example_comparisons.md),
on the log hazard ratio scale — the pattern of the real statin
literature, where a few landmark placebo trials carry the absolute
information while most comparisons are head-to-head. Anchors position
the relative surface against placebo = 0. Like the comparisons, the
numbers are **simulated**, not real trial results.

## Usage

``` r
example_anchors
```

## Format

A data frame with 3 rows and 5 columns, in the anchors input format
documented in
[directeffect_formats](https://ablack3.github.io/directeffect/reference/directeffect_formats.md):

- study_id:

  Simulated trial identifier.

- drug:

  The anchored statin.

- reference:

  Always `"placebo"`.

- estimate:

  Simulated log hazard ratio, drug vs placebo.

- std_error:

  Its standard error.

## Source

Simulated; see `data-raw/example-network.R`. The generating truth ships
as
[example_truth](https://ablack3.github.io/directeffect/reference/example_truth.md).

## See also

[directeffect_formats](https://ablack3.github.io/directeffect/reference/directeffect_formats.md)
for the full format reference.
