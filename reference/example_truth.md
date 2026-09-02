# Generating truth of the example network

The true direct effects the example data were simulated from, on the log
hazard ratio scale with `placebo = 0`. These values are **inventions**
calibrated to the magnitude of the statin literature (hazard ratios
roughly 0.7–0.85 versus placebo), not estimates of any real drug's
effect. Shipping the truth makes the example honest: recovery of these
values by
[`fit_surface()`](https://ablack3.github.io/directeffect/reference/fit_surface.md) +
[`anchor_surface()`](https://ablack3.github.io/directeffect/reference/anchor_surface.md)
on
[example_comparisons](https://ablack3.github.io/directeffect/reference/example_comparisons.md)
and
[example_anchors](https://ablack3.github.io/directeffect/reference/example_anchors.md)
is verified in the test suite, not assumed.

## Usage

``` r
example_truth
```

## Format

A data frame with 7 rows (placebo plus six statins) and 2 columns:

- drug:

  Drug name; the first row is `"placebo"`.

- theta:

  True direct effect versus placebo, log hazard ratio scale (placebo row
  is exactly 0).

## Source

Simulated; see `data-raw/example-network.R`.
