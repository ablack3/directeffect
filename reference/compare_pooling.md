# Compare bootstrapped pooling with the fitted surface

Runs
[`pool_bootstrap()`](https://ablack3.github.io/directeffect/reference/pool_bootstrap.md)
and lines its per-drug answer up against a
[`fit_surface()`](https://ablack3.github.io/directeffect/reference/fit_surface.md)
fit of the same network, on the same reference scale.

## Usage

``` r
compare_pooling(de, fit = NULL, ...)
```

## Arguments

- de:

  A `directeffect_network`.

- fit:

  A `directeffect_fit` from
  [`fit_surface()`](https://ablack3.github.io/directeffect/reference/fit_surface.md)
  on the same network. If `NULL`, one is fitted with the `netmeta`
  engine.

- ...:

  Passed to
  [`pool_bootstrap()`](https://ablack3.github.io/directeffect/reference/pool_bootstrap.md).

## Value

A data frame with one row per target drug: the pooled estimate, the
surface estimate, `comparator_offset`, and the difference in log HR.

## Details

The two are put on a common footing by re-expressing the surface
relative to each drug's own comparators: a drug's surface position minus
the evidence-weighted mean position of the drugs it was compared
against. That difference is what pooling implicitly assumes is zero – it
treats every comparator as if it sat at the network average.
`comparator_offset` reports how wrong that assumption is for each drug,
and is the quantity that explains where the two methods part company.
