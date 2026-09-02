# Changelog

## directeffect 0.1.0

First tracked release.

### New features

- New
  [`debias_surface()`](https://ablack3.github.io/directeffect/reference/debias_surface.md):
  corrects a fitted surface using a second surface fit on the same
  comparison network, but with the outcome replaced by something whose
  true effect is known to be null for every drug (typically each edge’s
  own negative-control panel mean). Every drug’s position on that second
  (“bias”) surface is a network-propagated estimate of its systematic
  bias – connectivity-aware, unlike a raw per-drug average of edge-level
  bias – and
  [`debias_surface()`](https://ablack3.github.io/directeffect/reference/debias_surface.md)
  subtracts it from the surface of interest, summing variances under an
  explicit (documented, not silently assumed) independence
  approximation.

  Motivated by a network where inverse-variance pooling in
  [`fit_surface()`](https://ablack3.github.io/directeffect/reference/fit_surface.md)
  concentrated a large, direction-consistent per-edge bias into a
  spuriously narrow, “significant” surface position for one drug (its
  edges had per-edge empirical calibration already applied, but the
  systematic error, EASE, from that calibration was never propagated
  into the pooling weights) – the fix was to calibrate the surface
  itself, not just its edges.
  [`debias_surface()`](https://ablack3.github.io/directeffect/reference/debias_surface.md)
  generalizes that fix as reusable package functionality: fit the same
  network on each edge’s negative-control panel, and subtract the
  resulting bias surface.

### Prior work carried into this release

- [`pool_bootstrap()`](https://ablack3.github.io/directeffect/reference/pool_bootstrap.md)
  /
  [`pool_meta()`](https://ablack3.github.io/directeffect/reference/pool_meta.md):
  bootstrapped-mean and meta-analytic pooling of comparative estimates,
  as alternatives to the fitted surface, with
  [`compare_pooling()`](https://ablack3.github.io/directeffect/reference/compare_pooling.md)
  to compare them
  ([\#1](https://github.com/ablack3/directeffect/issues/1)).
- [`anchor_surface()`](https://ablack3.github.io/directeffect/reference/anchor_surface.md),
  [`fit_surface()`](https://ablack3.github.io/directeffect/reference/fit_surface.md)
  (netmeta and Stan engines),
  [`direct_effect_network()`](https://ablack3.github.io/directeffect/reference/direct_effect_network.md),
  connectivity/diagnostics
  ([`check_connectivity()`](https://ablack3.github.io/directeffect/reference/check_connectivity.md),
  [`edge_residuals()`](https://ablack3.github.io/directeffect/reference/edge_residuals.md),
  [`cycle_consistency()`](https://ablack3.github.io/directeffect/reference/cycle_consistency.md)),
  and the netmeta-vs-Stan
  [`compare_engines()`](https://ablack3.github.io/directeffect/reference/compare_engines.md)
  – the package’s original surface estimation contract.
