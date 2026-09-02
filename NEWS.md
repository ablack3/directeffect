# directeffect 0.1.0

First tracked release.

## New features

* New `debias_surface()`: corrects a fitted surface using a second surface
  fit on the same comparison network, but with the outcome replaced by
  something whose true effect is known to be null for every drug (typically
  each edge's own negative-control panel mean). Every drug's position on
  that second ("bias") surface is a network-propagated estimate of its
  systematic bias -- connectivity-aware, unlike a raw per-drug average of
  edge-level bias -- and `debias_surface()` subtracts it from the surface of
  interest, summing variances under an explicit (documented, not silently
  assumed) independence approximation.

  Motivated by a network where inverse-variance pooling in `fit_surface()`
  concentrated a large, direction-consistent per-edge bias into a
  spuriously narrow, "significant" surface position for one drug (its
  edges had per-edge empirical calibration already applied, but the
  systematic error, EASE, from that calibration was never propagated into
  the pooling weights) -- the fix was to calibrate the surface itself, not
  just its edges. `debias_surface()` generalizes that fix as reusable
  package functionality: fit the same network on each edge's
  negative-control panel, and subtract the resulting bias surface.

## Prior work carried into this release

* `pool_bootstrap()` / `pool_meta()`: bootstrapped-mean and meta-analytic
  pooling of comparative estimates, as alternatives to the fitted surface,
  with `compare_pooling()` to compare them (#1).
* `anchor_surface()`, `fit_surface()` (netmeta and Stan engines),
  `direct_effect_network()`, connectivity/diagnostics
  (`check_connectivity()`, `edge_residuals()`, `cycle_consistency()`), and
  the netmeta-vs-Stan `compare_engines()` -- the package's original surface
  estimation contract.
