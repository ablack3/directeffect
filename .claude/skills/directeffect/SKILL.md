---
name: directeffect
description: Use the directeffect R package to estimate latent direct (vs-placebo) drug effects from a network of comparative-effect estimates. Use when the user wants to run directeffect, build a directeffect_network, fit or anchor a surface, triangulate placebo-anchored effects from head-to-head comparisons, or otherwise work with this repo's core R package (fit_surface, anchor_surface, direct_effect_network, compare_engines, check_connectivity, edge_residuals, plot_network, plot_sea_level).
---

# directeffect

This repo *is* the `directeffect` R package. It reconstructs each drug's
latent **direct** effect (vs placebo) from a network of **comparative**
(drug-vs-drug) effect estimates — most comparative-effectiveness evidence
only ever measures differences between drugs, not any drug's absolute
effect, so a null A-vs-B result does not mean either drug is safe; it only
means the two sit at the same height.

The method is deliberately two steps, never conflated:

1. **The surface** — `fit_surface()` estimates every drug's position
   *relative* to the others from comparisons alone (anchors are ignored at
   this step). Tells you whether the comparative evidence is internally
   coherent and what relative structure it implies.
2. **The sea level** — `anchor_surface()` uses placebo-controlled anchor
   estimates to position that whole surface against placebo = 0, carrying
   the anchors' own uncertainty forward. Tells you where the structure
   sits in absolute terms.

Two independent engines sit behind one interface and are cross-validated
against each other in the test suite: `"netmeta"` (frequentist, handles
multi-arm trials via `study_id`) and `"stan"` (Bayesian via rstan; refuses
multi-arm networks rather than mistreat correlated rows as independent).

## Data contract — everything is log scale

All estimates enter and leave the package on the **log scale** for
multiplicative measures (HR, RR, OR) — pass `log(1.02)`, never `1.02`.
Full validation rules: `?directeffect_formats`.

**Comparisons** (one row per drug-vs-drug estimate, required):
`study_id`, `target`, `comparator` (≠ target), `estimate` (log scale),
`std_error` (> 0). Extra columns are preserved; estimand-describing
columns that differ across rows for the same edge (`population`,
`database`, `design`, …) trigger a transportability warning — investigate
rather than silence it.

**Anchors** (optional, one row per placebo-controlled estimate):
`study_id`, `drug` (must already appear in a comparison), `reference`
(normally `"placebo"`), `estimate`, `std_error`.

**Output — `fit$effects`** (both engines, anchored or not): one row per
drug — `drug`, `estimate`, `std_error`, `lower`, `upper`, `scale` (always
`"log"`), `reference` (a fixed drug, or `"placebo"` once anchored),
`engine`. Stan fits additionally carry `median`, `mean`, `sd`, `q025`,
`q975`, `rhat`, `ess_bulk`, `ess_tail` — check `rhat`/`ess_*` before
trusting a Stan fit.

## Standard workflow

```r
library(directeffect)   # or devtools::load_all() inside the repo

# 1. Build a validated evidence network — no model fitted yet.
de <- direct_effect_network(
  comparisons,           # target/comparator/estimate/std_error rows
  anchors = anchors,      # optional placebo-controlled rows
  effect_measure = "HR"   # or "RR", "OR" — informs display, not the math
)
de                        # print method summarizes drugs/comparisons/components

# 2. Check identifiability BEFORE fitting — disconnected components can't
#    be compared to each other at all.
check_connectivity(de)

# 3. Fit the relative surface, then anchor it as a separate, explicit step.
fit <- fit_surface(de, engine = "netmeta")   # or engine = "stan"
absolute <- anchor_surface(fit)
absolute$effects

# 4. Diagnostics and plots.
edge_residuals(fit)        # leverage-aware residuals — coherence of the surface
cycle_consistency(de)      # closed-loop consistency check
plot_network(de)
plot_sea_level(absolute)
plot_effect_surface(fit)

# 5. Cross-validate engines when it matters (single-comparison-per-study
#    networks only — multi-arm studies are netmeta-only).
compare_engines(de)
```

## Rules that matter

- Never skip step 1→2: don't call `fit_surface()` on a network you haven't
  run `check_connectivity()` on first — a disconnected network silently
  gives you effects with no real meaning across components.
- Never anchor before fitting, and never let anchors leak into
  `fit_surface()` — the surface must be estimated from comparisons alone;
  anchoring is a strictly separate positioning step.
- Convert real-world estimates (HRs, RRs, ORs, and their CIs) to log scale
  before building the network — `estimate = log(hr)`, and derive
  `std_error` from the log-scale CI half-width (e.g.
  `(log(upper) - log(lower)) / (2 * 1.96)`), not from the linear-scale CI.
- `engine = "stan"` will refuse a network with multi-arm studies (more
  than one comparison sharing a `study_id`) — use `"netmeta"` for those,
  or split/exclude the offending studies if you need the Stan engine.
- Treat any transportability warning as a data-quality signal, not a
  formality — it means the same edge is being asserted from meaningfully
  different populations/designs/databases.

## Learn more in this repo

- `README.md` — the pitch and the full example on `example_comparisons` /
  `example_anchors` / `example_truth` (simulated statin network).
- `vignette("directeffect")` — getting-started walkthrough.
- `vignette("how-it-works")` — the method from scratch, graded against
  known simulated truth.
- `vignette("netmeta-vs-stan")` — same surface, both engines.
- `vignette("validating-the-surface")` — edge residuals, cycle
  consistency, simulation-based recovery.
- `?directeffect_formats` — the authoritative input/output column
  reference and validation rules.
- `.scratch/v0.1-core/` — spec and ticket breakdown for the current
  version.

## Working in real (non-simulated) applications

When applying this to a real drug-safety/effectiveness sweep (e.g. a
network of pairwise comparative estimates for many drugs against a
shared outcome), private/patient-level source data and any derived
results belong under `private/` in this repo (gitignored) — never commit
real study data or results to git, only the package code and simulated
example data.
