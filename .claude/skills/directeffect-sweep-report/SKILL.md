---
name: directeffect-sweep-report
description: Run the directeffect surface fit + OHDSI EmpiricalCalibration on a pragma/OHDSI drug-sweep output (target-vs-comparator Cox HRs with negative/positive controls, living under private/), and produce a standardized Excel workbook plus a self-contained HTML report. Use when the user asks to "run directeffect on a sweep," "calibrate and report a sweep," "produce the xlsx/report for <sweep>," or points at a pragma sweep's output data and wants triangulated, calibrated effect estimates written up.
---

# directeffect sweep → calibrated report

This is the repeatable version of the one-off pipeline built by hand for
`private/koa_tkr_sweep_FINAL_lean/` (see that folder's `SWEEP_OUTPUT_SPEC.md` for the
original data-format investigation this skill generalizes, and the `directeffect` skill for
the underlying package API). There is **no package function or CLI for this yet** — every
run means writing/adapting R scripts per this spec. If you're doing this for the third time
on the third sweep, that's the signal to actually build `report_sweep()` into the package
(see "Should this become a package function?" at the end) rather than copy-pasting again.

**Never commit any of this.** Sweep input data, generated `comparisons_used.csv`, `.rds`
fits, plots, and HTML reports are all real (or synthetic-but-sensitive) study data. Sweep
inputs live under `private/` in `directeffect` (gitignored); write ALL outputs there too
(a fresh subfolder such as `private/<sweep_name>/report/`), never elsewhere in the repo, and
never `git add` any of it. Per the user's own pattern in this session, generated report
bundles typically get moved out to `~/projects/rwe-dwas-results/` or similar afterward —
ask, don't assume, unless already instructed.

## 1. Detect the input format first — don't assume

Sweep output shows up in (at least) two incompatible shapes. Read a few real files before
writing any code:

**Shape A — xlsx + per-arm folder tree** (the `koa_tkr_sweep_FINAL_lean` shape):
- `<sweep>/<Name>_Sweep_Results.xlsx`, sheet `"All Results"`: one row per arm, columns
  `target`, `arm` (`armNN` = primary/pre-specified, `expNNN`/`chaoNNN`/`deepNNN` =
  exploratory), `cox_hr`/`cox_hr_lo`/`cox_hr_hi`, `shr`/`shr_lo`/`shr_hi`, `gate2_pass`,
  `outcome`, `population`, `study_design`.
- `<sweep>/results/<target>/<armNN_or_expNNN>_<comparator>/` per-arm folder, containing
  `cohort_metadata.csv`, `cox_unadjusted_results.csv`, `fine_gray_results.csv`,
  `competing_risks_summary.csv`, `covariate_balance.csv`, `negative_control_irr.csv`,
  `negative_control_null_summary.csv`, and (if run) `positive_control_irr.csv`.

**Shape B — flat comparisons/controls CSVs** (the `pragma` grid-engine shape, e.g.
`pragma/sweeps/koa_tkr_2026-08b/`):
- `comparisons.csv` — one row per target-vs-comparator pair, already close to
  `directeffect`'s contract.
- `negative_controls.csv` — the negative-control set(s), shared/referenced rather than
  duplicated per-arm.
- No per-arm folder tree; calibration inputs come from wherever the grid's calibration
  export lands (check `HANDOFF.md`/`README.md` in that specific sweep directory — this
  shape has changed across `pragma` sweeps, don't assume the `koa_tkr_2026-08b` layout is
  current).

**Before running anything**: `head`/`read.csv` the actual files in the target sweep
directory and confirm which shape (or a third, undocumented one) you're looking at. If it's
neither A nor B, stop and describe the format back to the user rather than guessing a
mapping.

## 2. Build the `directeffect` comparisons table

Target contract (`?directeffect_formats`): `study_id`, `target`, `comparator`, `estimate`
(log scale), `std_error`. From Shape A:

```r
comparisons$estimate  <- log(comparisons$cox_hr)
comparisons$std_error <- (log(comparisons$cox_hr_hi) - log(comparisons$cox_hr_lo)) /
                          (2 * qnorm(0.975))
```

From Shape B, `comparisons.csv` may already carry `estimate`/`std_error` in log scale —
check before re-deriving.

**Judgment calls to make explicitly, every time, and document in the run script's
comments** (do not silently default):
- Which arms/rows form the network edges? (e.g. Shape A: primary `armNN` rows only, unless
  the user wants exploratory arms included too — ask if it's not obvious.)
- Is there a placebo/vs-nothing anchor anywhere in this sweep? If yes, use
  `anchor_surface()`. If no (the common case for active-comparator sweeps), do **not**
  fabricate one — fit relative to a reference drug only (pick the most-connected drug,
  document why) and say explicitly in the report that these are relative, not absolute,
  effects.
- Multi-arm studies (shared `study_id` across >1 comparison) force `engine = "netmeta"` —
  `"stan"` will refuse them.

## 3. Calibrate with OHDSI `EmpiricalCalibration`

Requires the package installed (source at `private/EmpiricalCalibration/` if not already
on CRAN-installed in this environment: `devtools::install("private/EmpiricalCalibration")`).

Per arm/comparison:
1. Load that arm's negative controls (Shape A: `negative_control_irr.csv`, filtered to
   `control_type == "negative"`; Shape B: wherever that sweep's calibration export lives).
   Drop non-estimable rows (`NA`/zero/infinite `irr`, `ci_lo`, or `ci_hi`).
2. `logRr = log(irr)`, `seLogRr = (log(ci_hi) - log(ci_lo)) / (2*qnorm(0.975))`.
3. Require a minimum estimable-control count (5 is what was used for `koa_tkr_sweep_FINAL_lean`
   — adjust and document if this sweep's controls are sparser/denser) before fitting; skip
   and flag (`calibrated = FALSE`, `skip_reason`) rather than crash on arms below it.
4. `EmpiricalCalibration::fitNull()` on the negative controls →
   `convertNullToErrorModel()`.
5. **If `positive_control_irr.csv` (or equivalent) exists for this sweep**: fit the
   systematic-error model's slope from both negative AND positive controls together
   (check `EmpiricalCalibration`'s current function signatures for this — do not assume
   slope=1 when positive controls are actually available; that shortcut is only correct
   when they aren't).
6. `calibrateConfidenceInterval()` on that arm's own `estimate`/`std_error`.

Output: `raw_estimate`, `raw_std_error`, `calibrated_estimate`, `calibrated_std_error`,
`calibrated` (bool), `n_negative_controls_used`, `n_positive_controls_used` (0 if none),
`skip_reason` — added as new columns onto the comparisons table, not replacing the raw
columns (both need to survive into the report).

## 4. Fit the surface (raw AND calibrated)

Run `direct_effect_network()` → `fit_surface(engine = "netmeta")` (→ `anchor_surface()` only
if a real placebo anchor exists) **twice**: once on raw estimates, once on calibrated —
same reference drug, same engine, same edge set, so the two `$effects` tables are directly
comparable row-for-row. This before/after pair is the core deliverable; don't skip the raw
fit just because calibration ran.

## 5. Excel workbook — standardized sheet layout

One `.xlsx`, written with `writexl::write_xlsx()` (add as a dependency if not present) or
`openxlsx` if formatting/conditional highlighting is wanted (e.g. shading rows where
significance flips between raw and calibrated). Sheets, in this order:

1. **`Summary`** — sweep name, date run, n drugs, n comparisons, engine, reference drug,
   whether a placebo anchor was used, calibration status (n arms calibrated / skipped),
   one row per headline caveat (no placebo anchor / no positive controls / etc.).
2. **`Surface (raw)`** — `fit_surface()$effects` on raw estimates: drug, estimate, std_error,
   lower, upper, scale, reference, engine.
3. **`Surface (calibrated)`** — same shape, calibrated estimates.
4. **`Comparisons used`** — the full comparisons table actually fed into both fits, with
   raw + calibrated columns from step 3, plus `gate2_pass`/equivalent calibration-gate flag
   from the source sweep if present.
5. **`Connectivity`** — `check_connectivity()` output (flag any disconnected components —
   these drugs' estimates aren't comparable to the rest of the network).
6. **`Edge residuals`** — `edge_residuals()` output, for spotting individual
   incoherent-with-the-network comparisons.
7. Optional **`Calibration detail`** sheet if calibration ran: per-arm null-fit parameters
   (`null_mean_log_irr`, `null_sd_log_irr`, n controls used, slope if fit from positive
   controls) — the numbers behind step 3, not just the final calibrated estimate.

Name the file `<sweep_name>_directeffect_report_<YYYY-MM-DD>.xlsx`.

## 6. HTML report — standardized structure

Self-contained single file (all images as base64 data URIs, no external assets, no CDN
dependencies — the report may be opened offline or on a machine with no internet), light/
dark-theme-aware CSS (see `results.html`/`results_calibrated.html` in
`private/koa_tkr_sweep_FINAL_lean/` — since moved to
`~/projects/rwe-dwas-results/directeffect_calibration_2026-08-29/` — as the concrete style
reference: title + subtitle, caveat callout boxes, summary stat cards, plots with captions,
detail tables, appendix). Section order, standardized across sweeps:

1. Title + one-line subtitle (sweep name, n drugs/comparisons, engine).
2. Caveat callout(s) — anchor status (placebo vs relative-only), calibration status,
   any other sweep-specific caveats — **always near the top, never buried**.
3. Summary stat cards (n drugs, n significant raw, n significant calibrated, n arms
   calibrated/skipped).
4. Network plot (`plot_network()`).
5. Effect surface plot — raw, then calibrated, side by side or stacked, clearly labeled.
6. Full drug ranking table — raw vs calibrated estimate + CI + significance, flipped rows
   highlighted.
7. Per-drug drill-down section for any drug the user specifically asked about (e.g. the
   roflumilast calibration-before-after table/plot pattern from the prior run) — optional,
   only when there's a specific drug of interest, not generated for all 40+ drugs by
   default.
8. Limitations section — spell out every judgment call from step 2, the calibration
   minimum-controls threshold and how many arms it excluded, whether positive controls were
   available, and what confounding sources calibration cannot detect.
9. Appendix — links to any external literature cited in a drill-down section (see the
   roflumilast biological-mechanism appendix as the template), and the exact list of source
   files this report was built from (traceability).

Name the file `<sweep_name>_report_<YYYY-MM-DD>.html`.

## 7. Constraints, every run

- All source data reads and output writes stay under `private/` inside `directeffect`
  (or wherever the sweep's own private data root is) — never git-tracked.
- Actually run the R code end to end (`devtools::load_all()` the `directeffect` package,
  execute for real, verify real output) — don't hand back a script that "should work."
- State explicitly, in both the xlsx `Summary` sheet and the HTML caveat section: whether
  effects are absolute (placebo-anchored) or relative-to-a-reference-drug, and whether
  calibration used positive controls or negative-only (slope assumed vs fit).
- If the sweep is large enough that Cox/PSM computation itself needs (re-)running — as
  opposed to just being read and reported on — that's a different, much bigger job (see
  `private/koa_tkr_sweep_FINAL_lean/SUMMARY.md`'s account of extending a sweep): check with
  the user about compute location (local machine vs AWS VM vs Snowflake) and resource
  headroom before doing that. This skill assumes the sweep's raw Cox/PSM/negative-control
  output already exists on disk and only the reporting/calibration layer needs running.

## Should this become a package function?

Right now this is a spec for hand-written scripts, on purpose — the two known sweep shapes
(§1) aren't even unified yet, and `EmpiricalCalibration` isn't a `directeffect` dependency
(see `private/SUMMARY.md`). If this skill gets invoked a third time on a third sweep with a
third data shape, that's the signal to stop re-deriving the mapping and instead: add
`EmpiricalCalibration` as a real `Suggests`/`Imports` in `DESCRIPTION`, write an exported
`report_sweep(sweep_dir, out_dir, shape = c("auto", "xlsx_arms", "flat_csv"))` function with
tests against fixtures of both shapes, and retire this skill in favor of calling that
function. Don't make that call unilaterally — it's a real API-design decision for the
package owner.
