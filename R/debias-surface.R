#' Debias a fitted surface using a network-level bias estimate
#'
#' Per-edge calibration (e.g. empirical calibration against a negative-control
#' panel) corrects each comparison in isolation. When a drug's own comparisons
#' are consistently biased in the same direction across many edges, inverse-
#' variance pooling in [fit_surface()] treats that consistency as *precision*
#' rather than *shared bias* — it concentrates the bias into a spuriously
#' narrow surface position instead of cancelling it out.
#'
#' `debias_surface()` corrects for this by taking a second surface fit on the
#' *same* comparison network, but with the outcome of interest replaced by
#' something whose true effect is known to be null for every drug — typically
#' the mean of each edge's own negative-control panel. Because every drug's
#' true position on that second surface is 0, any nonzero fitted position is
#' a network-propagated estimate of that drug's systematic bias: connectivity-
#' aware, unlike a raw per-drug average of edge-level bias, which ignores how
#' bias can concentrate through well- versus poorly-connected nodes. Subtracting
#' it from the surface of interest is a network-level analogue of empirical
#' calibration — corrected once, for the whole surface, rather than edge by edge.
#'
#' @param fit A `directeffect_fit` for the outcome of interest, from
#'   [fit_surface()].
#' @param bias_fit A second `directeffect_fit`, built by calling
#'   [direct_effect_network()] and [fit_surface()] on the same comparisons
#'   (same drugs, same edges) but with `estimate`/`std_error` replaced by each
#'   edge's own bias signal (e.g. the mean and standard error of the mean of
#'   its negative-control panel's log effect). Every drug's fitted position
#'   here is an estimate of its systematic bias, not a real effect.
#' @param independent Are `fit` and `bias_fit`'s sampling errors independent?
#'   Default `TRUE`, in which case `var(debiased) = var(fit) + var(bias_fit)`.
#'   This likely *understates* the true variance somewhat, since both fits
#'   typically come from the same underlying cohorts and design and so share
#'   correlated error sources that this simple sum does not capture — there is
#'   currently no general way to estimate the covariance between two
#'   independently-fit surfaces, so that understatement is a known limitation,
#'   not something this function corrects for. `independent = FALSE` is not
#'   yet implemented and errors.
#' @param shrink Shrink `bias_fit`'s per-drug estimates toward a common target
#'   before subtracting, rather than taking each at face value. `"none"`
#'   (default) subtracts `bias_fit` as-is. `"zero"` and `"mean"` both apply
#'   DerSimonian--Laird empirical-Bayes shrinkage — the same between-group
#'   variance estimator [pool_meta()] uses for `method = "random"`, but here
#'   estimated *across drugs'* bias estimates rather than across one drug's
#'   trials: drug `i`'s raw bias estimate is pulled toward the target by
#'   `tau2 / (tau2 + se_i^2)`, where `tau2` is the between-drug variance of
#'   true bias estimated by the method of moments. A drug with a large,
#'   imprecise network-propagated bias estimate (few or weak negative
#'   controls, reached mostly through connectivity rather than direct
#'   evidence) gets pulled hard toward the target; a drug with a precise one
#'   is barely moved. This is not automatic goodness — it is a real judgment
#'   call about what negative controls are telling you on average:
#'   * `"zero"` shrinks toward 0, i.e. assumes negative controls should show
#'     no systematic bias for a typical drug, so any nonzero panel-wide
#'     signal is itself noise to shrink away.
#'   * `"mean"` shrinks toward the precision-weighted grand mean of all
#'     drugs' bias estimates, i.e. assumes some shared baseline bias (a
#'     database- or design-level effect common to every drug) is real and
#'     expected, and only *deviations* from that baseline are drug-specific
#'     signal worth keeping.
#'
#'   `bias_fit`'s own pinned reference drug (`estimate`/`std_error` fixed at
#'   exactly 0 by [fit_surface()]'s normalization, not measured) carries no
#'   information about between-drug bias variance and is excluded from
#'   estimating `tau2`; it is given weight 1 (trusted exactly, never pulled)
#'   rather than divide by its zero variance. Requires at least 2 *other*
#'   drugs in common between `fit` and `bias_fit` -- `tau2` is not estimable
#'   from fewer, and this errors rather than silently treat `tau2` as 0
#'   (which would force full shrinkage on every other drug for no real
#'   reason). The shrunk bias covariance uses `Cov[i, j] * sqrt(w_i * w_j)`
#'   (`w` the shrinkage weights above) as a simple, diagonal-consistent
#'   scaling; like `independent`, this does not account for the extra
#'   estimation uncertainty in `tau2` and the target themselves, so treat
#'   shrunk intervals as an approximation, not an exact posterior.
#'
#' @return A `directeffect_fit` restricted to the drugs `fit` and `bias_fit`
#'   have in common (a warning names any dropped), with `estimate` equal to
#'   `fit`'s estimate minus `bias_fit`'s (optionally shrunk) estimate per drug
#'   and `covariance` equal to the (assumed-independent) sum. `$bias_fit` on
#'   the result carries the original, unshrunk bias surface, for provenance.
#'   When `shrink != "none"`, `$shrinkage` additionally carries a list with
#'   `method`, `tau2`, `target`, and `weight` (the per-drug shrinkage weights,
#'   named by drug).
#'
#' @examples
#' comparisons <- data.frame(
#'   study_id   = c("S1", "S2", "S3"),
#'   target     = c("A", "A", "B"),
#'   comparator = c("B", "C", "C"),
#'   estimate   = c(0.5, 0.9, 0.4),
#'   std_error  = c(0.05, 0.05, 0.05)
#' )
#' bias <- data.frame(
#'   study_id   = c("S1", "S2", "S3"),
#'   target     = c("A", "A", "B"),
#'   comparator = c("B", "C", "C"),
#'   estimate   = c(0.0, 0.2, 0.2),
#'   std_error  = c(0.05, 0.05, 0.05)
#' )
#' if (requireNamespace("netmeta", quietly = TRUE)) {
#'   fit <- fit_surface(direct_effect_network(comparisons), engine = "netmeta")
#'   bias_fit <- fit_surface(direct_effect_network(bias), engine = "netmeta")
#'   debiased <- debias_surface(fit, bias_fit)
#'   debiased$effects
#'
#'   # Shrink bias estimates toward 0 before subtracting, rather than trusting
#'   # each drug's raw network-propagated bias estimate at face value.
#'   shrunk <- debias_surface(fit, bias_fit, shrink = "zero")
#'   shrunk$effects
#'   shrunk$shrinkage
#' }
#' @export
debias_surface <- function(fit, bias_fit, independent = TRUE,
                            shrink = c("none", "zero", "mean")) {
  assert_directeffect_fit(fit)
  assert_directeffect_fit(bias_fit)
  assert_fit_covariance(fit, "debias_surface()")
  assert_fit_covariance(bias_fit, "debias_surface()")
  if (!isTRUE(independent)) {
    stop("`debias_surface()` currently only supports `independent = TRUE`; ",
         "there is no general way to estimate the covariance between two ",
         "independently-fit surfaces.", call. = FALSE)
  }
  shrink <- match.arg(shrink)

  fit_drugs <- fit$effects$drug
  bias_drugs <- bias_fit$effects$drug
  common <- intersect(fit_drugs, bias_drugs)
  if (length(common) == 0) {
    stop("`fit` and `bias_fit` share no drugs in common.", call. = FALSE)
  }
  missing_from_bias <- setdiff(fit_drugs, bias_drugs)
  if (length(missing_from_bias) > 0) {
    warning(
      length(missing_from_bias), " of ", length(fit_drugs),
      " drug(s) in `fit` have no bias estimate in `bias_fit` and are dropped: ",
      paste(utils::head(missing_from_bias, 5), collapse = ", "),
      if (length(missing_from_bias) > 5) ", ...",
      call. = FALSE
    )
  }
  fit_pos <- match(common, fit$effects$drug)
  bias_pos <- match(common, bias_fit$effects$drug)

  bias_estimate <- bias_fit$effects$estimate[bias_pos]
  bias_covariance <- bias_fit$covariance[common, common, drop = FALSE]

  shrinkage <- NULL
  if (shrink != "none") {
    se2 <- diag(bias_covariance)
    # bias_fit's own pinned reference drug has estimate/std_error fixed at
    # exactly 0 by construction (fit_surface()'s normalization), not measured
    # -- it carries no information about between-drug bias variance and would
    # divide by zero in `w_fixed` below, so it is excluded from estimating
    # tau2 and given shrink weight 1 (trusted exactly, never pulled).
    estimable <- se2 > 0
    if (sum(estimable) < 2) {
      stop("`shrink` requires at least 2 drugs, other than `bias_fit`'s own ",
           "pinned reference, common to `fit` and `bias_fit`; tau^2 (the ",
           "between-drug bias variance) is not estimable from fewer.",
           call. = FALSE)
    }

    w_fixed <- 1 / se2[estimable]
    mu_fixed <- sum(w_fixed * bias_estimate[estimable]) / sum(w_fixed)

    # DerSimonian-Laird method-of-moments tau^2, across drugs' bias estimates
    # rather than across one drug's trials (contrast pool_meta()).
    Q <- sum(w_fixed * (bias_estimate[estimable] - mu_fixed)^2)
    df <- sum(estimable) - 1
    C <- sum(w_fixed) - sum(w_fixed^2) / sum(w_fixed)
    tau2 <- if (C > 0) max(0, (Q - df) / C) else 0

    target <- switch(shrink, zero = 0, mean = mu_fixed)
    weight <- ifelse(estimable, tau2 / (tau2 + se2), 1)

    shrinkage <- list(method = shrink, tau2 = tau2, target = target,
                       weight = stats::setNames(weight, common))

    bias_estimate <- target + weight * (bias_estimate - target)
    bias_covariance <- bias_covariance * outer(sqrt(weight), sqrt(weight))
  }

  estimate <- fit$effects$estimate[fit_pos] - bias_estimate
  covariance <- fit$covariance[common, common, drop = FALSE] + bias_covariance
  dimnames(covariance) <- list(common, common)

  z <- stats::qnorm(0.975)
  se <- sqrt(diag(covariance))
  debiased <- data.frame(
    drug = common,
    estimate = estimate,
    std_error = se,
    scale = fit$effects$scale[fit_pos],
    reference = fit$effects$reference[fit_pos],
    engine = fit$engine,
    row.names = NULL
  )
  debiased$lower <- debiased$estimate - z * debiased$std_error
  debiased$upper <- debiased$estimate + z * debiased$std_error
  debiased <- debiased[, c("drug", "estimate", "std_error", "lower", "upper",
                           "scale", "reference", "engine")]

  out <- new_directeffect_fit(
    effects = debiased,
    covariance = covariance,
    heterogeneity = fit$heterogeneity,
    engine = fit$engine,
    engine_fit = fit$engine_fit,
    network = fit$network
  )
  out$bias_fit <- bias_fit
  out$shrinkage <- shrinkage
  out
}
