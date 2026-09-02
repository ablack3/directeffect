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
#'
#' @return A `directeffect_fit` restricted to the drugs `fit` and `bias_fit`
#'   have in common (a warning names any dropped), with `estimate` equal to
#'   `fit`'s estimate minus `bias_fit`'s estimate per drug and `covariance`
#'   equal to the (assumed-independent) sum. `$bias_fit` on the result carries
#'   the bias surface used, for provenance.
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
#' }
#' @export
debias_surface <- function(fit, bias_fit, independent = TRUE) {
  assert_directeffect_fit(fit)
  assert_directeffect_fit(bias_fit)
  assert_fit_covariance(fit, "debias_surface()")
  assert_fit_covariance(bias_fit, "debias_surface()")
  if (!isTRUE(independent)) {
    stop("`debias_surface()` currently only supports `independent = TRUE`; ",
         "there is no general way to estimate the covariance between two ",
         "independently-fit surfaces.", call. = FALSE)
  }

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

  estimate <- fit$effects$estimate[fit_pos] - bias_fit$effects$estimate[bias_pos]
  covariance <- fit$covariance[common, common, drop = FALSE] +
    bias_fit$covariance[common, common, drop = FALSE]
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
  out
}
