#' Pool comparative estimates by bootstrapped mean across trials
#'
#' Averages a drug's comparative estimates directly, with a non-parametric
#' bootstrap over the estimates for the interval. This is the pooling used in
#' high-throughput target-trial emulation (Zang et al., *Nature Communications*
#' 2023;14:8180): per drug, take the sample mean of the adjusted hazard ratios
#' over the trials that balanced, bootstrap the trials 1,000 times for a 95%
#' interval, and test whether that mean is below 1.
#'
#' It is offered here as a deliberate contrast to [fit_surface()], not as an
#' alternative implementation of it. The two answer different questions:
#'
#' * `pool_bootstrap()` asks **"how did this drug fare against the comparators it
#'   happened to be tested against?"** Each estimate enters as-is, so the answer
#'   inherits whatever the comparators' own effects were. A drug compared mostly
#'   against harmful drugs looks protective, and nothing in the method can tell
#'   that apart from the drug being protective.
#' * [fit_surface()] asks **"where does this drug sit relative to every other
#'   drug?"**, solving the whole network at once so a comparator's own position
#'   is estimated rather than assumed to be neutral.
#'
#' Where the two disagree, the disagreement is informative: it measures how far
#' a drug's comparator set sits from the network average; [compare_pooling()]
#' reports that quantity per drug.
#'
#' @param de A `directeffect_network` from [direct_effect_network()], or a plain
#'   data frame of comparisons with `target`, `comparator` and `estimate`.
#' @param balance Optional data frame with `target`, `comparator` and
#'   `max_abs_smd`, used to keep only balanced trials. The source paper counts a
#'   trial as balanced when at most 2% of covariates exceed |SMD| 0.1; supply
#'   whatever balance summary the upstream pipeline produced and set
#'   `max_smd` to the threshold.
#' @param max_smd Balance threshold applied to `balance$max_abs_smd`. Trials
#'   above it are dropped. Ignored when `balance` is `NULL`.
#' @param min_balanced_frac A drug is reported only if at least this fraction of
#'   its trials survive the balance filter. The source paper uses 0.10. Ignored
#'   when `balance` is `NULL`.
#' @param R Bootstrap resamples. Default 1000, as in the source paper.
#' @param scale `"hr"` averages on the hazard-ratio (linear) scale, which is
#'   what the source paper does. `"log"` averages log hazard ratios and
#'   exponentiates, which is the scale the estimates are actually symmetric on
#'   and is usually the better choice; it is not what the paper did.
#' @param seed Optional integer seed. The caller's RNG state is restored on
#'   exit, so pooling never disturbs the session's reproducibility.
#'
#' @return A data frame with one row per drug that appears as a target:
#'   `drug`, `n_trials`, `n_balanced`, `frac_balanced`, `pooled_hr`,
#'   `ci_lower`, `ci_upper`, `p_value` (bootstrap, one-sided against
#'   HR < 1, as in the source paper), `scale`, and `reported` (whether it
#'   cleared `min_balanced_frac`).
#'
#' @references
#' Zang C, Zhang H, Xu J, et al. High-throughput target trial emulation for
#' Alzheimer's disease drug repurposing with real-world data.
#' *Nature Communications* 2023;14:8180.
#'
#' @seealso [fit_surface()] for the network alternative, and
#'   [compare_pooling()] to run both on one network and quantify the gap.
#'
#' @examples
#' de <- direct_effect_network(example_comparisons, effect_measure = "HR")
#' pool_bootstrap(de, R = 200, seed = 1)
#' @export
pool_bootstrap <- function(de,
                           balance = NULL,
                           max_smd = 0.1,
                           min_balanced_frac = 0.1,
                           R = 1000,
                           scale = c("hr", "log"),
                           seed = NULL) {
  scale <- match.arg(scale)
  comparisons <- if (inherits(de, "directeffect_network")) de$comparisons else de
  required <- c("target", "comparator", "estimate")
  missing <- setdiff(required, names(comparisons))
  if (length(missing)) {
    stop("`de` must supply columns ", paste0("`", missing, "`", collapse = ", "),
         ".", call. = FALSE)
  }

  comparisons <- comparisons[is.finite(comparisons$estimate), , drop = FALSE]

  # Balance filter. Kept separate from the estimates so a caller can pool
  # unfiltered -- the filter is the source paper's, not an intrinsic part of
  # averaging, and its effect is worth being able to isolate.
  if (!is.null(balance)) {
    need <- c("target", "comparator", "max_abs_smd")
    if (!all(need %in% names(balance))) {
      stop("`balance` must have columns `target`, `comparator`, `max_abs_smd`.",
           call. = FALSE)
    }
    key <- paste(comparisons$target, comparisons$comparator, sep = "\r")
    bkey <- paste(balance$target, balance$comparator, sep = "\r")
    smd <- balance$max_abs_smd[match(key, bkey)]
    # An estimate with no balance record cannot be shown to be balanced, so it
    # is dropped rather than silently admitted.
    comparisons$.balanced <- !is.na(smd) & smd <= max_smd
  } else {
    comparisons$.balanced <- TRUE
  }

  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }

  drugs <- sort(unique(comparisons$target))
  rows <- lapply(drugs, function(d) {
    trials <- comparisons[comparisons$target == d, , drop = FALSE]
    kept <- trials[trials$.balanced, , drop = FALSE]
    n_trials <- nrow(trials)
    n_bal <- nrow(kept)
    frac <- if (n_trials > 0) n_bal / n_trials else NA_real_
    reported <- is.na(frac) || frac >= min_balanced_frac

    if (n_bal == 0) {
      return(data.frame(
        drug = d, n_trials = n_trials, n_balanced = 0L, frac_balanced = frac,
        pooled_hr = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
        p_value = NA_real_, scale = scale, reported = FALSE,
        stringsAsFactors = FALSE
      ))
    }

    # The paper averages hazard ratios, not log hazard ratios. Reproduced
    # faithfully under scale = "hr"; scale = "log" is the better-behaved
    # alternative and is offered but not the default.
    vals <- if (scale == "hr") exp(kept$estimate) else kept$estimate
    point <- mean(vals)

    # Bootstrap resamples TRIALS, so the interval reflects spread between
    # comparisons and not the within-comparison standard error. That is the
    # source paper's choice and its main statistical cost: a 500-patient
    # emulation counts as much as a 500,000-patient one.
    if (n_bal == 1) {
      boot <- rep(point, R)
    } else {
      idx <- matrix(
        sample.int(n_bal, n_bal * R, replace = TRUE), nrow = R, ncol = n_bal
      )
      boot <- rowMeans(matrix(vals[idx], nrow = R))
    }

    ci <- stats::quantile(boot, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
    null_value <- if (scale == "hr") 1 else 0
    # One-sided bootstrap p, matching the paper's test of "is the mean < 1".
    p <- mean(boot >= null_value)

    if (scale == "log") {
      point <- exp(point); ci <- exp(ci)
    }

    data.frame(
      drug = d, n_trials = n_trials, n_balanced = n_bal, frac_balanced = frac,
      pooled_hr = point, ci_lower = ci[1], ci_upper = ci[2], p_value = p,
      scale = scale, reported = reported, stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


#' Compare bootstrapped pooling with the fitted surface
#'
#' Runs [pool_bootstrap()] and lines its per-drug answer up against a
#' [fit_surface()] fit of the same network, on the same reference scale.
#'
#' The two are put on a common footing by re-expressing the surface relative to
#' each drug's own comparators: a drug's surface position minus the
#' evidence-weighted mean position of the drugs it was compared against. That
#' difference is what pooling implicitly assumes is zero -- it treats every
#' comparator as if it sat at the network average. `comparator_offset` reports
#' how wrong that assumption is for each drug, and is the quantity that explains
#' where the two methods part company.
#'
#' @param de A `directeffect_network`.
#' @param fit A `directeffect_fit` from [fit_surface()] on the same network. If
#'   `NULL`, one is fitted with the `netmeta` engine.
#' @param ... Passed to [pool_bootstrap()].
#'
#' @return A data frame with one row per target drug: the pooled estimate, the
#'   surface estimate, `comparator_offset`, and the difference in log HR.
#' @export
compare_pooling <- function(de, fit = NULL, ...) {
  if (!inherits(de, "directeffect_network")) {
    stop("`de` must be a `directeffect_network`.", call. = FALSE)
  }
  if (is.null(fit)) fit <- fit_surface(de, engine = "netmeta")

  pooled <- pool_bootstrap(de, ...)
  eff <- fit$effects
  pos <- stats::setNames(eff$estimate, eff$drug)

  comparisons <- de$comparisons
  offset <- vapply(pooled$drug, function(d) {
    cs <- comparisons$comparator[comparisons$target == d]
    cs <- cs[!is.na(cs)]
    if (!length(cs)) return(NA_real_)
    mean(pos[cs], na.rm = TRUE)
  }, numeric(1))

  out <- data.frame(
    drug = pooled$drug,
    n_trials = pooled$n_trials,
    n_balanced = pooled$n_balanced,
    pooled_hr = pooled$pooled_hr,
    pooled_lower = pooled$ci_lower,
    pooled_upper = pooled$ci_upper,
    surface_log_hr = unname(pos[pooled$drug]),
    surface_hr = exp(unname(pos[pooled$drug])),
    comparator_offset = unname(offset),
    stringsAsFactors = FALSE
  )
  out$pooled_log_hr <- log(out$pooled_hr)
  # If pooling's assumption held, pooled log HR would equal the drug's position
  # minus its comparators' mean position. The residual is what is left over.
  out$implied_by_surface <- out$surface_log_hr - out$comparator_offset
  out$log_hr_diff <- out$pooled_log_hr - out$implied_by_surface
  rownames(out) <- NULL
  out
}
