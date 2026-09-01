#' Pool comparative estimates by inverse-variance meta-analysis
#'
#' Combines a drug's comparative estimates by weighting each by its precision,
#' the standard fixed-effect or random-effects (DerSimonian--Laird) meta-analysis.
#'
#' This is the efficient version of what [pool_bootstrap()] does crudely. Both
#' average a drug's estimates as they stand; meta-analysis weights each by
#' `1/se^2` instead of equally, and its interval comes from the pooled variance
#' rather than from resampling trials.
#'
#' **The two share an estimand, and therefore share a bias.** Both answer "how
#' did this drug fare against the comparators it happened to be tested against?"
#' Neither can distinguish a protective drug from one that was simply compared
#' against harmful drugs, because a comparator's own effect never enters the
#' calculation. Weighting fixes efficiency; it does not fix that. Only solving
#' the network ([fit_surface()]) estimates the comparators' positions rather
#' than assuming them to be zero.
#'
#' @param de A `directeffect_network`, or a data frame with `target`,
#'   `comparator`, `estimate` and `std_error`.
#' @param method `"fixed"` for a common-effect inverse-variance pool, or
#'   `"random"` for DerSimonian--Laird, which adds a between-comparison variance
#'   `tau^2` estimated by the method of moments.
#' @param balance,max_smd,min_balanced_frac Optional balance filter, applied
#'   exactly as in [pool_bootstrap()], so the two can be compared on identical
#'   subsets.
#'
#' @return A data frame with one row per drug appearing as a target: `drug`,
#'   `n_trials`, `n_balanced`, `frac_balanced`, `pooled_hr`, `ci_lower`,
#'   `ci_upper`, `log_hr`, `se_log_hr`, `p_value` (two-sided), `tau2`,
#'   `i_squared`, `q_p_value` (Cochran's Q test of homogeneity), `method`,
#'   `reported`.
#'
#' @seealso [pool_bootstrap()], [fit_surface()], [compare_pooling()].
#'
#' @examples
#' de <- direct_effect_network(example_comparisons, effect_measure = "HR")
#' pool_meta(de, method = "random")
#' @export
pool_meta <- function(de,
                      method = c("fixed", "random"),
                      balance = NULL,
                      max_smd = 0.1,
                      min_balanced_frac = 0.1) {
  method <- match.arg(method)
  comparisons <- if (inherits(de, "directeffect_network")) de$comparisons else de
  required <- c("target", "comparator", "estimate", "std_error")
  missing <- setdiff(required, names(comparisons))
  if (length(missing)) {
    stop("`de` must supply columns ", paste0("`", missing, "`", collapse = ", "),
         ".", call. = FALSE)
  }

  comparisons <- comparisons[
    is.finite(comparisons$estimate) & is.finite(comparisons$std_error) &
      comparisons$std_error > 0, , drop = FALSE
  ]

  if (!is.null(balance)) {
    need <- c("target", "comparator", "max_abs_smd")
    if (!all(need %in% names(balance))) {
      stop("`balance` must have columns `target`, `comparator`, `max_abs_smd`.",
           call. = FALSE)
    }
    key <- paste(comparisons$target, comparisons$comparator, sep = "\r")
    bkey <- paste(balance$target, balance$comparator, sep = "\r")
    smd <- balance$max_abs_smd[match(key, bkey)]
    comparisons$.balanced <- !is.na(smd) & smd <= max_smd
  } else {
    comparisons$.balanced <- TRUE
  }

  drugs <- sort(unique(comparisons$target))
  rows <- lapply(drugs, function(d) {
    trials <- comparisons[comparisons$target == d, , drop = FALSE]
    kept <- trials[trials$.balanced, , drop = FALSE]
    n_trials <- nrow(trials)
    k <- nrow(kept)
    frac <- if (n_trials > 0) k / n_trials else NA_real_
    reported <- is.na(frac) || frac >= min_balanced_frac

    empty <- data.frame(
      drug = d, n_trials = n_trials, n_balanced = k, frac_balanced = frac,
      pooled_hr = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
      log_hr = NA_real_, se_log_hr = NA_real_, p_value = NA_real_,
      tau2 = NA_real_, i_squared = NA_real_, q_p_value = NA_real_,
      method = method, reported = FALSE, stringsAsFactors = FALSE
    )
    if (k == 0) return(empty)

    y <- kept$estimate
    v <- kept$std_error^2
    w <- 1 / v
    mu_fixed <- sum(w * y) / sum(w)

    # Cochran's Q and DerSimonian--Laird tau^2, both from the fixed-effect fit.
    Q <- sum(w * (y - mu_fixed)^2)
    df <- k - 1
    tau2 <- 0
    i2 <- NA_real_
    q_p <- NA_real_
    if (df > 0) {
      q_p <- stats::pchisq(Q, df, lower.tail = FALSE)
      i2 <- max(0, (Q - df) / Q)
      C <- sum(w) - sum(w^2) / sum(w)
      if (C > 0) tau2 <- max(0, (Q - df) / C)
    }

    if (method == "random") {
      w <- 1 / (v + tau2)
    }
    mu <- sum(w * y) / sum(w)
    se <- sqrt(1 / sum(w))

    # A single comparison carries no between-comparison information; the
    # random-effects interval collapses to the fixed one, which is honest --
    # tau^2 is not estimable from one estimate.
    z <- mu / se
    p <- 2 * stats::pnorm(-abs(z))

    data.frame(
      drug = d, n_trials = n_trials, n_balanced = k, frac_balanced = frac,
      pooled_hr = exp(mu), ci_lower = exp(mu - 1.96 * se),
      ci_upper = exp(mu + 1.96 * se), log_hr = mu, se_log_hr = se,
      p_value = p, tau2 = tau2, i_squared = i2, q_p_value = q_p,
      method = method, reported = reported, stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
