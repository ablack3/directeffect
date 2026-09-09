test_that("debias_surface subtracts a known bias and sums the variance", {
  skip_if_not_installed("netmeta")

  comparisons <- data.frame(
    study_id   = c("S1", "S2", "S3"),
    target     = c("A", "A", "B"),
    comparator = c("B", "C", "C"),
    estimate   = c(0.5, 0.9, 0.4),
    std_error  = c(0.05, 0.05, 0.05)
  )
  # Every edge biased by +0.2 relative to `comparisons` above, except S1 which
  # carries no bias -- a null A-B edge biased at 0 confirms debiasing doesn't
  # require every edge to share the same shift.
  bias <- data.frame(
    study_id   = c("S1", "S2", "S3"),
    target     = c("A", "A", "B"),
    comparator = c("B", "C", "C"),
    estimate   = c(0.0, 0.2, 0.2),
    std_error  = c(0.05, 0.05, 0.05)
  )

  fit <- fit_surface(direct_effect_network(comparisons), engine = "netmeta")
  bias_fit <- fit_surface(direct_effect_network(bias), engine = "netmeta")
  debiased <- debias_surface(fit, bias_fit)

  expect_s3_class(debiased, "directeffect_fit")
  raw <- stats::setNames(fit$effects$estimate, fit$effects$drug)
  bias_est <- stats::setNames(bias_fit$effects$estimate, bias_fit$effects$drug)
  got <- stats::setNames(debiased$effects$estimate, debiased$effects$drug)

  for (drug in names(got)) {
    expect_equal(got[[drug]], raw[[drug]] - bias_est[[drug]], tolerance = 1e-8)
  }

  # Independent-sum variance: strictly larger than either input's on its own.
  se_fit <- stats::setNames(fit$effects$std_error, fit$effects$drug)
  se_bias <- stats::setNames(bias_fit$effects$std_error, bias_fit$effects$drug)
  se_out <- stats::setNames(debiased$effects$std_error, debiased$effects$drug)
  for (drug in names(se_out)) {
    expect_equal(se_out[[drug]]^2, se_fit[[drug]]^2 + se_bias[[drug]]^2,
                 tolerance = 1e-8)
  }
})

test_that("debias_surface drops drugs missing from bias_fit, with a warning", {
  skip_if_not_installed("netmeta")

  comparisons <- data.frame(
    study_id   = c("S1", "S2", "S3"),
    target     = c("A", "A", "B"),
    comparator = c("B", "C", "C"),
    estimate   = c(0.5, 0.9, 0.4),
    std_error  = c(0.05, 0.05, 0.05)
  )
  bias <- data.frame(
    study_id   = "S1",
    target     = "A",
    comparator = "B",
    estimate   = 0.0,
    std_error  = 0.05
  )

  fit <- fit_surface(direct_effect_network(comparisons), engine = "netmeta")
  bias_fit <- fit_surface(direct_effect_network(bias), engine = "netmeta")

  expect_warning(
    debiased <- debias_surface(fit, bias_fit),
    "have no bias estimate"
  )
  expect_setequal(debiased$effects$drug, c("A", "B"))
})

test_that("debias_surface(shrink) matches an independently-computed DL oracle", {
  skip_if_not_installed("netmeta")

  # A is the pinned reference in both fits (estimate/std_error exactly 0,
  # not measured). D's three edges into A/B/C share a +0.3 bias; E hangs off
  # D with no bias-network edge of its own (dropped, as tested above).
  comparisons <- data.frame(
    study_id   = c("S1", "S2", "S3", "S4", "S5", "S6", "S7"),
    target     = c("A",  "A",  "B",  "A",  "B",  "C",  "D"),
    comparator = c("B",  "C",  "C",  "D",  "D",  "D",  "E"),
    estimate   = c(0.0, -0.2, -0.2,  0.2,  0.2,  0.4,  0.2),
    std_error  = rep(0.05, 7)
  )
  bias <- data.frame(
    study_id   = c("S1", "S2", "S3", "S4", "S5", "S6"),
    target     = c("A",  "A",  "B",  "A",  "B",  "C"),
    comparator = c("B",  "C",  "C",  "D",  "D",  "D"),
    estimate   = c(0.0,  0.0,  0.0,  0.3,  0.3,  0.3),
    std_error  = rep(0.05, 6)
  )

  fit <- fit_surface(direct_effect_network(comparisons), engine = "netmeta")
  bias_fit <- fit_surface(direct_effect_network(bias), engine = "netmeta")

  # Independent oracle: DL tau2 across B, C, D only (A is the pinned
  # reference, se = 0, and must be excluded or this divides by zero).
  be <- stats::setNames(bias_fit$effects$estimate, bias_fit$effects$drug)
  se <- stats::setNames(bias_fit$effects$std_error, bias_fit$effects$drug)
  keep <- c("B", "C", "D")
  w <- 1 / se[keep]^2
  mu <- sum(w * be[keep]) / sum(w)
  Q <- sum(w * (be[keep] - mu)^2)
  C_stat <- sum(w) - sum(w^2) / sum(w)
  tau2_oracle <- max(0, (Q - length(keep) + 1) / C_stat)
  weight_oracle <- tau2_oracle / (tau2_oracle + se[keep]^2)

  expect_warning(zero <- debias_surface(fit, bias_fit, shrink = "zero"),
                 "have no bias estimate")
  expect_equal(zero$shrinkage$method, "zero")
  expect_equal(zero$shrinkage$tau2, tau2_oracle, tolerance = 1e-8)
  expect_equal(zero$shrinkage$target, 0)
  expect_equal(unname(zero$shrinkage$weight["A"]), 1)
  expect_equal(zero$shrinkage$weight[keep], weight_oracle, tolerance = 1e-8)

  shrunk_bias_zero <- 0 + zero$shrinkage$weight * (be[names(zero$shrinkage$weight)] - 0)
  raw <- stats::setNames(fit$effects$estimate, fit$effects$drug)
  got <- stats::setNames(zero$effects$estimate, zero$effects$drug)
  for (drug in names(shrunk_bias_zero)) {
    expect_equal(got[[drug]], raw[[drug]] - shrunk_bias_zero[[drug]], tolerance = 1e-8)
  }

  expect_warning(meanfit <- debias_surface(fit, bias_fit, shrink = "mean"),
                 "have no bias estimate")
  expect_equal(meanfit$shrinkage$target, mu, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(meanfit$effects$estimate, zero$effects$estimate)))

  # Shrinkage (weight < 1 on B, C, D above) pulls the raw bias estimate
  # partway back toward 0, so it subtracts LESS than the full unshrunk
  # correction -- the shrunk debiased D must sit strictly between D's raw
  # (still-biased) fit_surface() estimate and the fully unshrunk correction.
  raw_bias <- suppressWarnings(debias_surface(fit, bias_fit))
  d <- function(x) x$effects$estimate[x$effects$drug == "D"]
  expect_true(d(fit) < d(zero))
  expect_true(d(zero) < d(raw_bias))
})

test_that("debias_surface(shrink) shrinks variance and errors with too few estimable drugs", {
  skip_if_not_installed("netmeta")

  comparisons <- data.frame(
    study_id   = c("S1", "S2", "S3"),
    target     = c("A", "A", "B"),
    comparator = c("B", "C", "C"),
    estimate   = c(0.5, 0.9, 0.4),
    std_error  = c(0.05, 0.05, 0.05)
  )
  bias <- data.frame(
    study_id   = c("S1", "S2", "S3"),
    target     = c("A", "A", "B"),
    comparator = c("B", "C", "C"),
    estimate   = c(0.0, 0.2, 0.2),
    std_error  = c(0.05, 0.05, 0.05)
  )

  fit <- fit_surface(direct_effect_network(comparisons), engine = "netmeta")
  bias_fit <- fit_surface(direct_effect_network(bias), engine = "netmeta")

  shrunk <- debias_surface(fit, bias_fit, shrink = "zero")
  unshrunk <- debias_surface(fit, bias_fit)
  # Shrunk variance must not exceed the unshrunk (independent-sum) variance:
  # borrowing strength across drugs can only add information.
  expect_true(all(shrunk$effects$std_error <= unshrunk$effects$std_error + 1e-8))

  # Only two drugs total, A pinned as bias_fit's reference (se = 0) -- just
  # one estimable (non-reference) drug, not enough to fit tau2 on.
  comparisons2 <- data.frame(
    study_id = "S1", target = "A", comparator = "B",
    estimate = 0.5, std_error = 0.05
  )
  bias2 <- data.frame(
    study_id = "S1", target = "A", comparator = "B",
    estimate = 0.1, std_error = 0.05
  )
  fit2 <- fit_surface(direct_effect_network(comparisons2), engine = "netmeta")
  bias_fit2 <- fit_surface(direct_effect_network(bias2), engine = "netmeta")
  expect_error(
    debias_surface(fit2, bias_fit2, shrink = "mean"),
    "at least 2"
  )
})

test_that("debias_surface errors on independent = FALSE", {
  skip_if_not_installed("netmeta")

  comparisons <- data.frame(
    study_id   = c("S1"),
    target     = "A",
    comparator = "B",
    estimate   = 0.5,
    std_error  = 0.05
  )
  fit <- fit_surface(direct_effect_network(comparisons), engine = "netmeta")

  expect_error(
    debias_surface(fit, fit, independent = FALSE),
    "only supports"
  )
})
