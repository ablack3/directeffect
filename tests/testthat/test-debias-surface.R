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
