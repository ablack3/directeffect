# pool_meta() is checked against the inverse-variance formulae written out by
# hand, so a change in the arithmetic fails here rather than agreeing with
# itself.

meta_comparisons <- function() {
  data.frame(
    study_id   = paste0("S", 1:5),
    target     = c("A", "A", "A", "B", "C"),
    comparator = c("X", "Y", "Z", "X", "X"),
    estimate   = c(-0.20, -0.10, 0.30, 0.15, -0.05),
    std_error  = c(0.10, 0.20, 0.40, 0.10, 0.15)
  )
}

# Hand oracle for the fixed-effect pool.
iv_fixed <- function(y, se) {
  w <- 1 / se^2
  mu <- sum(w * y) / sum(w)
  list(mu = mu, se = sqrt(1 / sum(w)))
}

test_that("the fixed-effect pool matches the inverse-variance formula", {
  de <- direct_effect_network(meta_comparisons(), effect_measure = "HR")
  out <- pool_meta(de, method = "fixed")
  a <- out[out$drug == "A", ]

  o <- iv_fixed(c(-0.20, -0.10, 0.30), c(0.10, 0.20, 0.40))
  expect_equal(a$log_hr, o$mu)
  expect_equal(a$se_log_hr, o$se)
  expect_equal(a$pooled_hr, exp(o$mu))
  expect_equal(a$ci_lower, exp(o$mu - 1.96 * o$se))
  expect_equal(a$ci_upper, exp(o$mu + 1.96 * o$se))
})

test_that("tau^2 is zero when the spread is no more than the standard errors allow", {
  # A's estimates range over 0.5 log HR but carry SEs up to 0.40, so Q (1.57)
  # falls below its 2 degrees of freedom: no excess heterogeneity to model, and
  # random effects must reduce to the fixed-effect answer rather than invent
  # spread.
  de <- direct_effect_network(meta_comparisons(), effect_measure = "HR")
  fx <- pool_meta(de, method = "fixed")
  rn <- pool_meta(de, method = "random")
  a_fx <- fx[fx$drug == "A", ]; a_rn <- rn[rn$drug == "A", ]
  expect_equal(a_rn$tau2, 0)
  expect_equal(a_rn$i_squared, 0)
  expect_equal(a_rn$se_log_hr, a_fx$se_log_hr)
})

test_that("random effects widens the interval when the trials genuinely disagree", {
  # Same spread, far tighter standard errors: now Q greatly exceeds its df.
  het <- data.frame(
    study_id   = paste0("H", 1:3),
    target     = "A",
    comparator = c("X", "Y", "Z"),
    estimate   = c(-0.5, 0.0, 0.5),
    std_error  = c(0.05, 0.05, 0.05)
  )
  de <- direct_effect_network(het, effect_measure = "HR")
  fx <- pool_meta(de, method = "fixed")
  rn <- pool_meta(de, method = "random")

  expect_gt(rn$tau2, 0)
  expect_gt(rn$se_log_hr, fx$se_log_hr)
  expect_gt(rn$i_squared, 0.9)
  expect_lt(rn$q_p_value, 0.001)
})

test_that("Q, I^2 and tau^2 are undefined for a single trial, not fabricated", {
  de <- direct_effect_network(meta_comparisons(), effect_measure = "HR")
  out <- pool_meta(de, method = "random")
  b <- out[out$drug == "B", ]
  expect_identical(b$n_balanced, 1L)
  expect_true(is.na(b$i_squared))
  expect_true(is.na(b$q_p_value))
  expect_equal(b$tau2, 0)
  # With one trial the pooled estimate is that trial.
  expect_equal(b$log_hr, 0.15)
  expect_equal(b$se_log_hr, 0.10)
})

test_that("the balance filter behaves as it does in pool_bootstrap", {
  de <- direct_effect_network(meta_comparisons(), effect_measure = "HR")
  bal <- data.frame(
    target      = c("A", "A", "A", "B", "C"),
    comparator  = c("X", "Y", "Z", "X", "X"),
    max_abs_smd = c(0.01, 0.02, 0.90, 0.01, 0.01)
  )
  out <- pool_meta(de, balance = bal, max_smd = 0.1, method = "fixed")
  a <- out[out$drug == "A", ]
  expect_identical(a$n_balanced, 2L)
  o <- iv_fixed(c(-0.20, -0.10), c(0.10, 0.20))
  expect_equal(a$log_hr, o$mu)
})

test_that("malformed input is rejected", {
  expect_error(pool_meta(meta_comparisons()[, c("target", "comparator", "estimate")]),
               "std_error")
  de <- direct_effect_network(meta_comparisons(), effect_measure = "HR")
  expect_error(pool_meta(de, balance = data.frame(x = 1)), "max_abs_smd")
})

test_that("compare_pooling lines the surface up against pooling", {
  skip_if_not_installed("netmeta")
  de <- direct_effect_network(meta_comparisons(), effect_measure = "HR")
  cmp <- compare_pooling(de, R = 50, seed = 3)

  expect_true(all(c("drug", "pooled_hr", "surface_log_hr", "comparator_offset",
                    "implied_by_surface", "log_hr_diff") %in% names(cmp)))
  expect_identical(cmp$drug, c("A", "B", "C"))
  # The offset is the mean surface position of that drug's own comparators.
  fit <- fit_surface(de, engine = "netmeta")
  pos <- stats::setNames(fit$effects$estimate, fit$effects$drug)
  expect_equal(cmp$comparator_offset[cmp$drug == "A"],
               mean(pos[c("X", "Y", "Z")]))
  # And the reported difference is defined off it.
  expect_equal(cmp$log_hr_diff, cmp$pooled_log_hr - cmp$implied_by_surface)
  expect_error(compare_pooling(meta_comparisons()), "directeffect_network")
})
