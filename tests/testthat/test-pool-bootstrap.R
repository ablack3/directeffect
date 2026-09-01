# pool_bootstrap() is checked against arithmetic, not against another
# implementation: on a fixed set of trials the point estimate is a plain mean,
# so it can be written down by hand.

pool_comparisons <- function() {
  data.frame(
    study_id   = paste0("S", 1:6),
    target     = c("A", "A", "A", "B", "B", "C"),
    comparator = c("X", "Y", "Z", "X", "Y", "X"),
    estimate   = log(c(0.5, 0.8, 1.0, 1.2, 1.5, 0.9)),
    std_error  = c(0.10, 0.20, 0.40, 0.10, 0.20, 0.15)
  )
}

pool_balance <- function() {
  data.frame(
    target      = c("A", "A", "A", "B", "B", "C"),
    comparator  = c("X", "Y", "Z", "X", "Y", "X"),
    max_abs_smd = c(0.02, 0.05, 0.50, 0.03, 0.60, 0.01)
  )
}

test_that("the point estimate is the arithmetic mean on the requested scale", {
  de <- direct_effect_network(pool_comparisons(), effect_measure = "HR")
  hr  <- pool_bootstrap(de, R = 50, seed = 1, scale = "hr")
  lg  <- pool_bootstrap(de, R = 50, seed = 1, scale = "log")

  # A's three trials: HR 0.5, 0.8, 1.0.
  expect_equal(hr$pooled_hr[hr$drug == "A"], mean(c(0.5, 0.8, 1.0)))
  expect_equal(lg$pooled_hr[lg$drug == "A"], exp(mean(log(c(0.5, 0.8, 1.0)))))
  # Only drugs appearing as a target are reported.
  expect_identical(hr$drug, c("A", "B", "C"))
})

test_that("averaging ratios is >= averaging logs, which is Jensen's inequality", {
  de <- direct_effect_network(pool_comparisons(), effect_measure = "HR")
  hr <- pool_bootstrap(de, R = 50, seed = 1, scale = "hr")
  lg <- pool_bootstrap(de, R = 50, seed = 1, scale = "log")
  # Equality only when a drug has a single trial (C); strict otherwise.
  expect_true(all(hr$pooled_hr >= lg$pooled_hr - 1e-12))
  expect_gt(hr$pooled_hr[hr$drug == "A"], lg$pooled_hr[lg$drug == "A"])
  expect_equal(hr$pooled_hr[hr$drug == "C"], lg$pooled_hr[lg$drug == "C"])
})

test_that("the balance filter drops trials and gates reporting", {
  de <- direct_effect_network(pool_comparisons(), effect_measure = "HR")
  out <- pool_bootstrap(de, balance = pool_balance(), max_smd = 0.1,
                        min_balanced_frac = 0.5, R = 50, seed = 1, scale = "hr")

  a <- out[out$drug == "A", ]
  expect_identical(a$n_trials, 3L)
  expect_identical(a$n_balanced, 2L)          # Z is dropped at SMD 0.50
  expect_equal(a$pooled_hr, mean(c(0.5, 0.8)))
  expect_true(a$reported)                      # 2/3 clears the 0.5 bar

  b <- out[out$drug == "B", ]
  expect_identical(b$n_balanced, 1L)           # Y dropped at 0.60
  expect_equal(b$frac_balanced, 0.5)
  expect_true(b$reported)                      # 0.5 meets the 0.5 bar, inclusive

  # Raising the bar above the surviving fraction is what un-reports it.
  strict <- pool_bootstrap(de, balance = pool_balance(), max_smd = 0.1,
                           min_balanced_frac = 0.60, R = 50, seed = 1)
  expect_false(strict$reported[strict$drug == "B"])
  expect_true(strict$reported[strict$drug == "A"])   # 2/3 = 0.67 clears 0.60
})

test_that("a trial with no balance record is dropped rather than admitted", {
  de <- direct_effect_network(pool_comparisons(), effect_measure = "HR")
  bal <- pool_balance()[1, ]                  # only A-vs-X has a record
  out <- pool_bootstrap(de, balance = bal, max_smd = 0.1,
                        min_balanced_frac = 0, R = 50, seed = 1)
  expect_identical(out$n_balanced[out$drug == "A"], 1L)
  expect_identical(out$n_balanced[out$drug == "B"], 0L)
  expect_true(is.na(out$pooled_hr[out$drug == "B"]))
  expect_false(out$reported[out$drug == "B"])
})

test_that("a single trial gives a degenerate interval rather than an error", {
  de <- direct_effect_network(pool_comparisons(), effect_measure = "HR")
  out <- pool_bootstrap(de, R = 50, seed = 1, scale = "hr")
  c_row <- out[out$drug == "C", ]
  expect_identical(c_row$n_balanced, 1L)
  expect_equal(c_row$ci_lower, c_row$pooled_hr)
  expect_equal(c_row$ci_upper, c_row$pooled_hr)
})

test_that("seeding is reproducible and leaves the caller's RNG untouched", {
  de <- direct_effect_network(pool_comparisons(), effect_measure = "HR")
  set.seed(99)
  before <- .Random.seed
  a <- pool_bootstrap(de, R = 100, seed = 7)
  expect_identical(.Random.seed, before)
  b <- pool_bootstrap(de, R = 100, seed = 7)
  expect_equal(a$ci_lower, b$ci_lower)
})

test_that("it accepts a plain data frame and rejects a malformed one", {
  expect_s3_class(pool_bootstrap(pool_comparisons(), R = 20, seed = 1), "data.frame")
  expect_error(pool_bootstrap(pool_comparisons()[, c("target", "estimate")]),
               "comparator")
  de <- direct_effect_network(pool_comparisons(), effect_measure = "HR")
  expect_error(pool_bootstrap(de, balance = data.frame(x = 1)), "max_abs_smd")
})
