#!/usr/bin/env Rscript
# Which pooling method is right? With simulated data the truth is known, so this
# is answerable rather than arguable.
#
# Four estimators of each drug's latent direct effect:
#   1. bootstrap mean of HRs      (Zang et al. 2023)
#   2. bootstrap mean of log HRs  (same, on the symmetric scale)
#   3. inverse-variance meta-analysis, fixed and random effects
#   4. the network surface        (fit_surface)
#
# Scored on bias, RMSE and interval coverage of the KNOWN theta. The first three
# estimate "drug vs its own comparators"; only the fourth targets "drug vs the
# reference". Whether that distinction matters in practice is exactly what the
# simulation settles.
suppressPackageStartupMessages({ library(devtools) })
devtools::load_all("/Users/adam.black/projects/directeffect", quiet = TRUE)

args_env  <- Sys.getenv
N_REP     <- as.integer(args_env("SIM_REPS", "200"))
N_DRUGS   <- 25
N_COMP    <- 150
ENGINE    <- args_env("SIM_ENGINE", "netmeta")
TAG       <- args_env("SIM_TAG", "netmeta200")

run_one <- function(rep_id) {
  sim <- simulate_direct_effect_network(
    n_drugs = N_DRUGS, n_comparisons = N_COMP, n_anchors = 0,
    heterogeneity = 0.05, seed = 1000 + rep_id, effect_sd = 0.4
  )
  de <- sim$network
  truth <- setNames(sim$truth$theta, sim$truth$drug)

  # The simulator's truth table carries a `placebo` row that is not a node in
  # the network, so the reference has to come from the network's own treatments.
  ref <- de$treatments[1]
  fit <- tryCatch(
    if (ENGINE == "stan") {
      fit_surface(de, engine = "stan", reference = ref,
                  chains = 2, iter = 2000, seed = 7000 + rep_id, refresh = 0)
    } else {
      fit_surface(de, engine = "netmeta", reference = ref)
    },
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  pb_hr  <- pool_bootstrap(de, R = 400, seed = rep_id, scale = "hr")
  pb_log <- pool_bootstrap(de, R = 400, seed = rep_id, scale = "log")
  pm_fix <- pool_meta(de, method = "fixed")
  pm_ran <- pool_meta(de, method = "random")

  # Every estimator is put on the SAME footing: effect relative to the reference
  # drug, which is what `truth` is measured in. For the surface that is what it
  # already estimates. For the pooled estimators it is not -- they estimate the
  # drug minus its own comparators -- and that gap is the thing under test, so
  # it is left uncorrected rather than patched.
  theta_rel <- truth - truth[ref]

  grab <- function(df, est_col, lo, hi) {
    d <- df$drug
    data.frame(drug = d,
               est = log(df[[est_col]]),
               lo = log(df[[lo]]), hi = log(df[[hi]]),
               stringsAsFactors = FALSE)
  }
  out <- list(
    boot_hr   = grab(pb_hr,  "pooled_hr", "ci_lower", "ci_upper"),
    boot_log  = grab(pb_log, "pooled_hr", "ci_lower", "ci_upper"),
    meta_fix  = grab(pm_fix, "pooled_hr", "ci_lower", "ci_upper"),
    meta_ran  = grab(pm_ran, "pooled_hr", "ci_lower", "ci_upper"),
    surface   = data.frame(drug = fit$effects$drug, est = fit$effects$estimate,
                           lo = fit$effects$lower, hi = fit$effects$upper,
                           stringsAsFactors = FALSE)
  )
  pieces <- lapply(names(out), function(m) {
    x <- out[[m]]
    if (is.null(x) || !nrow(x)) return(NULL)
    x$truth  <- unname(theta_rel[x$drug])
    x$method <- rep(m, nrow(x))
    x$rep    <- rep_id
    x[is.finite(x$est) & is.finite(x$truth), , drop = FALSE]
  })
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (!length(pieces)) return(NULL)
  do.call(rbind, pieces)
}

cat(sprintf("[%s] running %d replications, engine=%s, %d drugs, %d comparisons ...\n",
            TAG, N_REP, ENGINE, N_DRUGS, N_COMP)); flush(stdout())
all <- do.call(rbind, lapply(seq_len(N_REP), function(i) {
  if (i %% 25 == 0) { cat("  rep", i, "\n"); flush(stdout()) }
  run_one(i)
}))

all$err <- all$est - all$truth
all$covered <- all$truth >= all$lo & all$truth <= all$hi

summ <- do.call(rbind, lapply(split(all, all$method), function(d) {
  data.frame(method = d$method[1],
             n = nrow(d),
             bias = mean(d$err),
             rmse = sqrt(mean(d$err^2)),
             coverage95 = mean(d$covered),
             median_ci_width = median(d$hi - d$lo),
             stringsAsFactors = FALSE)
}))
summ <- summ[order(summ$rmse), ]
cat("\n=== scored against known truth ===\n")
print(summ, row.names = FALSE, digits = 3)

summ$engine <- ENGINE; summ$reps <- N_REP; summ$tag <- TAG
base <- getwd()
write.csv(summ, file.path(base, paste0("pooling_sim_summary_", TAG, ".csv")), row.names = FALSE)
write.csv(all,  file.path(base, paste0("pooling_sim_raw_", TAG, ".csv")), row.names = FALSE)
cat(sprintf("\nwrote pooling_sim_summary_%s.csv\n", TAG))
