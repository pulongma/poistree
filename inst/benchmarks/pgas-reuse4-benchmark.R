# Run in fresh R processes against a baseline and optimized package library.
# Usage: Rscript pgas-reuse4-benchmark.R LIBRARY OUTPUT_DIRECTORY
# Extended fitting-time checks; these are not convergence or ESS benchmarks.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
.libPaths(c(normalizePath(args[1]), .libPaths()))
library(poistree)
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
cases <- list(
  exact_2d = list(n = 100L, d = 2L, iter = 300L, particles = 10L,
                  max_depth = 4L, update_gate = TRUE),
  exact_10d = list(n = 150L, d = 10L, iter = 200L, particles = 8L,
                   max_depth = 4L, update_gate = FALSE),
  mixed_10d = list(n = 400L, d = 10L, iter = 150L, particles = 8L,
                   max_depth = 4L, update_gate = TRUE))
results <- list()
for (name in names(cases)) {
  config <- cases[[name]]
  set.seed(241)
  x <- matrix(runif(config$n * config$d), config$n, config$d)
  x[, 1] <- qbeta(x[, 1], 2, 5)
  region <- cbind(rep(0, config$d), rep(1, config$d))
  common <- list(x = x, region = region, predict_at = x[1:12, , drop = FALSE],
                 gating = "soft", sampler = "pgas", a = 0.5, b = 0.5/config$n,
                 gate = 12, a_gate = 36, b_gate = 3, update_gate = config$update_gate,
                 particles = config$particles, iter = config$iter, burn = 10,
                 thin = 2, max_depth = config$max_depth, alpha = 0.9, eta = 1,
                 exact_max = 150, cut_candidates = 50, verbose = FALSE)
  for (repeat_id in 1:5) {
    elapsed <- system.time(fit <- do.call(ppt_fit, c(common, list(seed = 910 + repeat_id))))[["elapsed"]]
    results[[length(results) + 1L]] <- data.frame(
      case = name, repeat_id = repeat_id, n = config$n, d = config$d,
      iter = config$iter, particles = config$particles, cuts = 50,
      seconds = elapsed, mean_leaves = fit$posterior$mean_leaves,
      expanded_nodes = fit$diagnostics$expanded_nodes)
    saveRDS(list(prediction = fit$prediction, posterior = fit$posterior,
                 diagnostics = fit$diagnostics, rng = .Random.seed),
            file.path(args[2], paste0(name, "-", repeat_id, ".rds")))
    write.csv(do.call(rbind, results), file.path(args[2], "timing.csv"), row.names = FALSE)
    cat(name, repeat_id, elapsed, "seconds\n"); flush.console()
  }
}
writeLines(capture.output(sessionInfo()), file.path(args[2], "sessionInfo.txt"))
