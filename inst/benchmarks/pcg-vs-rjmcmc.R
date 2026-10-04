# Reproduce the PCG/RAM versus standard RJMCMC comparison.
# Usage: Rscript pcg-vs-rjmcmc.R output-directory [iterations] [burn] [chains]
# Requires the installed development version of poistree and package posterior.
# Both methods use identical data, priors, tree controls, and prediction points.
# Elapsed time includes warm-up, retained-state storage, and predictions.
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1L]] else "pcg-benchmark"
iterations <- if (length(args) >= 2L) as.integer(args[[2L]]) else 12000L
burn <- if (length(args) >= 3L) as.integer(args[[3L]]) else 4000L
chains <- if (length(args) >= 4L) as.integer(args[[4L]]) else 4L
stopifnot(iterations > burn, burn >= 0L, chains >= 2L)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stopifnot(requireNamespace("poistree", quietly = TRUE),
          requireNamespace("posterior", quietly = TRUE))

make_case <- function(name, d, n, seed, a_gate, b_gate) {
  set.seed(seed)
  group <- sample.int(3L, n, replace = TRUE, prob = c(.425, .425, .15))
  x <- matrix(runif(n * d), n, d)
  for (j in seq_len(d)) {
    left <- if (j %% 2L) 1L else 2L
    x[group == left, j] <- rbeta(sum(group == left), 3, 10)
    x[group == 3L - left, j] <- rbeta(sum(group == 3L - left), 10, 3)
  }
  colnames(x) <- paste0("x", seq_len(d))
  mon <- rbind(rep(.23, d), rep(.5, d), rep(.77, d))
  if (d > 1L) {
    mon[1L, seq.int(2L, d, by = 2L)] <- .77
    mon[3L, seq.int(2L, d, by = 2L)] <- .23
  }
  list(name = name, x = x, region = cbind(rep(0, d), rep(1, d)),
       mon = mon, a_gate = a_gate, b_gate = b_gate,
       gate = a_gate / b_gate, data_seed = seed)
}
cases <- list(
  make_case("1d_default_gate_prior", 1L, 160L, 7301L, 36, 3),
  make_case("2d_broad_gate_prior", 2L, 200L, 7302L, 4, .4),
  make_case("4d_default_gate_prior", 4L, 240L, 7303L, 36, 3)
)
# A smoother two-dimensional example supplements the separated-mode examples.
set.seed(7304L)
smooth <- make_case("2d_smooth_broad_prior", 2L, 120L, 7304L, 4, .4)
set.seed(7304L)
smooth$x <- matrix(rbeta(240L, 2, 3), ncol = 2L,
                   dimnames = list(NULL, c("x1", "x2")))
cases[[4L]] <- smooth
cases[[5L]] <- make_case("10d_default_gate_prior", 10L, 400L, 7305L, 36, 3)
for (i in seq_along(cases)) cases[[i]]$run_id <- i
selected <- Sys.getenv("PCG_BENCH_CASES", "")
if (nzchar(selected)) cases <- Filter(function(x) x$name %in%
  strsplit(selected, ",", fixed = TRUE)[[1L]], cases)
stopifnot(length(cases) > 0L)
saveRDS(cases, file.path(out, "data.rds"))

draw_matrix <- function(fit) {
  gate <- fit$posterior$state$gate
  colnames(gate) <- paste0("gamma", seq_len(ncol(gate)))
  lambda <- t(fit$prediction$draws)
  colnames(lambda) <- paste0("intensity", seq_len(ncol(lambda)))
  cbind(gate, leaves = fit$diagnostics$leaf_count_trace,
        integrated_intensity = fit$posterior$integrated_intensity_draws,
        lambda)
}
as_array <- function(mat, chain) {
  ids <- unique(chain)
  n <- sum(chain == ids[[1L]])
  stopifnot(all(table(chain) == n))
  arr <- array(NA_real_, c(n, length(ids), ncol(mat)),
               dimnames = list(iteration = seq_len(n), chain = ids,
                               variable = colnames(mat)))
  for (j in seq_along(ids)) arr[, j, ] <- mat[chain == ids[[j]], , drop = FALSE]
  posterior::as_draws_array(arr)
}

all_stats <- list()
runs <- list()
fits <- list()
for (case_index in seq_along(cases)) {
  dat <- cases[[case_index]]
  # Alternate method order to reduce a consistent timing-order effect.
  methods <- if (dat$run_id %% 2L) c("rjmcmc", "pcg") else c("pcg", "rjmcmc")
  for (method in methods) {
    key <- paste(dat$name, method, sep = "__")
    cat("Running", key, "\n")
    controls <- list(x = dat$x, region = dat$region, gating = "soft",
      sampler = method, gate_family = "logistic", gate_structure = "dimension",
      predict_at = dat$mon, a = .5, b = .5 / nrow(dat$x),
      gate = dat$gate, a_gate = dat$a_gate, b_gate = dat$b_gate,
      sd_gate = .07, update_gate = TRUE, max_depth = 4L, min_leaf_n = 1L,
      chains = chains, iter = iterations, burn = burn, thin = 1L,
      tree_moves = 3L, change_moves = 8L, cut_candidates = 50L,
      cut_proposal = "quantile",
      seed = 9100L + dat$run_id + if (method == "pcg") 1000L else 0L,
      verbose = FALSE)
    gc()
    elapsed <- system.time(fit <- do.call(poistree::ppt_fit, controls))[["elapsed"]]
    mat <- draw_matrix(fit)
    chain <- fit$posterior$state$chain
    draws <- as_array(mat, chain)
    stats <- as.data.frame(posterior::summarise_draws(draws, "mean", "sd",
      "mcse_mean", "ess_bulk", "ess_tail", "rhat"))
    stats$ess_bulk_per_second <- stats$ess_bulk / elapsed
    stats$case <- dat$name
    stats$sampler <- method
    stats$elapsed_seconds <- elapsed
    stats$lag1_mean <- vapply(seq_len(ncol(mat)), function(k) {
      mean(vapply(unique(chain), function(j) {
        z <- mat[chain == j, k]
        if (sd(z) == 0) return(NA_real_)
        cor(z[-length(z)], z[-1L])
      }, numeric(1)), na.rm = TRUE)
    }, numeric(1))
    all_stats[[key]] <- stats
    runs[[key]] <- data.frame(case = dat$name, sampler = method,
      elapsed_seconds = elapsed, chains = chains, iterations = iterations,
      burn = burn, retained_per_chain = iterations - burn,
      ram_failures = if (method == "pcg") sum(fit$diagnostics$ram$failures) else NA_real_)
    fits[[key]] <- fit
    saveRDS(list(fit = fit, elapsed = elapsed, controls = controls,
                 statistics = stats), file.path(out, paste0(key, ".rds")))
    write.csv(do.call(rbind, all_stats), file.path(out, "diagnostics.csv"), row.names = FALSE)
    write.csv(do.call(rbind, runs), file.path(out, "timing.csv"), row.names = FALSE)
    cat("Elapsed", elapsed, "seconds; max Rhat", max(stats$rhat, na.rm = TRUE), "\n")
  }
}
stats <- do.call(rbind, all_stats)
reference <- stats[stats$sampler == "rjmcmc", ]
candidate <- stats[stats$sampler == "pcg", ]
comparison <- merge(reference, candidate, by = c("case", "variable"),
                    suffixes = c("_rjmcmc", "_pcg"))
comparison$ess_per_second_ratio <- comparison$ess_bulk_per_second_pcg /
  comparison$ess_bulk_per_second_rjmcmc
comparison$mean_difference <- comparison$mean_pcg - comparison$mean_rjmcmc
comparison$combined_mcse <- sqrt(comparison$mcse_mean_pcg^2 +
                                  comparison$mcse_mean_rjmcmc^2)
comparison$mean_difference_in_mcse <- comparison$mean_difference / comparison$combined_mcse
write.csv(comparison, file.path(out, "comparison.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))

colors <- c("#176B9A", "#D66B28", "#459A71", "#925AB0")
for (dat in cases) {
  keys <- paste(dat$name, c("rjmcmc", "pcg"), sep = "__")
  png(file.path(out, paste0(dat$name, "-traces.png")), width = 1800,
      height = 1800, res = 170, type = "cairo")
  par(mfrow = c(3, 2), mar = c(4.6, 4.1, 2.6, 1))
  for (variable in c("gamma1", "leaves", "intensity1")) {
    range_y <- range(unlist(lapply(fits[keys], function(f) draw_matrix(f)[, variable])))
    for (key in keys) {
      fit <- fits[[key]]
      mat <- draw_matrix(fit)
      chain <- fit$posterior$state$chain
      seqs <- lapply(unique(chain), function(j) mat[chain == j, variable])
      plot(seqs[[1L]], type = "l", col = colors[[1L]], ylim = range_y,
           xlab = "Retained iteration", ylab = variable,
           main = paste(fit$model$sampler, variable))
      for (j in seq_along(seqs)[-1L]) lines(seqs[[j]], col = colors[j])
    }
  }
  dev.off()
}
cat("Results saved to", normalizePath(out), "\n")
