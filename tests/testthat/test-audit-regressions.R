audit_root_action <- function(tree) {
  root <- Filter(function(node) {
    !is.null(node) && node$depth == 0L
  }, tree)[[1L]]
  if (isTRUE(root$is_leaf)) "stop" else
    sprintf("%d:%.12g", root$J + 1L, root$L)
}

audit_log_q <- function(n, volume, a, b) {
  lgamma(n + a) - lgamma(a) + a * log(b) -
    (n + a) * log(b + volume)
}

audit_root_law <- function(x, region, cuts, a, b) {
  volume <- prod(region[, 2L] - region[, 1L])
  scores <- c(stop = log(0.5) + audit_log_q(nrow(x), volume, a, b))
  for (j in seq_len(ncol(x))) {
    for (cut in cuts[[j]]) {
      n_left <- sum(x[, j] < cut)
      v_left <- volume * (cut - region[j, 1L]) / diff(region[j, ])
      score <- log(0.5) - log(ncol(x)) - log(length(cuts[[j]])) +
        audit_log_q(n_left, v_left, a, b) +
        audit_log_q(nrow(x) - n_left, volume - v_left, a, b)
      scores[sprintf("%d:%.12g", j, cut)] <- score
    }
  }
  probability <- exp(scores - max(scores))
  probability / sum(probability)
}

test_that("SMC root support and probabilities agree on small boundary designs", {
  unit1 <- matrix(c(0, 1), nrow = 1L)
  unit2 <- matrix(c(0, 1, 0, 1), nrow = 2L, byrow = TRUE)
  cases <- list(

    list(x = matrix(c(0.1, 0.9), ncol = 1L), region = unit1,
         cuts = list(0.5), candidates = 30L),
    list(x = cbind(c(0.2, 0.2, 0.5, 0.5, 0.8, 0.8),
                   c(0.8, 0.8, 0.5, 0.5, 0.2, 0.2)),
         region = unit2, cuts = list(0.5, 0.5), candidates = 30L),

    list(x = rbind(c(0.2, 0.2), c(0.5, 1), c(1, 0.5), c(0.8, 0.8)),
         region = unit2, cuts = list(0.5, 0.5), candidates = 30L),

    list(x = matrix(seq(0.1, 0.9, length.out = 9L), ncol = 1L),
         region = unit1, cuts = list(0.5), candidates = 1L)
  )
  particles <- 2000L
  for (case in cases) {
    exact <- audit_root_law(case$x, case$region, case$cuts, 0.5, 0.1)
    for (engine in c("shared", "dense")) {
      fit <- ppt_fit(
        case$x, case$region, gating = "hard", sampler = "smc",
        engine = engine, max_depth = 1L, particles = particles,
        cut_candidates = case$candidates, a = 0.5, b = 0.1, seed = 270L
      )
      actions <- vapply(fit$posterior$tree_draws,
                        audit_root_action, character(1))
      expect_setequal(unique(actions), names(exact))
      observed <- as.numeric(table(factor(actions, levels = names(exact)))) /
        particles

      error_bound <- 6 * sqrt(exact * (1 - exact) / particles) + 1 / particles
      expect_true(all(abs(observed - exact) < error_bound))
      expect_equal(fit$posterior$particle_weights,
                   rep(1 / particles, particles), tolerance = 1e-12)
    }
  }
})

test_that("hard SMC and PG retain the same predictions on every split tie", {
  coordinates <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
  x <- as.matrix(expand.grid(coordinates, coordinates))

  grid_values <- sort(unique(c(coordinates,
                              as.numeric(outer(coordinates, coordinates, "+") / 2))))
  grid <- as.matrix(expand.grid(grid_values, grid_values))
  region <- matrix(c(0, 1, 0, 1), nrow = 2L, byrow = TRUE)
  for (method in c("shared", "dense", "pgas")) {
    args <- list(
      x = x, region = region, gating = "hard", predict_at = grid, test = x,
      a = 0.5, b = 0.1, max_depth = 2L, seed = 81L
    )
    if (method == "pgas") {
      args <- c(args, list(sampler = "pgas", particles = 12L,
                           chains = 1L, iter = 40L, burn = 10L,
                           verbose = FALSE))
    } else {
      args <- c(args, list(sampler = "smc", engine = method, particles = 100L))
    }
    fit <- do.call(ppt_fit, args)
    expect_true(any(fit$diagnostics$leaf_count_trace > 1L))
    recomputed <- ppt_lambda(fit, at = grid, type = "draws")
    expect_equal(recomputed$draws, fit$prediction$draws, tolerance = 1e-12)
    expect_true(all(recomputed$draws > 0))

    with_new_row <- ppt_predict(
      fit, newdata = rbind(grid, c(0.12345, 0.67891)), type = "mean"
    )
    expect_equal(with_new_row[seq_len(nrow(grid))],
                 ppt_predict(fit, newdata = grid, type = "mean"),
                 tolerance = 1e-12)
    expect_equal(as.numeric(ppt_lppd(fit, test = x)),
                 as.numeric(ppt_lppd(fit)), tolerance = 1e-12)

    for (tree in fit$posterior$tree_draws) {
      leaves <- Filter(function(node) !is.null(node) && isTRUE(node$is_leaf),
                       tree)
      counts <- vapply(leaves, function(leaf) {
        box <- as.matrix(leaf$region)
        inside <- rep(TRUE, nrow(x))
        for (j in seq_len(ncol(x))) {
          inside <- inside & x[, j] >= box[j, 1L] &
            (x[, j] < box[j, 2L] |
               (box[j, 2L] == region[j, 2L] & x[, j] <= region[j, 2L]))
        }
        sum(inside)
      }, numeric(1))
      stored_counts <- vapply(leaves, function(leaf) {
        if (!is.null(leaf$m)) as.numeric(leaf$m) else length(leaf$idx)
      }, numeric(1))
      expect_equal(stored_counts, counts)
      expect_equal(sum(counts), nrow(x))
    }
  }
})

test_that("a scalar hard marginal preserves all unequally weighted draws", {
  region <- matrix(c(0, 1), nrow = 1L)
  rates <- c(2, 7, 11)
  weights <- c(0.1, 0.25, 0.65)
  fit <- structure(list(
    model = list(gating = "hard", scales = "leaf"),
    data = list(x = matrix(0.2, ncol = 1L), dimension = 1L, region = region),
    posterior = list(
      tree_draws = lapply(rates, function(rate) {
        list(list(is_leaf = TRUE, region = region, lambda = rate))
      }),
      particle_weights = weights, state = list(mode = "leafbox")
    )
  ), class = "ppt")
  draws <- ppt_marginal(fit, 1L, grid = 0.3, type = "draws")
  expect_equal(dim(draws), c(1L, 3L))
  expect_equal(as.numeric(draws), rates)
  expect_equal(attr(draws, "weights"), weights)
  summary <- ppt_marginal(fit, 1L, grid = 0.3)
  expect_equal(summary$mean, sum(rates * weights))
  expect_equal(summary$median, 11)
  expect_equal(summary$lower, 2)
  expect_equal(summary$upper, 11)
})

test_that("hard SMC and PG logLik is the mean conditional Poisson likelihood", {
  x <- matrix(c(0, seq(0.05, 0.95, length.out = 30L), 1), ncol = 1L)
  region <- matrix(c(0, 1), nrow = 1L)
  for (method in c("shared", "dense", "pgas")) {
    args <- list(x = x, region = region, gating = "hard", a = 1, b = 0.2,
                 max_depth = 2L, seed = 38L)
    if (method == "pgas") {
      args <- c(args, list(sampler = "pgas", particles = 10L,
                           chains = 1L, iter = 30L, burn = 5L,
                           verbose = FALSE))
    } else {
      args <- c(args, list(sampler = "smc", engine = method, particles = 80L))
    }
    fit <- do.call(ppt_fit, args)
    evaluated <- ppt_lambda(fit, at = x, type = "draws")
    direct <- colSums(log(evaluated$draws)) - poistree:::ppt_integral(fit, "draws")
    expect_equal(as.numeric(poistree:::ppt_logLik(fit)),
                 sum(evaluated$weights * direct), tolerance = 1e-12)
  }
})

test_that("SMC target normalizers include the proper-prior root marginal", {
  x <- matrix(c(0.1, 0.4, 0.9), ncol = 1L)
  region <- matrix(c(0, 2), nrow = 1L)
  a <- 0.7
  b <- 0.4
  expected <- audit_log_q(nrow(x), 2, a, b)
  for (engine in c("shared", "dense")) {
    fit <- ppt_fit(x, region, gating = "hard", sampler = "smc",
                   engine = engine, particles = 10L, max_depth = 2L,
                   min_leaf_n = 2L, a = a, b = b, seed = 64L)
    expect_true(all(fit$diagnostics$leaf_count_trace == 1L))
    expect_equal(fit$posterior$log_target_normalizer, expected,
                 tolerance = 1e-12)
    expect_equal(fit$posterior$log_relative_normalizer, 0, tolerance = 1e-12)
    expect_true(is.na(fit$posterior$log_evidence))
    expect_output(print(summary(fit)), "Log target normalizer")
    improper <- ppt_fit(x, region, gating = "hard", sampler = "smc",
                        engine = engine, particles = 10L, max_depth = 2L,
                        min_leaf_n = 2L, a = a, b = 0, seed = 64L)
    expect_true(is.na(improper$posterior$log_evidence))
    expect_true(is.na(improper$posterior$log_target_normalizer))
    expect_equal(improper$posterior$log_relative_normalizer, 0,
                 tolerance = 1e-12)
    expect_output(print(summary(improper)), "Log relative normalizer")
  }
})

test_that("SMC normalizer preserves tree weights when an axis has no valid cut", {
  x <- cbind(c(0.1, 0.9), c(0, 0))
  region <- matrix(c(0, 1, 0, 1), nrow = 2L, byrow = TRUE)

  root <- audit_log_q(2, 1, 0.5, 0.1)
  scores <- c(log(0.5) + root,
              log(0.25) + 2 * audit_log_q(1, 0.5, 0.5, 0.1))
  expected <- max(scores) + log(sum(exp(scores - max(scores))))
  for (engine in c("shared", "dense")) {
    fit <- ppt_fit(x, region, gating = "hard", sampler = "smc",
                   engine = engine, particles = 30L, max_depth = 1L,
                   a = 0.5, b = 0.1, seed = 87L)
    expect_equal(fit$posterior$log_target_normalizer, expected,
                 tolerance = 1e-12)
    expect_equal(fit$posterior$log_relative_normalizer, expected - root,
                 tolerance = 1e-12)
    expect_true(is.na(fit$posterior$log_evidence))
  }
})

test_that("dense SMC uses the requested number of cut candidates", {
  x <- matrix(seq(0.1, 0.9, length.out = 9L), ncol = 1L)
  region <- matrix(c(0, 1), nrow = 1L)
  candidate_counts <- c(2L, 20L)
  expected_cuts <- list(x[8L, 1L], x[2:8, 1L])
  for (k in seq_along(candidate_counts)) {
    fit <- ppt_fit(x, region, gating = "hard", sampler = "smc",
                   engine = "dense", particles = 2000L, max_depth = 1L,
                   cut_candidates = candidate_counts[k],
                   a = 0.5, b = 0.1, seed = 31L)
    actions <- vapply(fit$posterior$tree_draws,
                      audit_root_action, character(1))
    expect_setequal(unique(actions),
                     c("stop", sprintf("1:%.12g", expected_cuts[[k]])))
    expect_equal(fit$control$cut_candidates, candidate_counts[k])
  }
})

test_that("soft diagnostics use the fitting backend's complete gate scan", {
  x <- rbind(c(0.2, 0.3), c(0.7, 0.8))
  region <- matrix(c(0, 1, 0, 1), nrow = 2L, byrow = TRUE)
  for (shared in c(0L, 1L)) {
    args <- list(
      X = x, region = region, a = 0.5, b = 0.1,
      gate = c(3, 3), a_gate = 2, b_gate = 1,
      sd_gate = 0.15, gate_min = 0, gate_shared = shared,
      alpha = 0.5, eta = 2, Dmax = 0L, nmin = 1L,
      iters = 100L, burn = 10L, thin = 1L, nmove = 0L, ncc = 0L,
      cut_mode = 0L, ncand = 3L, update_gate = 1L, gate_family = 0L
    )
    set.seed(93L)
    diagnostic <- do.call(poistree:::ppstree_diag, c(args, list(mon = x)))
    set.seed(93L)
    fitting <- do.call(poistree:::ppstree_multi, c(args, list(
      grid = x, Xtest = matrix(numeric(), nrow = 0L, ncol = 2L),
      chains = 1L, verbose = 0L
    )))
    expect_equal(diagnostic$gate, fitting$state_gate, tolerance = 1e-12)
    expect_equal(diagnostic$mon, t(fitting$draws), tolerance = 1e-12)
    expect_equal(as.numeric(diagnostic$logdens),
                 as.numeric(fitting$loglik_draws), tolerance = 1e-12)
    expect_equal(diagnostic$accept, fitting$accept, tolerance = 1e-12)
    expect_equal(diagnostic$gate_accept, fitting$gate_accept, tolerance = 1e-12)
    expect_true(all(apply(diagnostic$gate, 2L, stats::var) > 0))
    if (shared == 1L) {
      expect_equal(diagnostic$gate[, 1L], diagnostic$gate[, 2L])
    }
  }
})
