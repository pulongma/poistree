# Posterior state store: exact post-hoc reproduction of every retained
# intensity draw at arbitrary locations, and the features built on it
# (ppt_lambda, post-hoc ppt_predict / ppt_lppd / ppt_marginal).

make_data <- function(seed, n = 80) {
  set.seed(seed)
  x <- matrix(runif(2 * n), ncol = 2)
  colnames(x) <- c("x", "y")
  x
}
region2 <- matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE)

fit_args <- function(x, grid, ...) {
  c(list(
    x = x, region = region2, predict_at = grid,
    max_depth = 3, min_leaf_n = 3,
    chains = 2, iter = 60, burn = 20, thin = 2,
    cut_candidates = 5, verbose = FALSE
  ), list(...))
}

expect_state_matches_prediction <- function(fit, grid, tol = 1e-10) {
  draws <- poistree:::.ppt_state_eval(fit, grid)
  stored <- as.matrix(fit$prediction$draws)
  expect_equal(dim(draws), dim(stored))
  expect_equal(draws, stored, tolerance = tol, ignore_attr = TRUE)
}

test_that("heap state reproduces stored draws: soft leaf, logistic", {
  x <- make_data(101)
  grid <- matrix(runif(24), ncol = 2)
  fit <- do.call(ppt_fit, fit_args(
    x, grid,
    gating = "soft", scales = "leaf", sampler = "rjmcmc",
    gate = 12, update_gate = TRUE
  ))
  expect_identical(fit$posterior$state$mode, "heap")
  expect_length(fit$posterior$state$chain, fit$posterior$draws)
  expect_state_matches_prediction(fit, grid)
})

test_that("heap state reproduces stored draws: soft leaf, compact", {
  x <- make_data(102)
  grid <- matrix(runif(24), ncol = 2)
  fit <- do.call(ppt_fit, fit_args(
    x, grid,
    gating = "soft", scales = "leaf", sampler = "rjmcmc",
    gate_family = "compact", gate = 3, update_gate = TRUE
  ))
  expect_state_matches_prediction(fit, grid)
})

test_that("heap state reproduces stored draws: hard multiscale", {
  x <- make_data(103)
  grid <- matrix(runif(24), ncol = 2)
  fit <- do.call(ppt_fit, fit_args(
    x, grid,
    gating = "hard", scales = "multiscale", sampler = "rjmcmc"
  ))
  expect_identical(fit$posterior$state$mode, "heap")
  expect_state_matches_prediction(fit, grid)
})

test_that("heap state reproduces stored draws: soft multiscale", {
  x <- make_data(104)
  grid <- matrix(runif(24), ncol = 2)
  fit <- do.call(ppt_fit, fit_args(
    x, grid,
    gating = "soft", scales = "multiscale", sampler = "rjmcmc",
    gate = 12, update_gate = TRUE
  ))
  expect_state_matches_prediction(fit, grid)
})

test_that("leaf-box state reproduces stored draws: hard leaf RJ and SMC", {
  x <- make_data(105)
  grid <- matrix(runif(24), ncol = 2)
  rj <- ppt_fit(
    x, region2,
    gating = "hard", scales = "leaf", sampler = "rjmcmc",
    predict_at = grid, max_depth = 3, min_leaf_n = 3,
    chains = 2, iter = 60, burn = 20, cut_candidates = 5,
    prediction_draws = 15, seed = 7, verbose = FALSE
  )
  expect_identical(rj$posterior$state$mode, "leafbox")
  expect_length(rj$posterior$state$chain, ncol(rj$prediction$draws))
  expect_state_matches_prediction(rj, grid)

  smc <- ppt_fit(
    x, region2,
    gating = "hard", scales = "leaf", sampler = "smc",
    predict_at = grid, max_depth = 3, min_leaf_n = 3,
    particles = 30, seed = 8
  )
  expect_identical(smc$posterior$state$mode, "leafbox")
  expect_state_matches_prediction(smc, grid)
})

test_that("ppt_lambda summarizes on a grid and matches raw draws", {
  x <- make_data(106)
  grid <- matrix(runif(20), ncol = 2)
  fit <- do.call(ppt_fit, fit_args(
    x, grid,
    gating = "soft", scales = "leaf", sampler = "rjmcmc",
    gate = 12, update_gate = FALSE
  ))
  surface <- ppt_lambda(fit, n = 6)
  expect_s3_class(surface, "data.frame")
  expect_equal(nrow(surface), 36L)
  expect_named(surface, c("x", "y", "mean", "median", "lower", "upper"))
  expect_true(all(is.finite(surface$mean)))
  expect_true(all(surface$lower <= surface$median + 1e-12))
  expect_true(all(surface$median <= surface$upper + 1e-12))

  raw <- ppt_lambda(fit, at = as.matrix(surface[, 1:2]), type = "draws")
  expect_equal(as.numeric(raw$draws %*% raw$weights), surface$mean,
               tolerance = 1e-12)

  # one-dimensional grid path
  x1 <- matrix(runif(40), ncol = 1)
  fit1 <- ppt_fit(
    x1, matrix(c(0, 1), nrow = 1),
    gating = "hard", scales = "multiscale", sampler = "rjmcmc",
    max_depth = 3, min_leaf_n = 3, chains = 1, iter = 40, burn = 10,
    cut_candidates = 4, verbose = FALSE
  )
  line <- ppt_lambda(fit1, n = 15)
  expect_equal(nrow(line), 15L)
  expect_error(ppt_lambda(fit1, n = 1), "at least 2")
})

test_that("ppt_predict evaluates locations outside predict_at", {
  x <- make_data(107)
  grid <- matrix(runif(20), ncol = 2)
  fit <- do.call(ppt_fit, fit_args(
    x, grid,
    gating = "soft", scales = "multiscale", sampler = "rjmcmc",
    gate = 12, update_gate = FALSE
  ))
  fresh <- matrix(runif(10), ncol = 2)
  mean_new <- ppt_predict(fit, newdata = fresh, type = "mean")
  expect_length(mean_new, nrow(fresh))
  expect_true(all(is.finite(mean_new)))
  lam <- ppt_lambda(fit, at = fresh)
  expect_equal(as.numeric(mean_new), lam$mean, tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(as.numeric(ppt_predict(fit, newdata = fresh)), lam$median,
               tolerance = 1e-12, ignore_attr = TRUE)
  interval <- ppt_predict(fit, newdata = fresh, type = "interval")
  expect_equal(interval$lower, lam$lower, tolerance = 1e-12,
               ignore_attr = TRUE)
  # rows present in predict_at still take the stored fast path
  expect_equal(
    as.numeric(ppt_predict(fit, newdata = grid[c(4, 2), , drop = FALSE])),
    as.numeric(ppt_predict(fit)[c(4, 2)]),
    tolerance = 1e-12
  )
})

test_that("post-hoc ppt_lppd(test) equals fit-time lppd", {
  x <- make_data(108)
  test <- matrix(runif(16), ncol = 2)
  grid <- matrix(runif(12), ncol = 2)

  soft <- do.call(ppt_fit, fit_args(
    x, grid,
    gating = "soft", scales = "leaf", sampler = "rjmcmc",
    gate = 12, update_gate = TRUE, test = test
  ))
  expect_equal(as.numeric(ppt_lppd(soft, test = test)),
               as.numeric(ppt_lppd(soft)), tolerance = 1e-8)

  hard <- ppt_fit(
    x, region2,
    gating = "hard", scales = "leaf", sampler = "rjmcmc",
    predict_at = grid, test = test, max_depth = 3, min_leaf_n = 3,
    chains = 2, iter = 60, burn = 20, cut_candidates = 5,
    prediction_draws = 15, seed = 9, verbose = FALSE
  )
  expect_equal(as.numeric(ppt_lppd(hard, test = test)),
               as.numeric(ppt_lppd(hard)), tolerance = 1e-8)

  expect_error(ppt_lppd(soft, test = matrix(2, 1, 2)), "region")
})

test_that("plot.ppt renders from state without predict_at", {
  x <- make_data(110)
  fit <- ppt_fit(
    x, region2,
    gating = "soft", scales = "multiscale", sampler = "rjmcmc",
    gate = 12, update_gate = FALSE,
    max_depth = 3, min_leaf_n = 3, chains = 1, iter = 40, burn = 10,
    cut_candidates = 4, verbose = FALSE
  )
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(fit, n = 12, points = TRUE))
  expect_invisible(plot(fit, dims = 2L, n = 12, points = TRUE))
  expect_invisible(plot(fit, type = "mean", dims = c(2L, 1L), n = 8))
  expect_error(plot(fit, dims = c(1L, 1L)), "distinct")
  expect_error(plot(fit, n = 1), "at least 2")
})

test_that("ppt_marginal works without pre-listed reference rows", {
  x <- make_data(109)
  grid <- matrix(runif(12), ncol = 2)
  fit <- do.call(ppt_fit, fit_args(
    x, grid,
    gating = "soft", scales = "leaf", sampler = "rjmcmc",
    gate = 12, update_gate = FALSE
  ))
  reference <- matrix(runif(30), ncol = 2)
  values <- seq(0.1, 0.9, length.out = 5)
  marg <- ppt_marginal(fit, "x", grid = values, reference = reference)
  expect_s3_class(marg, "ppt_marginal")
  expect_equal(marg$value, values)
  expect_true(all(is.finite(marg$mean)))

  # agrees with a manual average of state-evaluated substitution designs
  draws <- ppt_marginal(fit, "x", grid = values, reference = reference,
                        type = "draws")
  substituted <- reference
  substituted[, 1L] <- values[3L]
  manual <- colMeans(poistree:::.ppt_state_eval(fit, substituted))
  expect_equal(as.numeric(draws[3L, ]), as.numeric(manual),
               tolerance = 1e-10)
})
