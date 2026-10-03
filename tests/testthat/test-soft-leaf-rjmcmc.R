test_that("defaults select S-PPT; hard gating selects PPT", {
  x <- matrix(seq(0.1, 0.9, length.out = 20), ncol = 1)
  region <- matrix(c(0, 1), nrow = 1)
  soft <- ppt_fit(
    x, region, chains = 1, iter = 20, burn = 5, thin = 2, max_depth = 1,
    tree_moves = 0, change_moves = 0, cut_candidates = 3,
    update_gate = FALSE, verbose = FALSE
  )
  hard <- ppt_fit(
    x, region, gating = "hard", chains = 1, iter = 20, burn = 5,
    max_depth = 1, cut_candidates = 3, prediction_draws = 5, verbose = FALSE
  )
  expect_identical(soft$model$label, "S-PPT")
  expect_identical(hard$model$label, "PPT")
  expect_identical(soft$model$sampler, "rjmcmc")
  expect_identical(hard$model$sampler, "rjmcmc")
  expect_identical(soft$model$scales, "leaf")
  expect_identical(hard$model$scales, "leaf")
  expect_identical(soft$control$min_leaf_n, 1L)
  expect_identical(hard$control$min_leaf_n, 1L)

  ## The default S-PPT stores enough state for post-hoc evaluation, even
  ## though no `predict_at` locations were supplied during fitting.
  surface <- ppt_lambda(soft, n = 7)
  expect_equal(nrow(surface), 7L)
  expect_true(all(is.finite(surface$mean)))
})

test_that("every fitting backend defaults to one observation per leaf", {
  backend_names <- c(
    ".ppt_fit_hard_leaf_smc",
    ".ppt_fit_hard_leaf_rjmcmc",
    ".ppt_fit_hard_leaf_pgas",
    ".ppt_fit_soft_leaf_rjmcmc"
  )
  defaults <- vapply(backend_names, function(name) {
    formals(getFromNamespace(name, "poistree"))$min_leaf_n
  }, integer(1))
  expect_identical(unname(defaults), rep(1L, length(backend_names)))
})

test_that("soft leaf RJ-MCMC returns a ppt object", {
  set.seed(10)
  x <- matrix(runif(80), ncol = 2)
  colnames(x) <- c("x", "y")
  region <- matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE)
  grid <- matrix(runif(30), ncol = 2)

  fit <- ppt_fit(
    x, region,
    gating = "soft", scales = "leaf", sampler = "rjmcmc",
    predict_at = grid,
    gate = 12, update_gate = FALSE,
    max_depth = 3, min_leaf_n = 3,
    chains = 1, iter = 40, burn = 10, thin = 2,
    tree_moves = 1, change_moves = 1, cut_candidates = 5,
    verbose = FALSE
  )

  expect_s3_class(fit, "ppt")
  expect_identical(fit$model$label, "S-PPT")
  expect_identical(fit$model$gate_family, "logistic")
  expect_identical(fit$prior$gate$family, "logistic")
  expect_length(ppt_predict(fit), nrow(grid))
  expect_equal(nrow(ppt_predict(fit, type = "interval")), nrow(grid))
  expect_s3_class(ppt_logLik(fit), "logLik")
  expect_s3_class(ppt_summary(fit), "summary.ppt")
  expect_s3_class(ppt_diagnostics(fit), "ppt_diagnostics")
  expect_true(all(is.finite(ppt_predict(fit))))
  expect_true(is.finite(fit$posterior$mean_max_depth))
  expect_equal(dim(fit$prediction$draws),
               c(nrow(grid), fit$posterior$draws))
})

test_that("default logistic children partition their parent exactly", {
  region <- matrix(c(0, 1, -1, 2), ncol = 2, byrow = TRUE)
  points <- cbind(
    seq(0, 1, length.out = 101),
    seq(-1, 2, length.out = 101)
  )
  empty_i <- integer()
  empty_d <- numeric()

  root <- poistree:::ppstree_geometry(
    empty_i, empty_d, empty_d, empty_i, points, region, gate = 9
  )
  left <- poistree:::ppstree_geometry(
    1L, 0.4, 1, -1L, points, region, gate = 9
  )
  right <- poistree:::ppstree_geometry(
    1L, 0.4, 1, 1L, points, region, gate = 9
  )

  expect_equal(left$phi + right$phi, root$phi, tolerance = 1e-12)
  expect_equal(left$H + right$H, root$H, tolerance = 1e-12)
})

test_that("compact gate remains available as an explicit option", {
  x <- matrix(seq(0.05, 0.95, length.out = 20), ncol = 1)
  region <- matrix(c(0, 1), nrow = 1)
  fit <- ppt_fit(
    x, region,
    gating = "soft", scales = "leaf", sampler = "rjmcmc",
    gate_family = "compact", update_gate = FALSE,
    max_depth = 1, min_leaf_n = 2,
    chains = 1, iter = 20, burn = 5, thin = 2,
    tree_moves = 1, change_moves = 1, cut_candidates = 4,
    verbose = FALSE
  )
  expect_identical(fit$model$gate_family, "compact")
})

test_that("all unified models use the same output schema", {
  set.seed(12)
  x <- matrix(runif(80), ncol = 2)
  colnames(x) <- c("x", "y")
  region <- matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE)
  grid <- matrix(runif(16), ncol = 2)

  hard_leaf <- ppt_fit(
    x, region,
    gating = "hard", scales = "leaf", sampler = "smc",
    predict_at = grid, max_depth = 2, min_leaf_n = 3,
    particles = 10, seed = 13
  )
  soft_leaf <- ppt_fit(
    x, region,
    gating = "soft", scales = "leaf", sampler = "rjmcmc",
    predict_at = grid, gate = 10, update_gate = FALSE,
    max_depth = 2, min_leaf_n = 3,
    chains = 1, iter = 20, burn = 5, thin = 2,
    tree_moves = 1, change_moves = 1, cut_candidates = 4,
    seed = 14, verbose = FALSE
  )
  fits <- list(hard_leaf, soft_leaf)

  posterior_names <- names(hard_leaf$posterior)
  diagnostic_names <- names(hard_leaf$diagnostics)
  expect_true(all(vapply(
    fits, function(fit) identical(names(fit$posterior), posterior_names),
    logical(1)
  )))
  expect_true(all(vapply(
    fits, function(fit) identical(names(fit$diagnostics), diagnostic_names),
    logical(1)
  )))
  expect_true(all(vapply(
    fits, function(fit) is.finite(fit$posterior$mean_max_depth), logical(1)
  )))
  expect_true(all(vapply(
    fits, function(fit) is.finite(ppt_integral(fit)), logical(1)
  )))
  expect_true(all(vapply(
    fits,
    function(fit) length(ppt_integral(fit, "draws")) == fit$posterior$draws,
    logical(1)
  )))
  expect_true(is.na(hard_leaf$posterior$mean_gate))
})

test_that("soft-tree native identifiers put S before T", {
  native_names <- ls(asNamespace("poistree"), all.names = TRUE)
  expect_true(any(grepl("ppstree", native_names, fixed = TRUE)))
  forbidden_order <- "ppts"
  expect_false(any(vapply(
    forbidden_order,
    function(pattern) any(grepl(pattern, native_names, ignore.case = TRUE)),
    logical(1)
  )))
})

test_that("roundoff at a region boundary is snapped safely", {
  x <- matrix(c(0.2, 0.8, 1 + 5e-14), ncol = 1)
  region <- matrix(c(0, 1), nrow = 1)
  fit <- ppt_fit(
    x, region,
    gating = "soft", scales = "leaf", sampler = "rjmcmc",
    max_depth = 1, min_leaf_n = 1,
    chains = 1, iter = 20, burn = 5, thin = 1,
    tree_moves = 0, change_moves = 0, cut_candidates = 3,
    update_gate = FALSE, verbose = FALSE
  )
  expect_equal(max(fit$data$x), 1)
})

test_that("ppt_sim produces points inside the region", {
  region <- matrix(c(0, 1, -1, 1), ncol = 2, byrow = TRUE)
  simulated <- ppt_sim(
    function(z) rep(5, nrow(z)),
    region = region,
    lambda_max = 5,
    seed = 2
  )
  expect_equal(ncol(simulated), 2)
  expect_true(all(simulated[, 1] >= 0 & simulated[, 1] <= 1))
  expect_true(all(simulated[, 2] >= -1 & simulated[, 2] <= 1))
})
