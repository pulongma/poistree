test_that("all fitting backends default to fifty candidate cuts", {
  backend_names <- unname(poistree:::.ppt_backend_registry)
  defaults <- vapply(backend_names, function(name) {
    formals(getFromNamespace(name, "poistree"))$cut_candidates
  }, integer(1))
  expect_identical(unname(defaults), rep(50L, length(backend_names)))
  expect_identical(formals(poistree:::PPT_fit_PG)$cut_grid_n, 50L)
  expect_identical(formals(poistree:::PPT_fit_SMC)$cut_grid_n, 50L)
  expect_identical(formals(poistree:::PPT_fit_MCMC)$cut_grid_n, 50L)
  expect_identical(formals(poistree:::PPT_fit_IMCMC)$cut_grid_n, 50L)
})

test_that("informed hard and soft fits return ordinary posterior draws", {
  x <- matrix(c(0.1, 0.15, 0.3, 0.6, 0.85, 0.9), ncol = 1L)
  region <- matrix(c(0, 1), nrow = 1L)
  for (gating in c("hard", "soft")) {
    args <- list(
      x = x, region = region, gating = gating, sampler = "irjmcmc",
      a = 0.5, b = 0.2, max_depth = 1L, cut_candidates = 4L,
      chains = 1L, iter = 30L, burn = 5L, seed = 507L, verbose = FALSE
    )
    if (gating == "soft") {
      args <- c(args, list(thin = 1L, tree_moves = 1L, change_moves = 1L,
                           update_gate = FALSE))
    } else {
      args <- c(args, list(prediction_draws = 25L))
    }
    fit <- do.call(ppt_fit, args)
    expect_s3_class(fit, "ppt")
    expect_identical(fit$model$sampler, "irjmcmc")
    expect_identical(fit$model$algorithm, "Informed MH")
    expect_identical(fit$control$cut_candidates, 4L)
    expect_length(fit$posterior$particle_weights, 0L)
    evaluated <- ppt_lambda(fit, at = x, type = "draws")
    expect_equal(evaluated$draws, fit$prediction$draws, tolerance = 1e-12)
    expect_equal(evaluated$weights,
                 rep(1 / ncol(evaluated$draws), ncol(evaluated$draws)))
    expect_equal(ppt_predict(fit, type = "mean"),
                 rowMeans(evaluated$draws), tolerance = 1e-12)
    expect_true(all(is.finite(evaluated$draws)))
    expect_true(all(evaluated$draws > 0))
    expect_true(is.finite(as.numeric(ppt_logLik(fit))))
    expect_true(all(is.finite(ppt_integral(fit, "draws"))))
    expect_length(fit$diagnostics$leaf_count_trace, fit$posterior$draws)
    expect_length(fit$diagnostics$max_depth_trace, fit$posterior$draws)
    expect_equal(mean(fit$diagnostics$leaf_count_trace), fit$posterior$mean_leaves)
    expect_equal(mean(fit$diagnostics$max_depth_trace), fit$posterior$mean_max_depth)
    expect_identical(ppt_diagnostics(fit)$sampler, "irjmcmc")
    expect_output(print(fit), "Informed MH")
  }
})

test_that("omitted cut counts match explicit fifty for every sampler", {
  x <- matrix(c(0.1, 0.2, 0.4, 0.7, 0.8, 0.9), ncol = 1L)
  region <- matrix(c(0, 1), nrow = 1L)
  for (key in names(poistree:::.ppt_backend_registry)) {
    components <- strsplit(key, ":", fixed = TRUE)[[1L]]
    args <- list(
      x = x, region = region, gating = components[1L],
      sampler = components[3L], a = 0.5, b = 0.2,
      max_depth = 1L, seed = 65L
    )
    if (components[3L] == "smc") {
      args <- c(args, list(particles = 8L))
    } else {
      args <- c(args, list(chains = 1L, iter = 12L, burn = 2L,
                           verbose = FALSE))
      if (components[3L] == "pgas") {
        args <- c(args, list(particles = 4L))
      } else if (components[1L] == "soft") {
        args <- c(args, list(thin = 1L, tree_moves = 1L, change_moves = 1L))
      }
    }
    implicit <- do.call(ppt_fit, args)
    explicit <- do.call(ppt_fit, c(args, list(cut_candidates = 50L)))
    expect_identical(implicit$control$cut_candidates, 50L)
    expect_identical(implicit$prediction, explicit$prediction)
    expect_identical(implicit$posterior, explicit$posterior)
    expect_identical(implicit$diagnostics, explicit$diagnostics)
  }
})

test_that("hard Particle Gibbs passes its explicit candidate count", {
  x <- matrix(seq(0.1, 0.9, length.out = 9L), ncol = 1L)
  region <- matrix(c(0, 1), nrow = 1L)
  fit <- ppt_fit(
    x, region, gating = "hard", sampler = "pgas",
    a = 0.5, b = 0.1, cut_candidates = 2L, max_depth = 1L,
    chains = 1L, particles = 5L, iter = 120L, burn = 20L,
    seed = 74L, verbose = FALSE
  )
  expect_identical(fit$control$cut_candidates, 2L)
  split_roots <- Filter(function(tree) !isTRUE(tree[[1L]]$is_leaf),
                         fit$posterior$tree_draws)
  expect_gt(length(split_roots), 0L)
  expect_true(all(vapply(split_roots, function(tree) {
    abs(tree[[1L]]$L - x[8L, 1L]) < 1e-12
  }, logical(1))))
  expect_error(ppt_fit(x, region, gating = "hard", sampler = "pgas",
                       cut_candidates = 0L), "Invalid Particle-Gibbs")
})
