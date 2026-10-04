test_that("hard leaf SMC uses the unified ppt API", {
  set.seed(20)
  x <- matrix(runif(100), ncol = 2)
  colnames(x) <- c("x", "y")
  region <- matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE)
  grid <- matrix(runif(30), ncol = 2)
  test <- matrix(runif(10), ncol = 2)

  fit <- ppt_fit(
    x, region,
    gating = "hard", scales = "leaf", sampler = "smc",
    predict_at = grid, test = test,
    a = 0.6, b = 0.02, resample_thresh = 0.7,
    max_depth = 3, min_leaf_n = 3, particles = 40, seed = 21
  )

  expect_s3_class(fit, "ppt")
  expect_identical(fit$model$label, "PPT")
  expect_identical(fit$model$sampler, "smc")
  expect_length(ppt_predict(fit), nrow(grid))
  expect_length(
    ppt_predict(fit, newdata = grid[c(3, 1), , drop = FALSE]),
    2L
  )
  expect_equal(nrow(ppt_predict(fit, type = "interval")), nrow(grid))
  expect_true(all(is.finite(ppt_predict(fit))))
  expect_true(is.finite(as.numeric(ppt_logLik(fit))))
  expect_true(is.finite(as.numeric(ppt_lppd(fit))))
  expect_true(is.na(fit$posterior$log_evidence))
  expect_true(is.finite(fit$posterior$log_target_normalizer))
  expect_output(print(ppt_summary(fit)), "Log target normalizer")

  diagnostics <- ppt_diagnostics(fit)
  expect_s3_class(diagnostics, "ppt_diagnostics")
  expect_true(is.finite(diagnostics$particle_ess))
  expect_gte(diagnostics$particle_ess, 1)
  expect_lte(diagnostics$particle_ess, fit$control$particles)
  expect_true(length(diagnostics$ess_history) > 0L)
  expect_true(is.finite(diagnostics$unique_trees))
  expect_equal(dim(fit$prediction$draws),
               c(nrow(grid), fit$control$particles))
  expect_length(fit$posterior$tree_draws, fit$control$particles)
  expect_length(diagnostics$leaf_count_trace, fit$control$particles)
  expect_length(diagnostics$max_depth_trace, fit$control$particles)
  expect_true(length(diagnostics$log_evidence_increment) > 0L)
  expect_true(length(diagnostics$log_evidence_running) > 0L)
  expect_true(is.finite(fit$posterior$mean_leaves))
  expect_true(is.finite(fit$posterior$mean_max_depth))
  expect_true(is.finite(ppt_integral(fit)))
  expect_length(ppt_integral(fit, "draws"), fit$control$particles)
  expect_equal(fit$prior$intensity, list(shape = 0.6, rate = 0.02))
  expect_equal(fit$control$resample_thresh, 0.7)

  marginal <- ppt_marginal(fit, "x", grid = c(0.2, 0.8))
  marginal_draws <- ppt_marginal(
    fit, "x", grid = c(0.2, 0.8), type = "draws"
  )
  expect_s3_class(marginal, "ppt_marginal")
  expect_named(
    marginal,
    c("variable", "value", "mean", "median", "lower", "upper")
  )
  expect_equal(dim(marginal_draws), c(2L, fit$posterior$draws))
  expect_equal(
    marginal$mean,
    as.numeric(marginal_draws %*% attr(marginal_draws, "weights"))
  )

  marginal_integral <- ppt_marginal(
    fit, "x", grid = c(0.2, 0.8), average = FALSE, type = "draws"
  )
  y_width <- diff(region[2L, ])
  expect_equal(as.numeric(marginal_integral),
               as.numeric(marginal_draws) * y_width)
})

test_that("hard leaf SMC validates its controls", {
  x <- matrix(seq(0.1, 0.9, length.out = 20), ncol = 1)
  region <- matrix(c(0, 1), nrow = 1)

  expect_error(
    ppt_fit(
      x, region,
      gating = "hard", scales = "leaf", sampler = "smc",
      particles = 1
    ),
    "Invalid SMC or tree controls"
  )
  expect_error(
    ppt_fit(
      x, region,
      gating = "hard", scales = "leaf", sampler = "smc",
      resample_thresh = 0
    ),
    "resample_thresh"
  )
})

test_that("hard leaf SMC retains the requested resampling threshold", {
  x <- matrix(seq(0.05, 0.95, length.out = 20L), ncol = 1L)
  region <- matrix(c(0, 1), nrow = 1L)
  raw <- poistree:::PPT_fit_SMC(
    x, matrix(0.5, ncol = 1L), region,
    max_depth = 2L, P = 8L, min_leaf_n = 1L,
    resample_thresh = 0.73, a = 0.5, b = 0,
    max_aspect_ratio = Inf
  )
  expect_identical(as.numeric(raw$resample_thresh), 0.73)
})

test_that("hard leaf SMC retains a supplied nonunit region", {
  x <- matrix(seq(2.05, 3.95, length.out = 20L), ncol = 1L)
  region <- matrix(c(2, 4), nrow = 1L)
  raw <- poistree:::PPT_fit_SMC(
    x, matrix(3, ncol = 1L), region,
    max_depth = 1L, P = 6L, min_leaf_n = 1L,
    resample_thresh = 0.5, a = 0.5, b = 0,
    max_aspect_ratio = Inf
  )
  expect_equal(as.matrix(raw$particle[[1L]][[1L]]$region), region)
})

test_that("legacy PPT.SMC wrapper is absent", {
  namespace <- asNamespace("poistree")
  expect_false("PPT.SMC" %in% getNamespaceExports("poistree"))
  expect_false(exists("PPT.SMC", envir = namespace, inherits = FALSE))
  expect_false("SMC.diagnostics" %in% getNamespaceExports("poistree"))
  expect_false("SMC.stability" %in% getNamespaceExports("poistree"))
  expect_false(exists("SMC.diagnostics", envir = namespace, inherits = FALSE))
  expect_false(exists("SMC.stability", envir = namespace, inherits = FALSE))
  expect_true(exists("PPT_fit_SMC", envir = namespace, inherits = FALSE))
})
