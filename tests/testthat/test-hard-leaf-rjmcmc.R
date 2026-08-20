test_that("hard leaf RJ-MCMC uses the unified ppt API", {
  set.seed(10)
  x <- matrix(runif(80), ncol = 2)
  colnames(x) <- c("x", "y")
  region <- matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE)
  grid <- matrix(runif(30), ncol = 2)
  test <- matrix(runif(10), ncol = 2)

  fit <- ppt_fit(
    x, region,
    gating = "hard", scales = "leaf", sampler = "rjmcmc",
    predict_at = grid, test = test,
    max_depth = 3, min_leaf_n = 3,
    chains = 1, iter = 50, burn = 10,
    cut_candidates = 5, prediction_draws = 10,
    seed = 11, verbose = FALSE
  )

  expect_s3_class(fit, "ppt")
  expect_identical(fit$model$label, "PPT")
  expect_identical(fit$model$gating, "hard")
  expect_identical(fit$model$scales, "leaf")
  expect_identical(fit$model$sampler, "rjmcmc")
  expect_length(ppt_predict(fit), nrow(grid))
  expect_equal(nrow(ppt_predict(fit, type = "interval")), nrow(grid))
  expect_true(all(is.finite(ppt_predict(fit))))
  expect_true(is.finite(as.numeric(ppt_logLik(fit))))
  expect_true(is.finite(as.numeric(ppt_lppd(fit))))
  expect_true(is.na(ppt_evidence(fit, warn = FALSE)))
  expect_true(is.finite(fit$posterior$mean_leaves))
  expect_true(is.finite(fit$posterior$mean_max_depth))
  expect_true(is.finite(ppt_integral(fit)))
  expect_length(ppt_integral(fit, "draws"), fit$posterior$draws)
  expect_named(
    fit$diagnostics$tree_acceptance,
    c("grow", "prune", "change")
  )
})

test_that("hard leaf PGAS uses the unified ppt API", {
  set.seed(20)
  x <- matrix(runif(40), ncol = 2)
  colnames(x) <- c("x", "y")
  region <- matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE)
  grid <- matrix(runif(12), ncol = 2)
  test <- matrix(runif(6), ncol = 2)

  fit <- ppt_fit(
    x, region,
    gating = "hard", scales = "leaf", sampler = "pgas",
    predict_at = grid, test = test,
    max_depth = 2, min_leaf_n = 2,
    particles = 10, chains = 1, iter = 6, burn = 2,
    seed = 21, verbose = FALSE
  )

  expect_s3_class(fit, "ppt")
  expect_identical(fit$model$label, "PPT")
  expect_identical(fit$model$gating, "hard")
  expect_identical(fit$model$scales, "leaf")
  expect_identical(fit$model$sampler, "pgas")
  expect_length(ppt_predict(fit), nrow(grid))
  expect_true(all(is.finite(ppt_predict(fit))))
  expect_true(is.finite(as.numeric(ppt_logLik(fit))))
  expect_true(is.finite(as.numeric(ppt_lppd(fit))))
  expect_true(is.na(ppt_evidence(fit, warn = FALSE)))
  expect_true(is.finite(fit$posterior$mean_leaves))
  expect_true(is.finite(fit$posterior$mean_max_depth))
  expect_true(is.finite(ppt_integral(fit)))
  expect_length(ppt_integral(fit, "draws"), fit$posterior$draws)
  expect_true(is.finite(fit$diagnostics$particle_ess))
})

test_that("all PPT samplers return the same unified schema", {
  x <- matrix(seq(0.05, 0.95, length.out = 20), ncol = 1)
  region <- matrix(c(0, 1), nrow = 1)
  grid <- matrix(seq(0.1, 0.9, length.out = 5), ncol = 1)

  smc <- ppt_fit(
    x, region, gating = "hard", scales = "leaf", sampler = "smc",
    predict_at = grid, max_depth = 2, min_leaf_n = 2,
    particles = 10, seed = 12
  )
  rjmcmc <- ppt_fit(
    x, region, gating = "hard", scales = "leaf", sampler = "rjmcmc",
    predict_at = grid, max_depth = 2, min_leaf_n = 2,
    chains = 1, iter = 30, burn = 10,
    cut_candidates = 5, prediction_draws = 5,
    seed = 13, verbose = FALSE
  )
  pgas <- ppt_fit(
    x, region, gating = "hard", scales = "leaf", sampler = "pgas",
    predict_at = grid, max_depth = 2, min_leaf_n = 2,
    particles = 8, chains = 1, iter = 5, burn = 2,
    seed = 14, verbose = FALSE
  )

  expect_identical(names(rjmcmc), names(smc))
  expect_identical(names(pgas), names(smc))
  expect_identical(names(rjmcmc$posterior), names(smc$posterior))
  expect_identical(names(pgas$posterior), names(smc$posterior))
  expect_identical(names(rjmcmc$diagnostics), names(smc$diagnostics))
  expect_identical(names(pgas$diagnostics), names(smc$diagnostics))
})

test_that("non-PPT component combinations are reserved", {
  x <- matrix(seq(0.1, 0.9, length.out = 20), ncol = 1)
  region <- matrix(c(0, 1), nrow = 1)

  expect_error(
    ppt_fit(
      x, region, gating = "soft", scales = "leaf",
      sampler = "rjmcmc"
    ),
    "reserved for future model backends"
  )
  expect_error(
    ppt_fit(
      x, region, gating = "hard", scales = "multiscale",
      sampler = "rjmcmc"
    ),
    "reserved for future model backends"
  )
  expect_error(
    ppt_fit(
      x, region, gating = "soft", scales = "leaf",
      sampler = "pgas"
    ),
    "available only for the hard-gated, terminal-leaf PPT"
  )
  expect_error(
    ppt_fit(
      x, region, gating = "hard", scales = "multiscale",
      sampler = "pgas"
    ),
    "available only for the hard-gated, terminal-leaf PPT"
  )
})

test_that("soft-gating native functions are absent", {
  native_names <- ls(asNamespace("poistree"), all.names = TRUE)
  expect_false(any(grepl("ppstree|soft", native_names, ignore.case = TRUE)))
  expect_false("PPT.MCMC" %in% getNamespaceExports("poistree"))
  expect_false("PPT.lppd" %in% getNamespaceExports("poistree"))
  expect_false(exists("PPT.MCMC", envir = asNamespace("poistree"), inherits = FALSE))
  expect_false(exists("PPT.lppd", envir = asNamespace("poistree"), inherits = FALSE))
  expect_false("PPT.PG" %in% getNamespaceExports("poistree"))
  expect_false(exists("PPT.PG", envir = asNamespace("poistree"), inherits = FALSE))
  expect_true(any(grepl("PPT_fit_MCMC", native_names, fixed = TRUE)))
  expect_true(any(grepl("PPT_fit_PG", native_names, fixed = TRUE)))
})

test_that("hard RJ-MCMC snaps boundary roundoff safely", {
  x <- matrix(c(0.2, 0.8, 1 + 5e-14), ncol = 1)
  region <- matrix(c(0, 1), nrow = 1)
  fit <- ppt_fit(
    x, region, gating = "hard", scales = "leaf", sampler = "rjmcmc",
    max_depth = 1, min_leaf_n = 1,
    chains = 1, iter = 20, burn = 5,
    cut_candidates = 3, prediction_draws = 5,
    verbose = FALSE
  )
  expect_equal(max(fit$data$x), 1)
})
