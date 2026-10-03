test_that("all MPPT cut proposals enforce minimum child counts", {
  x <- matrix(seq(0.05, 0.95, length.out = 10), ncol = 1)
  region <- matrix(c(0, 1), nrow = 1)

  fit <- ppt_fit(
    x, region,
    gating = "hard", scales = "multiscale",
    scale_prior = "independent", sampler = "rjmcmc",
    cut_proposal = "uniform",
    max_depth = 2, min_leaf_n = 10,
    chains = 1, iter = 30, burn = 10, thin = 2,
    tree_moves = 2, change_moves = 1,
    update_hyper = FALSE, update_gate = FALSE,
    verbose = FALSE
  )
  soft_fit <- ppt_fit(
    x, region,
    gating = "soft", scales = "multiscale",
    scale_prior = "independent", sampler = "rjmcmc",
    cut_proposal = "uniform",
    max_depth = 2, min_leaf_n = 10,
    chains = 1, iter = 30, burn = 10, thin = 2,
    tree_moves = 2, change_moves = 1,
    update_hyper = FALSE, update_gate = FALSE,
    verbose = FALSE
  )

  expect_equal(fit$posterior$mean_leaves, 1)
  expect_equal(fit$posterior$mean_max_depth, 0)
  for (candidate in list(fit, soft_fit)) {
    expect_named(
      candidate$diagnostics$tree_acceptance,
      c("grow", "prune", "change")
    )
    expect_equal(
      unname(candidate$diagnostics$tree_acceptance),
      rep(0, 3)
    )
  }
})

test_that("data cut proposals do not split tied observations", {
  x <- matrix(rep(0.5, 20), ncol = 1)
  region <- matrix(c(0, 1), nrow = 1)

  hard_multiscale <- ppt_fit(
    x, region,
    gating = "hard", scales = "multiscale",
    scale_prior = "independent", sampler = "rjmcmc",
    cut_proposal = "data",
    max_depth = 2, min_leaf_n = 10,
    chains = 1, iter = 30, burn = 10, thin = 2,
    tree_moves = 2, change_moves = 1,
    update_hyper = FALSE, update_gate = FALSE,
    verbose = FALSE
  )
  soft_leaf <- ppt_fit(
    x, region,
    gating = "soft", scales = "leaf", sampler = "rjmcmc",
    cut_proposal = "data",
    max_depth = 2, min_leaf_n = 10,
    chains = 1, iter = 30, burn = 10, thin = 2,
    tree_moves = 2, change_moves = 1,
    update_gate = FALSE, verbose = FALSE
  )

  expect_equal(hard_multiscale$posterior$mean_leaves, 1)
  expect_equal(soft_leaf$posterior$mean_leaves, 1)
})

test_that("soft MPPT gate updates are independent of scale hyper-updates", {
  x <- matrix(seq(0.05, 0.95, length.out = 30), ncol = 1)
  region <- matrix(c(0, 1), nrow = 1)
  common <- list(
    x = x, region = region,
    gating = "soft", scales = "multiscale",
    scale_prior = "independent", sampler = "rjmcmc",
    gate = 1, a_gate = 100, b_gate = 5, sd_gate = 0.3,
    max_depth = 0, min_leaf_n = 2,
    chains = 1, iter = 80, burn = 20, thin = 2,
    tree_moves = 0, change_moves = 0,
    update_hyper = FALSE, verbose = FALSE
  )

  fixed_gate <- do.call(
    ppt_fit, c(common, list(update_gate = FALSE, seed = 41))
  )
  moving_gate <- do.call(
    ppt_fit, c(common, list(update_gate = TRUE, seed = 42))
  )

  expect_equal(fixed_gate$posterior$mean_gate, 1)
  expect_equal(fixed_gate$diagnostics$gate_acceptance[[1]], 0)
  expect_gt(moving_gate$diagnostics$gate_acceptance[[1]], 0)
  expect_false(isTRUE(all.equal(
    moving_gate$posterior$mean_gate, 1, tolerance = 1e-8
  )))
})

test_that("default independent scale rate retains the unit-region convention", {
  x <- matrix(seq(0.05, 1.95, length.out = 20), ncol = 1)
  region <- matrix(c(0, 2), nrow = 1)

  fit <- ppt_fit(
    x, region,
    gating = "hard", scales = "multiscale",
    scale_prior = "independent", sampler = "rjmcmc",
    a_xi = 2,
    max_depth = 0, min_leaf_n = 2,
    chains = 1, iter = 20, burn = 5, thin = 2,
    tree_moves = 0, change_moves = 0,
    update_hyper = FALSE, update_gate = FALSE,
    verbose = FALSE
  )

  expect_equal(fit$prior$scales$rate, 2 / 20)
})
