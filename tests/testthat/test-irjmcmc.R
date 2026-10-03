test_that("iRJ-MCMC routes MPPT and S-MPPT through a common output schema", {
  x <- cbind(
    x1 = seq(0.04, 0.96, length.out = 24),
    x2 = rep(c(0.15, 0.85), 12)
  )
  region <- matrix(c(0, 1, 0, 1), nrow = 2, byrow = TRUE)
  common <- list(
    x = x, region = region, scales = "multiscale",
    scale_prior = "independent", sampler = "irjmcmc",
    max_depth = 2, min_leaf_n = 1,
    chains = 1, iter = 40, burn = 10, thin = 2,
    tree_moves = 2, change_moves = 1, cut_candidates = 6,
    update_hyper = FALSE, update_gate = FALSE,
    verbose = FALSE
  )

  fits <- list(
    MPPT = do.call(
      ppt_fit,
      c(common, list(gating = "hard", seed = 51))
    ),
    `S-MPPT logistic` = do.call(
      ppt_fit,
      c(common, list(
        gating = "soft", gate_family = "logistic", gate = 10,
        seed = 53
      ))
    ),
    `S-MPPT compact` = do.call(
      ppt_fit,
      c(common, list(
        gating = "soft", gate_family = "compact", gate = 10,
        seed = 54
      ))
    )
  )

  for (fit in fits) {
    expect_s3_class(fit, "ppt")
    expect_identical(fit$model$sampler, "irjmcmc")
    expect_identical(fit$model$algorithm, "iRJ-MCMC")
    expect_true(fit$control$informed_cut)
    expect_true(all(is.finite(fit$prediction$mean)))
    expect_true(is.finite(fit$posterior$mean_leaves))
    expect_true(is.finite(fit$posterior$mean_max_depth))
    expect_named(
      fit$diagnostics$tree_acceptance,
      c("grow", "prune", "change")
    )
    expect_length(
      fit$posterior$state$nodes,
      fit$control$chains * fit$control$retained_per_chain
    )
  }

  expect_identical(fits$MPPT$model$label, "MPPT")
  expect_identical(
    fits[["S-MPPT logistic"]]$model$label,
    "S-MPPT"
  )
  expect_identical(
    fits[["S-MPPT compact"]]$model$label,
    "S-MPPT"
  )
  expect_true(fits$MPPT$control$conditional_labels)
  expect_true(fits[["S-MPPT logistic"]]$control$conditional_labels)
  expect_true(fits[["S-MPPT compact"]]$control$conditional_labels)
  expect_true(all(vapply(fits, function(fit) {
    is.null(fit$control$proposal_kernel)
  }, logical(1))))
  expect_identical(
    ppt_diagnostics(fits$MPPT)$algorithm,
    "iRJ-MCMC"
  )
})

test_that("iRJ-MCMC always uses the combined kernel and validates its scope", {
  x <- matrix(seq(0.05, 0.95, length.out = 12), ncol = 1)
  region <- matrix(c(0, 1), nrow = 1)
  controls <- list(
    x = x, region = region, gating = "hard", scales = "multiscale",
    scale_prior = "independent", sampler = "irjmcmc",
    max_depth = 1, chains = 1, iter = 20, burn = 5, thin = 2,
    tree_moves = 1, change_moves = 0, cut_candidates = 4,
    update_hyper = FALSE, update_gate = FALSE,
    verbose = FALSE, seed = 61
  )
  fit <- do.call(ppt_fit, controls)
  expect_null(fit$control$proposal_kernel)
  expect_true(fit$control$informed_cut)
  expect_true(fit$control$conditional_labels)

  markov <- controls
  markov$scale_prior <- "markov"
  expect_error(do.call(ppt_fit, markov), "independent-scale")

  leaf <- controls
  leaf$scales <- "leaf"
  expect_error(do.call(ppt_fit, leaf), "independent-scale")

  obsolete_option <- controls
  obsolete_option$proposal_kernel <- "combined"
  expect_error(
    do.call(ppt_fit, obsolete_option),
    "unused argument"
  )
})

test_that("iRJ-MCMC depth-one grow and prune ratios satisfy detailed balance", {
  x <- cbind(
    seq(0.05, 0.95, length.out = 12),
    rep(c(0.2, 0.8), 6)
  )
  region <- matrix(c(0, 1, 0, 1), nrow = 2, byrow = TRUE)
  configurations <- list(
    list(soft = 0L, gate_mode = 0L),
    list(soft = 1L, gate_mode = 1L),
    list(soft = 1L, gate_mode = 2L)
  )

  for (configuration in configurations) {
    check <- mppstree_irj_balance_check(
      x, region, configuration$soft, rep(12, 2),
      configuration$gate_mode,
      0.7, 2, 0.2, 0.8, 1, 1L, 7L
    )
    expect_lt(abs(check$grow_direct_error), 1e-10)
    expect_lt(abs(check$grow_prune_cycle_error), 1e-10)
    expect_gt(check$cut_probability, 0)
    expect_lte(check$cut_probability, 1)
    expect_gt(check$label_probability, 0)
    expect_lte(check$label_probability, 1)
  }
})
