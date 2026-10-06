# PGAS retains its box-free cut prior and root-scaled logistic gates. A local
# gate width would change both the valid cut support and ancestor weights.
# Reject that combination until those changes have their own posterior test.
test_that("soft PGAS rejects node scaling before any RNG draws", {
  x <- matrix(c(0.12, 0.25, 0.38, 0.56, 0.72, 0.91), ncol = 1)
  region <- matrix(c(0, 1), nrow = 1)
  for (ancestor_sampling in c(TRUE, FALSE)) {
    set.seed(803)
    before <- .Random.seed
    expect_error(
      ppt_fit(x, region, sampler = "pgas", gate_scale = "node",
              ancestor_sampling = ancestor_sampling, verbose = FALSE),
      'not supported by soft PGAS.*pcg.*rjmcmc.*irjmcmc'
    )
    expect_identical(.Random.seed, before)
  }
})

test_that("explicit root scaling preserves the PGAS chain and state evaluator", {
  x <- matrix(c(0.12, 0.25, 0.38, 0.56, 0.72, 0.91), ncol = 1)
  region <- matrix(c(0, 1), nrow = 1)
  at <- matrix(c(0.15, 0.45, 0.85), ncol = 1)
  heldout <- matrix(c(0.21, 0.63), ncol = 1)
  arguments <- list(
    x = x, region = region, sampler = "pgas", predict_at = at, test = heldout,
    chains = 1, particles = 3, iter = 16, burn = 4, thin = 2,
    max_depth = 2, cut_grid = list(c(0.3, 0.5, 0.7)),
    seed = 803, verbose = FALSE
  )
  for (proposal in c("laplace", "hard")) {
    arguments$proposal_score <- proposal
    arguments$exact_max <- if (proposal == "hard") 0 else 150
    old_default <- do.call(ppt_fit, arguments)
    end_seed <- .Random.seed
    explicit_root <- do.call(ppt_fit, c(arguments, list(gate_scale = "root")))
    expect_identical(.Random.seed, end_seed)
    expect_identical(explicit_root$posterior, old_default$posterior)
    expect_identical(explicit_root$prediction, old_default$prediction)
    expect_equal(explicit_root$posterior$state$gate_mode, 2L)

    # Version 0.4 objects have a gate mode but no gate-scale metadata. The
    # persisted numeric mode must continue to determine state evaluation.
    legacy <- old_default
    legacy$model$gate_scale <- NULL
    legacy$control$gate_scale <- NULL
    legacy$posterior$state$gate_scale <- NULL
    new_locations <- matrix(c(0.07, 0.47, 0.97), ncol = 1)
    expect_identical(ppt_predict(legacy, newdata = new_locations),
                     ppt_predict(explicit_root, newdata = new_locations))
    expect_identical(ppt_lppd(legacy, scale = 0.25),
                     ppt_lppd(explicit_root, scale = 0.25))
    expect_identical(
      ppt_marginal(legacy, variable = 1, grid = new_locations[, 1]),
      ppt_marginal(explicit_root, variable = 1, grid = new_locations[, 1])
    )
  }
})
