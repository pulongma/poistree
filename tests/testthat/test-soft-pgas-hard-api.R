test_that("soft PGAS records and validates the hard-proposal controls", {
  x <- matrix(c(0.1, 0.2, 0.4, 0.7, 0.8, 0.9), ncol = 1)
  args <- list(x = x, region = matrix(c(0, 1), 1), gating = "soft", sampler = "pgas",
               particles = 3, iter = 12, burn = 2, max_depth = 2,
               cut_candidates = 3, update_gate = FALSE, exact_max = 0,
               verbose = FALSE)
  default <- do.call(ppt_fit, args)
  expect_identical(default$control$proposal_score, "laplace")
  hard <- do.call(ppt_fit, c(args, list(proposal_score = "hard")))
  expect_identical(hard$control$proposal_score, "hard")
  expect_equal(hard$control$proposal_temperature, 0.5)
  expect_equal(hard$control$proposal_defensive, 0.1)
  expect_true(all(is.finite(hard$prediction$mean)))
  expect_length(hard$diagnostics$log_likelihood_trace, hard$posterior$draws)
  expect_true(all(is.finite(hard$diagnostics$log_likelihood_trace)))
  for (bad in list(0, -1, 1.1, NA_real_, Inf, c(.2, .4), "0.5")) {
    expect_error(do.call(ppt_fit, c(args, list(proposal_temperature = bad))),
                 "proposal_temperature")
  }
  for (bad in list(0, -1, 1, NA_real_, Inf, c(.2, .4), "0.1")) {
    expect_error(do.call(ppt_fit, c(args, list(proposal_defensive = bad))),
                 "proposal_defensive")
  }
  expect_error(do.call(ppt_fit, c(args, list(proposal_score = "invalid"))), "arg")
})

test_that("hard score controls leave the exact small-node algorithm unchanged", {
  x <- matrix(c(0.1, 0.2, 0.4, 0.7, 0.8, 0.9), ncol = 1)
  args <- list(x = x, region = matrix(c(0, 1), 1), gating = "soft", sampler = "pgas",
               particles = 4, iter = 40, burn = 5, max_depth = 2,
               cut_candidates = 3, update_gate = FALSE, exact_max = 150,
               seed = 425, verbose = FALSE)
  base <- do.call(ppt_fit, args)
  hard <- do.call(ppt_fit, c(args, list(proposal_score = "hard",
                 proposal_temperature = .2, proposal_defensive = .3)))
  expect_identical(base$posterior, hard$posterior)
  expect_identical(base$prediction, hard$prediction)
  expect_identical(base$diagnostics, hard$diagnostics)
})
