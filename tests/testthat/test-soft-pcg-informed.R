informed_pcg_args <- function() {
  set.seed(7)
  x <- cbind(c(runif(28, 0, .5), runif(12, .5, 1)), runif(40))
  list(x = x, region = matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE),
       sampler = "pcg", a = .5, b = .5 / 40, gate = 8, max_depth = 2L,
       cut_candidates = 3L, chains = 1L, iter = 30L, burn = 10L, thin = 2L,
       tree_moves = 2L, change_moves = 2L, seed = 31L, verbose = FALSE,
       predict_at = x[1:2, ])
}

informed_pcg_fit <- function(...) {
  do.call(ppt_fit, utils::modifyList(informed_pcg_args(), list(...)))
}

test_that("informed PCG validates its proposal controls", {
  expect_error(informed_pcg_fit(informed = NA), "`informed`")
  expect_error(informed_pcg_fit(informed = TRUE, proposal_temperature = 0),
               "proposal_temperature")
  expect_error(informed_pcg_fit(informed = TRUE, proposal_temperature = 1.5),
               "proposal_temperature")
  expect_error(informed_pcg_fit(informed = TRUE, proposal_defensive = 1),
               "proposal_defensive")
  expect_error(informed_pcg_fit(informed = TRUE, proposal_defensive = -.1),
               "proposal_defensive")
})

test_that("informed PCG records its controls and runs for every gate family", {
  for (family in c("logistic", "compact")) {
    for (structure in c("dimension", "shared")) {
      fit <- informed_pcg_fit(informed = TRUE, gate_family = family,
                              gate_structure = structure)
      expect_s3_class(fit, "ppt")
      expect_identical(fit$model$algorithm,
                       "Partially collapsed Gibbs (RAM, informed tree proposals)")
      expect_true(fit$control$informed)
      expect_identical(fit$control$proposal_temperature, 0.5)
      expect_identical(fit$control$proposal_defensive, 0.1)
      expect_true(all(is.finite(fit$prediction$draws)))
    }
  }
  plain <- informed_pcg_fit()
  expect_false(plain$control$informed)
  expect_null(plain$control$proposal_temperature)
})

test_that("the standard PCG path ignores the proposal controls", {
  a <- informed_pcg_fit()
  b <- informed_pcg_fit(proposal_temperature = .9, proposal_defensive = .5)
  a$call <- b$call <- NULL
  expect_identical(a, b)
  informed <- informed_pcg_fit(informed = TRUE)
  expect_false(identical(informed$prediction$draws, a$prediction$draws))
})

test_that("informed PCG proposals leave the tree-size posterior unchanged", {
  skip_on_cran()
  size_distribution <- function(...) {
    fit <- informed_pcg_fit(iter = 40000L, burn = 4000L, thin = 4L,
                            update_gate = FALSE, ...)
    leaves <- fit$diagnostics$leaf_count_trace
    vapply(1:4, function(k) mean(leaves == k), numeric(1))
  }
  standard <- size_distribution(seed = 1L)
  informed <- size_distribution(seed = 2L, informed = TRUE)

  expect_true(all(abs(standard - informed) < 0.04))
  expect_equal(sum(informed), 1)
})
