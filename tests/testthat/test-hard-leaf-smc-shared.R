test_that("shared-path SMC reproduces the exact one-step posterior at the root", {
  set.seed(3)
  n <- 60
  x <- cbind(c(runif(n / 2, 0, 0.4), runif(n / 2, 0.4, 1)), runif(n))
  region <- matrix(c(0, 1, 0, 1), 2, byrow = TRUE)
  key <- function(J, L) ifelse(is.na(J), "stop", sprintf("%d:%.6f", J, L))

  exact <- poistree:::PPT_transition_probabilities(x, region, 1L, 30L, Inf)
  p_exact <- stats::setNames(exact$probability, key(exact$axis - 1L, exact$cut))

  P <- 8000
  fit <- ppt_fit(x, region, gating = "hard", sampler = "smc", particles = P,
                 max_depth = 1, engine = "shared", cut_candidates = 30L, seed = 1)
  decision <- vapply(fit$posterior$tree_draws, function(nodes) {
    root <- nodes[[1]]
    if (isTRUE(root$is_leaf)) "stop" else key(root$J, root$L)
  }, character(1))
  freq <- as.numeric(table(factor(decision, levels = names(p_exact)))) / P
  z <- (freq - p_exact) / sqrt(p_exact * (1 - p_exact) / P)
  expect_true(all(decision %in% names(p_exact)))
  expect_lt(mean(z^2), 2)
  expect_lt(max(abs(z)), 4.5)
})

test_that("shared and dense SMC engines agree and share the API", {
  set.seed(5)
  x <- matrix(runif(300), ncol = 2)
  region <- matrix(c(0, 1, 0, 1), 2, byrow = TRUE)
  grid <- matrix(runif(20), ncol = 2)
  fits <- lapply(c("shared", "dense"), function(engine) {
    ppt_fit(x, region, gating = "hard", sampler = "smc", predict_at = grid,
            particles = 200, max_depth = 4, min_leaf_n = 2, engine = engine, seed = 7)
  })
  expect_identical(names(fits[[1]]$posterior), names(fits[[2]]$posterior))
  expect_identical(names(fits[[1]]$diagnostics), names(fits[[2]]$diagnostics))
  expect_identical(fits[[1]]$backend, "PPT_fit_SMC_shared")
  expect_identical(fits[[2]]$backend, "PPT_fit_SMC")
  expect_true(is.finite(fits[[1]]$diagnostics$expanded_nodes))
  expect_true(is.na(fits[[2]]$diagnostics$expanded_nodes))
  expect_true(all(is.finite(ppt_predict(fits[[1]]))))
  expect_length(fits[[1]]$diagnostics$ess_history, 2^4 - 1)

  expect_lt(abs(fits[[1]]$posterior$log_relative_normalizer -
                  fits[[2]]$posterior$log_relative_normalizer), 6)
  expect_true(is.na(fits[[1]]$posterior$log_evidence))
  expect_true(is.na(fits[[2]]$posterior$log_evidence))

  leaf_m <- unlist(lapply(fits[[1]]$posterior$tree_draws[[1]], function(node) if (isTRUE(node$is_leaf)) node$m))
  expect_equal(sum(leaf_m), nrow(x))
})
