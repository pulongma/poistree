test_that("shared ancestor prediction agrees with independent leaf probabilities", {
  locations <- rbind(as.matrix(expand.grid(
    x = c(-2, -1, -0.3, 0, 1.8, 3), y = c(4, 6, 7, 12))),
    c(-0.30000001, 7.00000001))
  tree <- prediction_tree()

  nodes <- list(tree, tree[c(9, 3, 6, 1, 8, 4, 7, 2, 5), ])
  gates <- rbind(c(8, 13), c(1e6, 1e3))
  for (mode in c(1L, 2L, 3L)) {
    depth <- if (mode == 1L) 1.25 else 0
    expected <- prediction_reference(nodes, gates, locations,
                                     prediction_region, mode, depth)
    actual <- poistree:::ppt_eval_state(nodes, gates, prediction_region,
                                        locations, mode, depth)
    expect_equal(actual, expected, tolerance = 2e-12, ignore_attr = TRUE)
    expect_true(all(is.finite(actual)))

    equal_tree <- tree
    equal_tree[equal_tree[, 2] < 0, 4] <- 7
    constant <- poistree:::ppt_eval_state(list(equal_tree), gates[1, , drop = FALSE],
                                          prediction_region, locations, mode, depth)
    expect_equal(as.numeric(constant), rep(7, nrow(locations)), tolerance = 2e-12)
  }
})

test_that("root-only, invalid-rate and deep skew states retain legacy semantics", {
  points <- rbind(c(-2, 4), c(-2 + 1e-6, 8), c(-1.999, 7), c(0, 9), c(3, 12))
  gates <- matrix(c(120, 9), 1)
  tree <- prediction_tree()
  tree[tree[, 2] < 0, 4] <- c(0, -2, Inf, NA_real_, 6)
  root <- matrix(c(1, -1, NA_real_, 4.25), 1)
  skew <- matrix(numeric(), 0, 4)
  width <- diff(prediction_region[1, ])
  for (depth in 0:24) {
    id <- 2^depth
    cut <- prediction_region[1, 1] + width / 2
    skew <- rbind(skew, c(id, 0, cut, 0),
                  c(2 * id + 1, -1, NA_real_, 1 + depth / 7))
    width <- width / 2
  }
  skew <- rbind(skew, c(2^25, -1, NA_real_, 3))
  for (mode in c(1L, 2L, 3L)) {
    for (state in list(tree, root, skew[nrow(skew):1, ])) {
      expected <- prediction_reference(list(state), gates, points,
                                       prediction_region, mode, 0.75)
      actual <- poistree:::ppt_eval_state(list(state), gates, prediction_region,
                                          points, mode, 0.75)
      expect_equal(actual, expected, tolerance = 3e-12, ignore_attr = TRUE)
    }
    empty <- poistree:::ppt_eval_state(list(root), gates, prediction_region,
                                       matrix(numeric(), 0, 2), mode, 0)
    expect_identical(dim(empty), c(0L, 1L))
  }
})

test_that("prediction row matching distinguishes adjacent doubles and accepts signed zero", {

  one_up <- 1 + .Machine$double.eps
  one_down <- 1 - .Machine$double.eps / 2
  locations <- rbind(c(0, 5), c(1, 5), c(one_up, 5), c(one_down, 5), c(1, 5))
  query <- rbind(c(-0, 5), c(one_up, 5), c(1, 5), c(one_down, 5),
                 c(1 + 2 * .Machine$double.eps, 5), c(0, 5))
  expect_identical(poistree:::.ppt_match_prediction_rows(query, locations),
                   c(1L, 3L, 2L, 4L, NA_integer_, 1L))
  expect_identical(poistree:::.ppt_match_prediction_rows(matrix(numeric(), 0, 2), locations),
                   integer())
  expect_identical(poistree:::.ppt_match_prediction_rows(query, matrix(numeric(), 0, 2)),
                   rep(NA_integer_, nrow(query)))
})

test_that("stored posterior draws are reused for full, repeated and mixed locations", {
  fit <- prediction_fixture()
  locations <- fit$prediction$locations
  fresh <- rbind(c(-0.6, 10), c(2.5, 5))

  fit$prediction$draws <- fit$prediction$draws + 100
  for (index in list(1:4, c(4L, 2L, 2L, 1L))) {
    raw <- ppt_lambda(fit, at = locations[index, , drop = FALSE], type = "draws")
    expect_equal(raw$draws, fit$prediction$draws[index, , drop = FALSE], ignore_attr = TRUE)
    expect_identical(dim(raw$draws), c(length(index), 4L))
    expect_equal(raw$weights, c(0, 0.2, 0.5, 0.3))
  }
  mixed <- rbind(locations[3, ], fresh[1, ], locations[3, ], fresh[2, ])
  fresh_draws <- prediction_reference(fit$posterior$state$nodes,
                                    fit$posterior$state$gate, fresh,
                                    prediction_region, 3L)
  expected <- rbind(fit$prediction$draws[3, ], fresh_draws[1, ],
                     fit$prediction$draws[3, ], fresh_draws[2, ])
  expect_equal(ppt_lambda(fit, at = mixed, type = "draws")$draws,
               expected, tolerance = 1e-12, ignore_attr = TRUE)

  direct <- poistree:::.ppt_state_eval(fit, locations)
  expect_equal(direct, fit$prediction$draws - 100,
               tolerance = 1e-12, ignore_attr = TRUE)
  fit$prediction$draws <- NULL
  expect_equal(ppt_lambda(fit, at = fresh, type = "draws")$draws,
               fresh_draws, tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(ppt_lambda(fit, at = locations, type = "draws")$draws,
               direct, tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("only missing rows require native posterior evaluation", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat"), inherits = FALSE))
  fit <- prediction_fixture()
  fresh <- matrix(c(2.6, 10), 1)
  evaluated <- list()
  local_mocked_bindings(.ppt_state_eval = function(object, at) {
    evaluated[[length(evaluated) + 1L]] <<- at
    prediction_reference(object$posterior$state$nodes,
                         object$posterior$state$gate, at,
                         object$data$region, object$posterior$state$gate_mode)
  }, .package = "poistree")
  ppt_lambda(fit, at = fit$prediction$locations, type = "draws")
  expect_length(evaluated, 0)
  mixed <- rbind(fit$prediction$locations[2, ], fresh,
                 fit$prediction$locations[1, ], fresh)
  ppt_lambda(fit, at = mixed, type = "draws")
  expect_length(evaluated, 1)
  expect_true(all(evaluated[[1]][, 1] == fresh[1, 1]))
  expect_true(all(evaluated[[1]][, 2] == fresh[1, 2]))

  expect_true(nrow(evaluated[[1]]) %in% 1:2)
})

test_that("weighted pointwise summaries and public prediction conventions are preserved", {
  fit <- prediction_fixture()
  at <- rbind(fit$prediction$locations[2, ], c(1.1, 9), fit$prediction$locations[2, ])
  base <- prediction_reference(list(prediction_tree()), matrix(c(8, 13), 1),
                               at, prediction_region, 3L)[, 1]
  summary <- ppt_lambda(fit, at = at, level = 0.95)
  expect_s3_class(summary, "ppt_lambda")
  expect_named(summary, c("longitude", "season", "mean", "median", "lower", "upper"))
  expect_equal(summary$mean, 1.23 * base, tolerance = 1e-12)
  expect_equal(summary$median, base, tolerance = 1e-12)
  expect_equal(summary$lower, 0.8 * base, tolerance = 1e-12)
  expect_equal(summary$upper, 1.9 * base, tolerance = 1e-12)
  expect_equal(attr(summary, "level"), 0.95)
  expect_equal(as.numeric(ppt_predict(fit, newdata = at, type = "mean")), summary$mean)
  expect_equal(as.numeric(ppt_predict(fit, newdata = at)), summary$median)
  interval <- ppt_predict(fit, newdata = at, type = "interval")
  expect_named(interval, c("median", "lower", "upper"))
  expect_equal(interval$lower, summary$lower)
  expect_equal(interval$upper, summary$upper)
  expect_equal(as.numeric(ppt_predict(fit)), fit$prediction$median)
  expect_equal(ppt_predict(fit, type = "interval")$lower, fit$prediction$lower)
  empty <- matrix(numeric(), 0, 2)
  expect_length(ppt_predict(fit, newdata = empty), 0)
  expect_length(ppt_predict(fit, newdata = empty, type = "mean"), 0)
  expect_equal(nrow(ppt_predict(fit, newdata = empty, type = "interval")), 0)
  expect_error(ppt_lambda(fit, at = empty), "empty|observation|row", ignore.case = TRUE)
  expect_error(poistree:::.ppt_state_eval(fit, empty), "empty|observation|row", ignore.case = TRUE)
})

test_that("posterior and plug-in lppd retain the exact scaled Poisson likelihood", {
  fit <- prediction_fixture(root_only = TRUE)
  test <- rbind(fit$prediction$locations[c(2, 1, 2), ], c(-1.1, 10), c(2.7, 6))
  rates <- c(2, 2, 5, 9)
  weights <- c(0, 0.2, 0.5, 0.3)
  for (scale in c(1, 0.25, 1.7)) {
    lp <- nrow(test) * log(scale * rates) - scale * 40 * rates
    expected_posterior <- max(lp) + log(sum(weights * exp(lp - max(lp))))
    expected_plugin <- nrow(test) * log(scale * sum(weights * rates)) -
      scale * 40 * sum(weights * rates)
    posterior <- ppt_lppd(fit, test = test, scale = scale)
    plugin <- ppt_lppd(fit, test = test, scale = scale, type = "plugin")
    expect_equal(as.numeric(posterior), expected_posterior, tolerance = 1e-12)
    expect_equal(as.numeric(plugin), expected_plugin, tolerance = 1e-12)
    expect_identical(attr(posterior, "n_test"), nrow(test))
    expect_true(attr(posterior, "joint"))
    expect_false(attr(plugin, "joint"))
    expect_equal(attr(posterior, "scale"), scale)
  }
})

test_that("mean prediction avoids sorting posterior credible intervals", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat"), inherits = FALSE))
  fit <- prediction_fixture()
  fresh <- matrix(c(-0.7, 10), 1)
  expected <- prediction_reference(fit$posterior$state$nodes,
                                   fit$posterior$state$gate, fresh,
                                   prediction_region, 3L) %*% c(0, 0.2, 0.5, 0.3)
  local_mocked_bindings(.ppt_weighted_quantiles = function(...) {
    stop("Mean-only prediction must not compute quantiles.")
  }, .package = "poistree")
  expect_equal(as.numeric(ppt_predict(fit, newdata = fresh, type = "mean")),
               as.numeric(expected), tolerance = 1e-12)
})

test_that("one-location and one-draw posterior predictions keep their shapes", {
  fit <- prediction_fixture(root_only = TRUE)
  fit$posterior$state$nodes <- fit$posterior$state$nodes[1]
  fit$posterior$state$gate <- fit$posterior$state$gate[1, , drop = FALSE]
  fit$posterior$particle_weights <- numeric()
  fit$prediction$draws <- fit$prediction$draws[, 1, drop = FALSE]
  at <- fit$prediction$locations[2, , drop = FALSE]
  raw <- ppt_lambda(fit, at = at, type = "draws")
  expect_identical(dim(raw$draws), c(1L, 1L))
  expect_identical(raw$weights, 1)
  summary <- ppt_lambda(fit, at = at)
  expect_identical(dim(summary), c(1L, 6L))
  expect_equal(unname(unlist(summary[3:6])), rep(2, 4))
})

test_that("inconsistent stored-draw counts fall back to complete posterior states", {
  fit <- prediction_fixture()

  fit$prediction$draws <- fit$prediction$draws[, 1:2, drop = FALSE]
  for (at in list(fit$prediction$locations,
                 rbind(fit$prediction$locations[2, ], c(1.1, 9)))) {
    expected <- prediction_reference(fit$posterior$state$nodes,
                                     fit$posterior$state$gate, at,
                                     prediction_region, 3L)
    raw <- ppt_lambda(fit, at = at, type = "draws")
    expect_identical(dim(raw$draws), c(nrow(at), 4L))
    expect_equal(raw$draws, expected, tolerance = 1e-12, ignore_attr = TRUE)
    expect_equal(raw$weights, c(0, 0.2, 0.5, 0.3))
  }
})
