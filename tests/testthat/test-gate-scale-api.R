# gate_structure controls parameter sharing; gate_scale controls the width.
gate_scale_test_args <- function() {
  list(
    x = cbind(c(.08, .21, .34, .62, .79, .93),
              c(.12, .71, .25, .86, .43, .61)),
    region = matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE),
    gating = "soft", sampler = "pcg", a = .5, b = .2,
    max_depth = 2L, cut_candidates = 3L, chains = 1L,
    iter = 18L, burn = 8L, thin = 2L,
    tree_moves = 1L, change_moves = 1L, seed = 112L, verbose = FALSE
  )
}

gate_scale_test_fit <- function(...) {
  do.call(ppt_fit, utils::modifyList(gate_scale_test_args(), list(...)))
}

test_that("gate scales validate independently of gate sharing", {
  resolve <- poistree:::.ppt_resolve_gate_scale
  expect_identical(resolve(NULL, "logistic"), "root")
  expect_identical(resolve(NULL, "compact"), "node")
  expect_identical(resolve("node", "logistic"), "node")
  expect_identical(resolve("root", "logistic"), "root")
  expect_identical(resolve("node", "compact"), "node")
  expect_error(resolve("root", "compact"), "Compact gates require")
  for (invalid in list(NA_character_, NA, 1, character(), c("root", "node"),
                       "dimension", "local", "n")) {
    expect_error(resolve(invalid), "gate_scale")
  }
  expect_error(gate_scale_test_fit(gate_structure = "node"), "arg")
  expect_error(gate_scale_test_fit(gate_family = "compact", gate_scale = "root"),
               "Compact gates require")
})

test_that("omitting gate_scale preserves historical logistic and compact fits", {
  for (sampler in c("rjmcmc", "irjmcmc", "pcg")) {
    for (family in c("logistic", "compact")) {
      old_default <- gate_scale_test_fit(sampler = sampler, gate_family = family)
      scale <- if (family == "logistic") "root" else "node"
      explicit <- gate_scale_test_fit(sampler = sampler, gate_family = family,
                                      gate_scale = scale)
      expect_identical(old_default$prediction, explicit$prediction)
      expect_identical(old_default$posterior, explicit$posterior)
      expect_identical(old_default$diagnostics, explicit$diagnostics)
      expect_identical(old_default$model$gate_scale, scale)
      expect_identical(old_default$prior$gate$scale, scale)
      expect_identical(old_default$control$gate_scale, scale)
      # Serialized fits made before gate_scale was introduced have no new
      # metadata fields. Their stored gate_mode remains authoritative.
      legacy <- old_default
      legacy$model$gate_scale <- legacy$prior$gate$scale <-
        legacy$control$gate_scale <- NULL
      at <- matrix(c(.15, .45, .55, .85), ncol = 2L)
      expect_identical(ppt_predict(old_default, newdata = at),
                       ppt_predict(legacy, newdata = at))
    }
  }
})

test_that("all soft MCMC samplers retain node-scale state for both structures", {
  for (sampler in c("rjmcmc", "irjmcmc", "pcg")) {
    for (structure in c("dimension", "shared")) {
      fit <- gate_scale_test_fit(sampler = sampler, gate_scale = "node",
                                 gate_structure = structure)
      expect_identical(fit$model$gate_scale, "node")
      expect_identical(fit$prior$gate$scale, "node")
      expect_identical(fit$control$gate_scale, "node")
      expect_identical(fit$prior$gate$structure, structure)
      expect_identical(fit$posterior$state$gate_mode, 3L)
      expect_true(all(is.finite(fit$prediction$draws)))
      expect_true(all(fit$prediction$draws > 0))
      # Call the state evaluator directly so this is not the precomputed path.
      evaluated <- poistree:::.ppt_state_eval(fit, fit$prediction$locations)
      expect_equal(evaluated, fit$prediction$draws, tolerance = 1e-10)
      restored <- unserialize(serialize(fit, NULL))
      at <- matrix(c(.15, .45, .55, .85), ncol = 2L)
      expect_identical(ppt_predict(fit, newdata = at),
                       ppt_predict(restored, newdata = at))
      if (structure == "shared") {
        expect_equal(fit$posterior$state$gate[, 1],
                     fit$posterior$state$gate[, 2])
      }
    }
  }
})

# A repeated split on x has parent width .8 rather than the root width 2.
# No inference kernel is involved, so the check does not depend on visited trees.
node_scale_known_fit <- function(mode = 3L) {
  nodes <- cbind(
    heap_id = 1:7, axis = c(0, 0, 1, -1, -1, -1, -1),
    cut = c(.8, .3, .2, NA, NA, NA, NA),
    lambda = c(0, 0, 0, 2, 3, 5, 7), xi = NA_real_, m = 0
  )
  x <- matrix(c(.2, -.8, 1.5, 1.2), ncol = 2L, byrow = TRUE,
               dimnames = list(NULL, c("x", "y")))
  structure(list(
    model = list(gating = "soft", scales = "leaf", gate_family = "logistic"),
    data = list(x = x, dimension = 2L,
                region = matrix(c(0, 2, -1, 2), 2L, byrow = TRUE)),
    prediction = list(locations = matrix(numeric(), 0L, 2L),
                      draws = matrix(numeric(), 0L, 1L)),
    posterior = list(tree_draws = list(), particle_weights = numeric(),
                     state = list(mode = "heap", nodes = list(nodes),
                                  gate = matrix(c(6, 4), 1L),
                                  gate_mode = mode, gate_depth = 0))
  ), class = "ppt")
}

node_scale_known_intensity <- function(x, y, local_width = .8) {
  root <- plogis(6 * (x - .8) / 2)
  left <- plogis(6 * (x - .3) / local_width)
  right <- plogis(4 * (y - .2) / 3)
  (1 - root) * (2 * (1 - left) + 3 * left) +
    root * (5 * (1 - right) + 7 * right)
}

test_that("serialized node gates use the splitting parent's local width", {
  fit <- node_scale_known_fit()
  at <- cbind(c(.1, .4, .7, 1.7), c(-.8, 0, .9, 1.8))
  expected <- node_scale_known_intensity(at[, 1], at[, 2])
  observed <- poistree:::.ppt_state_eval(fit, at)
  expect_equal(as.numeric(observed), expected, tolerance = 1e-12)
  root_fit <- node_scale_known_fit(2L)
  root_expected <- node_scale_known_intensity(at[, 1], at[, 2], local_width = 2)
  expect_equal(as.numeric(poistree:::.ppt_state_eval(root_fit, at)),
               root_expected, tolerance = 1e-12)
  expect_gt(max(abs(expected - root_expected)), .01)
})

test_that("node-scale marginal integrals agree with independent R quadrature", {
  fit <- node_scale_known_fit()
  grid <- c(-1, -.3, .2, 1.1, 2)
  expected <- vapply(grid, function(y) {
    integrate(function(x) node_scale_known_intensity(x, y), 0, 2,
              rel.tol = 1e-11, abs.tol = 1e-12)$value / 2
  }, numeric(1))
  result <- ppt_marginal(fit, "y", grid = grid, type = "draws")
  expect_equal(as.numeric(result), expected, tolerance = 2e-9)
  expect_identical(attr(result, "method"), "adaptive posterior-state integration")
  projected <- ppt_marginal(fit, "y", grid = grid, average = FALSE,
                            type = "draws")
  expect_equal(as.numeric(projected), 2 * expected, tolerance = 2e-9)
})
