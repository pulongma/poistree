# Independent R geometry: each row is (axis, cut, side, parent width).
node_gate_test_phi <- function(points, path, gate) {
  value <- rep(1, nrow(points))
  for (k in seq_len(nrow(path))) {
    j <- as.integer(path[k, 1L])
    value <- value * plogis(path[k, 3L] * gate[j] *
                            (points[, j] - path[k, 2L]) / path[k, 4L])
  }
  value
}

node_gate_test_exposure <- function(path, region, gate) {
  prod(vapply(seq_len(nrow(region)), function(j) {
    selected <- path[path[, 1L] == j, , drop = FALSE]
    integrate(function(x) {
      value <- rep(1, length(x))
      for (k in seq_len(nrow(selected))) {
        value <- value * plogis(selected[k, 3L] * gate[j] *
                                (x - selected[k, 2L]) / selected[k, 4L])
      }
      value
    }, region[j, 1L], region[j, 2L], rel.tol = 1e-11,
    abs.tol = 1e-13, subdivisions = 2000L)$value
  }, numeric(1)))
}

test_that("node logistic paths use local widths and full-domain integration", {
  region <- rbind(c(-2, 4), c(10, 18))
  points <- cbind(seq(-2, 4, length.out = 31),
                  seq(10, 18, length.out = 31))
  points <- rbind(points, c(1.82, 14))
  path <- rbind(c(1, 1, 1, 6), c(2, 15, -1, 8),
                c(1, 2.5, -1, 3), c(1, 1.8, 1, 1.5))
  for (gate in list(c(3, 8), c(20, 2), c(500, 35))) {
    args <- list(as.integer(path[, 1]), path[, 2], path[, 4],
                 as.integer(path[, 3]), points, region, gate)
    local <- do.call(poistree:::ppstree_geometry,
                     c(args, list(gate_scale = 1L)))
    expect_equal(as.vector(local$phi), node_gate_test_phi(points, path, gate),
                 tolerance = 1e-12)
    expect_equal(local$H, node_gate_test_exposure(path, region, gate),
                 tolerance = 2e-9)
    root <- do.call(poistree:::ppstree_geometry, args)
    explicit_root <- do.call(poistree:::ppstree_geometry,
                             c(args, list(gate_scale = 0L)))
    expect_identical(root, explicit_root)
    expect_gt(max(abs(root$phi - local$phi)), 1e-7)
  }
})

test_that("node logistic integration resolves very thin local transitions", {
  # The final two gates create a narrow leaf that whole-domain quadrature
  # can miss. Use independent Gaussian quadrature on small physical panels.
  path <- rbind(c(1, .5, 1, 1), c(1, .500001, -1, .5),
                c(1, .5000005, 1, 1e-6), c(1, .5000008, -1, 5e-7))
  region <- matrix(c(0, 1), nrow = 1)
  gate <- 500
  breaks <- sort(unique(c(0, 1, path[, 2],
                           path[, 2] - 40 * path[, 4] / gate,
                           path[, 2] + 40 * path[, 4] / gate)))
  breaks <- breaks[breaks >= 0 & breaks <= 1]
  reference <- sum(vapply(seq_len(length(breaks) - 1L), function(k) {
    integrate(function(x) node_gate_test_phi(matrix(x, ncol = 1), path, gate),
              breaks[k], breaks[k + 1L], rel.tol = 1e-9,
              abs.tol = 1e-20, subdivisions = 1000L)$value
  }, numeric(1)))
  result <- poistree:::ppstree_geometry(
    as.integer(path[, 1]), path[, 2], path[, 4], as.integer(path[, 3]),
    matrix(c(.5000006, .5000007), ncol = 1), region, gate, gate_scale = 1L
  )
  expect_gt(reference, 1e-8)
  expect_lt(abs(result$H - reference) / reference, 1e-7)
  # Equal-width, very steep gates must also resolve the transition.
  sharp <- poistree:::ppstree_geometry(
    1L, .2, 1, 1L, matrix(.2, ncol = 1), region, 1e10, gate_scale = 1L
  )
  expect_lt(abs(sharp$H - .8), 1e-11)
})

test_that("node logistic PCG targets use the same basis and exposure", {
  region <- rbind(c(-2, 4), c(10, 18))
  x <- cbind(c(-1.8, -0.3, 1.5, 3.8), c(10.4, 17.2, 12.5, 16))
  splits <- rbind(c(1, 0, 1), c(2, 0, -0.5), c(3, 1, 15))
  paths <- list(
    rbind(c(1, 1, -1, 6), c(1, -0.5, -1, 3)),
    rbind(c(1, 1, -1, 6), c(1, -0.5, 1, 3)),
    rbind(c(1, 1, 1, 6), c(2, 15, -1, 8)),
    rbind(c(1, 1, 1, 6), c(2, 15, 1, 8))
  )
  lambda <- c(.3, .7, 1.1, 1.8)
  for (shared in c(FALSE, TRUE)) {
    gate <- if (shared) c(6, 6) else c(5, 9)
    shape <- c(2, 2)
    rate <- c(.2, .2)
    phi <- vapply(paths, function(path) node_gate_test_phi(x, path, gate),
                  numeric(nrow(x)))
    exposure <- vapply(paths, node_gate_test_exposure, numeric(1),
                       region = region, gate = gate)
    result <- poistree:::ppstree_pcg_inspect(
      x, region, splits, gate, lambda, shape, rate, c(0, 0),
      as.integer(shared), 2L
    )
    p <- if (shared) 1L else seq_along(gate)
    prior <- sum((shape[p] - 1) * log(gate[p]) - rate[p] * gate[p])
    weighted <- sweep(phi, 2L, lambda, "*")
    target <- prior - sum(lambda * exposure) + sum(log(rowSums(weighted)))
    expect_equal(as.integer(result$leaf_ids), 4:7)
    expect_equal(result$phi, phi, tolerance = 1e-12)
    expect_equal(as.vector(result$exposure), exposure, tolerance = 2e-9)
    expect_equal(rowSums(phi), rep(1, nrow(x)), tolerance = 1e-12)
    expect_equal(sum(exposure), prod(region[, 2] - region[, 1]),
                 tolerance = 2e-9)
    expect_equal(result$log_target, target, tolerance = 2e-9)
    expect_equal(result$log_scale_target, target + sum(log(gate[p])),
                 tolerance = 2e-9)
    expect_equal(result$allocation_prob, weighted / rowSums(weighted),
                 tolerance = 1e-12)
  }
})

test_that("node logistic informed proposals obey the augmented target", {
  x <- matrix(c(.1, .25, .4, .6, .75, .9), ncol = 1)
  region <- matrix(c(0, 1), nrow = 1)
  splits <- rbind(c(1, 0, .6), c(2, 0, .25), c(3, 0, .75))
  labels <- c(4L, 5L, 6L, 7L, 5L, 6L)
  # Reuse only the independent cut-support and tree-prior reference. Recompute
  # every basis and integral below using the new node-relative R formula.
  reference <- informed_exact(x, soft = TRUE, ncand = 4L, geometry_only = TRUE)
  decode <- function(splits) {
    tree <- list()
    for (k in seq_len(nrow(splits))) {
      tree[[as.character(splits[k, 1])]] <- c(splits[k, 2] + 1, splits[k, 3])
    }
    tree
  }
  key <- function(splits, labels) reference$key(decode(splits), as.integer(labels))
  target <- function(splits, labels) {
    info <- reference$detail(decode(splits))
    paths <- lapply(info$leaves, `[[`, "path")
    phi <- vapply(paths, function(path) node_gate_test_phi(x, path, 5),
                  numeric(nrow(x)))
    exposure <- vapply(paths, node_gate_test_exposure, numeric(1),
                       region = region, gate = 5)
    index <- match(labels, names(info$leaves))
    counts <- tabulate(index, length(paths))
    info$lp + sum(lgamma(.5 + counts) - lgamma(.5) + .5 * log(.1) -
                    (.5 + counts) * log(.1 + exposure)) +
      sum(log(phi[cbind(seq_len(nrow(x)), index)]))
  }
  law <- function(splits, labels, kind) poistree:::ppstree_informed_transition(
    x, region, splits, as.integer(labels), 5, a = .5, b = .1,
    alpha = .65, eta = 1, Dmax = 2L, nmin = 1L, cut_mode = 1L,
    ncand = 4L, gate_family = 2L, kind = kind
  )
  original_key <- key(splits, labels)
  original_target <- target(splits, labels)
  for (kind in 0:1) {
    current <- law(splits, labels, kind)
    expect_equal(sum(vapply(current$neighbors, `[[`, numeric(1), "probability")),
                 1, tolerance = 1e-12)
    nonself <- Filter(function(s) key(s$splits, s$labels) != original_key,
                      current$neighbors)
    for (next_state in head(nonself, 8L)) {
      reverse <- law(next_state$splits, next_state$labels, kind)
      back <- Filter(function(s) key(s$splits, s$labels) == original_key,
                     reverse$neighbors)
      expect_length(back, 1L)
      difference <- target(next_state$splits, next_state$labels) - original_target
      expect_equal(next_state$log_ratio,
                   difference + back[[1]]$log_q0 - next_state$log_q0,
                   tolerance = 2e-8)
      forward <- next_state$probability *
        min(1, exp(current$log_normalizer - reverse$log_normalizer))
      backward <- back[[1]]$probability *
        min(1, exp(reverse$log_normalizer - current$log_normalizer))
      expect_equal(forward, exp(difference) * backward, tolerance = 2e-8)
    }
  }
})

test_that("node logistic MCMC fits reproduce retained predictions and integrals", {
  x <- matrix(c(-1.8, -1.7, -1.5, -.8, .1, .4, 2.5, 3.5, 3.6, 3.8), ncol = 1)
  region <- matrix(c(-2, 4), nrow = 1)
  for (sampler in c("rjmcmc", "irjmcmc", "pcg")) {
    for (structure in c("shared", "dimension")) {
      fit <- ppt_fit(
        x, region, gating = "soft", sampler = sampler,
        gate_family = "logistic", gate_scale = "node",
        gate_structure = structure, gate = 5, a = .5, b = .2,
        alpha = .8, eta = .1, max_depth = 3L, cut_candidates = 4L,
        chains = 1L, iter = 50L, burn = 20L, thin = 5L,
        tree_moves = 2L, change_moves = 1L, seed = 917L, verbose = FALSE
      )
      expect_s3_class(fit, "ppt")
      evaluated <- ppt_lambda(fit, at = x, type = "draws")
      expect_equal(evaluated$draws, fit$prediction$draws, tolerance = 2e-10)
      expect_true(all(is.finite(evaluated$draws)))
      expect_true(all(evaluated$draws > 0))
      native_integrals <- as.vector(ppt_integral(fit, "draws"))
      for (s in seq_len(min(3L, fit$posterior$draws))) {
        numerical <- integrate(function(t) {
          ppt_lambda(fit, at = matrix(t, ncol = 1), type = "draws")$draws[, s]
        }, -2, 4, rel.tol = 1e-9, subdivisions = 1000L)$value
        expect_equal(native_integrals[s], numerical, tolerance = 2e-8)
      }
    }
  }
})
