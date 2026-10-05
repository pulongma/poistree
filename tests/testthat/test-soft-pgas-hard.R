hard_test_lse <- function(x) {
  top <- max(x)
  top + log(sum(exp(x - top)))
}

hard_test_axis_integral <- function(path, axis, gate, region, lower, upper,
                                    cut = NULL, side = NULL) {
  use <- path[path[, 1] == axis, , drop = FALSE]
  if (!is.null(cut)) use <- rbind(use, c(axis, cut, side))
  value <- function(x) {
    out <- rep(1, length(x))
    if (nrow(use)) for (k in seq_len(nrow(use))) {
      z <- gate[axis] * (x - use[k, 2]) / diff(region[axis, ])
      out <- out * plogis(use[k, 3] * z)
    }
    out
  }
  integrate(value, lower, upper, rel.tol = 1e-11, abs.tol = 1e-13,
            subdivisions = 1000)$value
}

test_that("hard histograms preserve cut order and strict boundary comparisons", {
  x <- cbind(c(-1, 0.5, 0.5, 1, 2.5), c(10.2, 11, 11, 12, 13.7))
  region <- rbind(c(-2, 3), c(10, 14))
  # Unsorted native grids and repeated cuts exercise the mapping back to the
  # original action order; observations exactly on cuts belong to the right.
  grids <- list(c(1, 0.5, 0.5, -0.5), c(12, 11, 13))
  path <- matrix(numeric(), 0, 3)
  out <- poistree:::SPPT_hard_node_probe(x, region, grids, c(9, 3.2), path)
  expected <- vapply(seq_along(out$cuts), function(c) {
    sum(x[, out$axis[c] + 1] < out$cuts[c])
  }, integer(1))
  expect_identical(out$count, expected)
  expect_equal(out$gate_evaluations_before, 0)
  expect_equal(out$gate_evaluations_after, 0)
  expect_equal(out$true_exposure_builds, 0)
  expect_equal(out$dense_cells, 0)
  expect_equal(out$bin_builds, 1)
  expect_true(all(is.nan(out$logHL)))
  expect_true(all(is.nan(out$logHR)))
})

test_that("hard exposure intervals preserve the original logistic scale", {
  x <- cbind(c(-1, 0.5, 0.5, 1, 2.5), c(10.2, 11, 11, 12, 13.7))
  region <- rbind(c(-2, 3), c(10, 14))
  grids <- list(c(1, 0.5, -0.5), c(12, 11, 13))
  paths <- list(matrix(numeric(), 0, 3),
                rbind(c(1, 0.4, -1), c(2, 11.8, 1)),
                rbind(c(1, 0.4, -1), c(1, 0.400000001, 1), c(2, 11.8, 1)))
  for (gate in list(c(0.03, 0.04), c(9, 3.2), c(80, 25))) for (path in paths) {
    out <- poistree:::SPPT_hard_node_probe(x, region, grids, gate, path, selected = 1L)
    parent_axis <- vapply(1:2, function(j) {
      hard_test_axis_integral(path, j, gate, region, region[j, 1], region[j, 2])
    }, numeric(1))
    expect_equal(exp(out$logH), prod(parent_axis), tolerance = 1e-9)
    for (c in seq_along(out$cuts)) {
      j <- out$axis[c] + 1
      left <- hard_test_axis_integral(path, j, gate, region, region[j, 1], out$cuts[c])
      right <- hard_test_axis_integral(path, j, gate, region, out$cuts[c], region[j, 2])
      other <- prod(parent_axis[-j])
      expect_equal(exp(out$proxy_logHL[c]), other * left, tolerance = 1e-8)
      expect_equal(exp(out$proxy_logHR[c]), other * right, tolerance = 1e-8)
    }
    chosen <- 2L
    j <- out$axis[chosen] + 1L
    true_left <- hard_test_axis_integral(path, j, gate, region,
      region[j, 1], region[j, 2], out$cuts[chosen], -1)
    true_right <- hard_test_axis_integral(path, j, gate, region,
      region[j, 1], region[j, 2], out$cuts[chosen], 1)
    expect_equal(exp(out$logHL[chosen]), prod(parent_axis[-j]) * true_left, tolerance = 1e-8)
    expect_equal(exp(out$logHR[chosen]), prod(parent_axis[-j]) * true_right, tolerance = 1e-8)
    expect_true(all(is.nan(out$logHL[-chosen])))
    expect_equal(out$true_exposure_builds, 1)
    expect_equal(out$gate_evaluations_before, 0)
    expect_equal(out$gate_evaluations_after, nrow(x))
    expect_equal(out$dense_cells, 0)
  }
})

test_that("tempered hard scores normalize with full prior-mixture support", {
  x <- matrix(c(0.02, 0.2, 0.5, 0.5, 0.8, 0.95), ncol = 1)
  region <- matrix(c(0, 1), 1)
  cuts <- seq(0.01, 0.99, length.out = 50)
  logQ <- function(k, H) 0.5 * log(0.15) - lgamma(0.5) + lgamma(0.5 + k) -
    (0.5 + k) * log(0.15 + H)
  for (temperature in c(0.1, 0.5, 1)) for (epsilon in c(0.001, 0.1, 0.8)) {
    out <- poistree:::SPPT_hard_node_probe(x, region, list(cuts), 0.001,
      matrix(numeric(), 0, 3), temperature, epsilon)
    gain <- logQ(out$count, exp(out$proxy_logHL)) +
      logQ(nrow(x) - out$count, exp(out$proxy_logHR)) - logQ(nrow(x), exp(out$logH))
    expected_score <- out$logprior + temperature * c(0, gain)
    expected_q <- (1 - epsilon) * exp(expected_score - hard_test_lse(expected_score)) +
      epsilon * exp(out$logprior)
    expect_equal(out$score, expected_score, tolerance = 1e-12)
    expect_equal(exp(out$logq), expected_q, tolerance = 1e-12)
    expect_equal(sum(exp(out$logq)), 1, tolerance = 1e-12)
    expect_true(all(exp(out$logq) >= epsilon * exp(out$logprior) - 1e-15))
    expect_length(out$count, 50L)
    expect_equal(out$gate_evaluations_after, 0)
  }
})

test_that("hard root cache respects gates, scorer changes, and new mixture controls", {
  gates <- rbind(c(4, 7), c(4, 7), c(5, 7), c(5, 7), c(5, 7), c(5, 7), c(5, 8))
  temperature <- c(0.5, 0.3, 0.3, 0.9, 0.5, 0.7, 0.7)
  epsilon <- c(0.1, 0.2, 0.2, 0.02, 0.1, 0.05, 0.05)
  hard <- c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE)
  cached <- poistree:::SPPT_hard_cache_probe(gates, temperature, epsilon, hard, TRUE)
  fresh <- poistree:::SPPT_hard_cache_probe(gates, temperature, epsilon, hard, FALSE)
  expect_identical(cached$output, fresh$output)
  expect_equal(cached$builds, c(2, 2, 3, 3, 5, 7, 8))
  expect_equal(cached$hits, c(0, 2, 3, 5, 5, 5, 6))
  expect_equal(cached$bins, rep(1, nrow(gates)))
  expect_equal(cached$true_builds, rep(0, nrow(gates)))
})

test_that("sparse gate entries invalidate per coordinate and remain numerically identical", {
  gates <- rbind(c(4, 7), c(4, 7), c(5, 7), c(5, 8))
  query <- rbind(c(1L, 1L), c(2L, 2L), c(3L, 4L), c(1L, 1L))
  requests <- rep(list(query), nrow(gates))
  lazy <- poistree:::SPPT_hard_gate_probe(gates, requests, TRUE)
  eager <- poistree:::SPPT_hard_gate_probe(gates, requests, FALSE)
  expect_identical(lazy$values, eager$values)
  expect_equal(lazy$evaluations, c(3, 3, 5, 6))
  expect_equal(eager$evaluations, c(30, 30, 48, 60))
  expect_equal(lazy$entries, rep(3, 4))
  # Tiny/disabled caches exercise eviction without stale last-entry reuse.
  for (capacity in c(0L, 1L, 2L)) {
    small <- poistree:::SPPT_hard_gate_probe(gates, requests, TRUE, capacity)
    expect_identical(small$values, eager$values)
    expect_true(all(small$entries <= 2 * capacity))
  }
  extreme <- matrix(c(2000, 1800), 1)
  all_points <- as.matrix(expand.grid(1:6, 1:5))
  expect_identical(
    poistree:::SPPT_hard_gate_probe(extreme, list(all_points), TRUE)$values,
    poistree:::SPPT_hard_gate_probe(extreme, list(all_points), FALSE)$values
  )
})

test_that("lazy and eager gate evaluation give identical complete PGAS trajectories", {
  for (hard in c(FALSE, TRUE)) for (exact_max in c(0L, 3L, 150L)) {
    for (config in list(c(FALSE, FALSE), c(TRUE, FALSE), c(TRUE, TRUE))) {
      set.seed(93051)
      lazy <- poistree:::SPPT_hard_gate_chain(hard, exact_max, config[1], config[2], TRUE)
      set.seed(93051)
      eager <- poistree:::SPPT_hard_gate_chain(hard, exact_max, config[1], config[2], FALSE)
      expect_identical(lazy, eager)
    }
  }
})

test_that("forced hard actions use true soft targets and sequential proposal weights", {
  logQ <- function(k, H) 0.5 * log(0.15) - lgamma(0.5) + lgamma(0.5 + k) -
    (0.5 + k) * log(0.15 + H)
  for (gate in list(c(0.03, 0.04), c(4, 7), c(800, 500))) {
    for (bits in list(rep(0L, 6), rep(1L, 6), c(0L, 1L, 0L, 1L, 0L, 1L))) {
      set.seed(93052)
      before <- .Random.seed
      out <- poistree:::SPPT_forced_hard_probe(bits, gate, 0.5, 0.1)
      expect_identical(.Random.seed, before)
      logits <- gate[1] * (0.2 - out$x[, 1])
      ll <- plogis(logits, log.p = TRUE)
      lr <- plogis(logits, lower.tail = FALSE, log.p = TRUE)
      H_left <- exp(out$logHL)
      H_right <- exp(out$logHR)
      # Independent analytic logistic integrals at a noncentral cut. Broad
      # gates have true exposures near 1/2, whereas the hard proxy uses 0.2/0.8.
      softplus <- function(z) pmax(z, 0) + log1p(exp(-abs(z)))
      expected_left <- (softplus(0.2 * gate[1]) - softplus(-0.8 * gate[1])) / gate[1]
      expected_right <- (softplus(0.8 * gate[1]) - softplus(-0.2 * gate[1])) / gate[1]
      expect_equal(H_left, expected_left, tolerance = 1e-12)
      expect_equal(H_right, expected_right, tolerance = 1e-12)
      expect_equal(exp(out$proxy_logHL), 0.2, tolerance = 1e-12)
      expect_equal(exp(out$proxy_logHR), 0.8, tolerance = 1e-12)
      if (gate[1] < 1) expect_gt(abs(H_left - exp(out$proxy_logHL)), 0.1)
      q_bits <- 0
      left <- right <- 0L
      for (i in seq_along(bits)) {
        w <- c(ll[i] + log(0.5 + left) - log(0.15 + H_left),
               lr[i] + log(0.5 + right) - log(0.15 + H_right))
        q_bits <- q_bits + w[if (bits[i]) 1 else 2] - hard_test_lse(w)
        left <- left + bits[i]
        right <- right + 1L - bits[i]
      }
      expected <- out$logprior[2] + sum(ifelse(bits == 1, ll, lr)) +
        logQ(sum(bits), H_left) + logQ(length(bits) - sum(bits), H_right) -
        out$logQA - out$logq[2] - q_bits
      expect_equal(out$logw, expected, tolerance = 1e-10)
      expect_equal(out$true_exposure_builds, 1)
      expect_equal(out$gate_evaluations, 6)
      expect_true(all(exp(out$logq) >= 0.1 * exp(out$logprior) - 1e-15))
    }
  }
})
