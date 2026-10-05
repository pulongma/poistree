# Independent enumeration oracle: sum over all 2^m allocations, without a
# Poisson-binomial recurrence or an exponential-tilting identity.
cut_tilt_lse <- function(x) {
  if (!length(x) || all(x == -Inf)) return(-Inf)
  anchor <- max(x)
  anchor + log(sum(exp(x - anchor)))
}

cut_tilt_enumerate <- function(x, cut, gate, width, H_left, H_right, a, b) {
  m <- length(x)
  allocations <- if (m == 0L) matrix(numeric(), 1L, 0L) else
    as.matrix(expand.grid(rep(list(0:1), m)))
  counts <- rowSums(allocations)
  logits <- gate * (cut - x) / width
  log_left <- plogis(logits, log.p = TRUE)
  log_right <- plogis(logits, lower.tail = FALSE, log.p = TRUE)
  log_weights <- vapply(seq_len(nrow(allocations)), function(i) {
    sum(ifelse(allocations[i, ] == 1, log_left, log_right))
  }, numeric(1))
  logQ <- function(n, H) {
    a * log(b) - lgamma(a) + lgamma(a + n) - (a + n) * log(b + H)
  }
  list(
    allocations = allocations,
    log_joint = log_weights,
    log_pmf = vapply(0:m, function(k) cut_tilt_lse(log_weights[counts == k]), numeric(1)),
    log_score = cut_tilt_lse(log_weights + logQ(counts, H_left) + logQ(m - counts, H_right))
  )
}

test_that("exact cut tilting agrees with independently enumerated allocations", {
  set.seed(3107)
  cuts <- c(1.8, -1.1, 0.35, 2.6, -0.2)
  # Deliberately unequal exposures and a nonunit domain width test both the
  # gate scaling and the collapsed marginal factors rather than symmetry.
  H_left <- c(0.01, 0.3, 1.9, 0.65, 0.02)
  H_right <- c(2.3, 0.06, 0.1, 0.9, 0.001)
  for (m in c(0L, 1L, 4L, 8L)) {
    x <- runif(m, -1.5, 3.2)
    for (gate in c(0, 1e-10, 3.2, 50, 500)) {
      actual <- poistree:::SPPT_exact_cut_scores(
        x, cuts, gate, 4.7, H_left, H_right, 0.5, 0.17
      )
      expected <- lapply(seq_along(cuts), function(c) {
        cut_tilt_enumerate(x, cuts[c], gate, 4.7, H_left[c], H_right[c], 0.5, 0.17)
      })
      expected_pmf <- do.call(rbind, lapply(expected, `[[`, "log_pmf"))
      expected_scores <- vapply(expected, `[[`, numeric(1), "log_score")
      expect_lte(max(abs(actual$log_pmf - expected_pmf)), 1e-9)
      expect_lte(max(abs(actual$log_scores - expected_scores)), 1e-9)
      expect_lte(max(abs(actual$direct_log_pmf - expected_pmf)), 1e-9)
      expect_lte(max(abs(actual$direct_log_scores - expected_scores)), 1e-9)
      expect_lte(max(abs(actual$original_log_scores - expected_scores)), 1e-9)
      for (c in seq_along(cuts)) {
        count <- 0:m
        logQ <- function(n, H) 0.5 * log(0.17) - lgamma(0.5) +
          lgamma(0.5 + n) - (0.5 + n) * log(0.17 + H)
        expected_count <- expected[[c]]$log_pmf + logQ(count, H_left[c]) +
          logQ(m - count, H_right[c]) - expected_scores[c]
        expect_lte(max(abs(actual$log_allocation_count[c, ] - expected_count)), 1e-9)
      }
      expect_equal(rowSums(exp(actual$log_pmf)), rep(1, length(cuts)), tolerance = 1e-12)
    }
  }
})

test_that("cut tilting preserves very small count-probability tails", {
  x <- c(0, 0, 0.1, 0.4, 0.4, 0.9, 1, 1)
  cuts <- c(-0.8, 0, 0.5, 1, 1.8)
  for (gate in c(500, 2000)) {
    actual <- poistree:::SPPT_exact_cut_scores(
      x, cuts, gate, 1, c(0, 1e-12, 0.01, 0.2, 0.999),
      c(1, 0.999, 0.3, 1e-12, 0), 0.5, 0.0001
    )
    expect_true(all(is.finite(actual$log_pmf)))
    expect_true(all(is.finite(actual$direct_log_pmf)))
    expect_true(any(actual$log_pmf < log(.Machine$double.xmin)))
    expect_lte(max(abs(actual$log_pmf - actual$direct_log_pmf)), 1e-9)
    expect_lte(max(abs(actual$log_scores - actual$direct_log_scores)), 1e-9)
    for (c in seq_along(cuts)) {
      expected <- cut_tilt_enumerate(
        x, cuts[c], gate, 1, c(0, 1e-12, 0.01, 0.2, 0.999)[c],
        c(1, 0.999, 0.3, 1e-12, 0)[c], 0.5, 0.0001
      )
      expect_lte(max(abs(actual$log_pmf[c, ] - expected$log_pmf)), 1e-9)
      expect_lte(abs(actual$log_scores[c] - expected$log_score), 1e-9)
    }
  }
})

test_that("repeated locations recover the binomial count distribution", {
  x <- rep(0.27, 12)
  cuts <- c(0.1, 0.3, 0.75)
  gate <- 9
  width <- 2.5
  actual <- poistree:::SPPT_exact_cut_scores(
    x, cuts, gate, width, rep(0.2, 3), rep(0.8, 3), 2, 0.3
  )
  expected <- t(vapply(cuts, function(cut) {
    dbinom(0:length(x), length(x), plogis(gate * (cut - x[1]) / width), log = TRUE)
  }, numeric(length(x) + 1L)))
  expect_lte(max(abs(actual$log_pmf - expected)), 1e-11)
})

test_that("a 50-cut node matches separate exact recursions at the default threshold", {
  set.seed(411)
  x <- runif(150, -2, 3)
  cuts <- seq(-1.95, 2.95, length.out = 50)
  H_left <- seq(0.003, 1.8, length.out = 50)
  H_right <- rev(seq(0.003, 1.8, length.out = 50))
  for (gate in c(2, 50, 500)) {
    actual <- poistree:::SPPT_exact_cut_scores(x, cuts, gate, 5, H_left, H_right, 0.5, 0.1)
    expect_true(all(is.finite(actual$log_scores)))
    expect_lte(max(abs(actual$log_pmf - actual$direct_log_pmf)), 1e-8)
    expect_lte(max(abs(actual$log_scores - actual$direct_log_scores)), 1e-9)
    expect_equal(rowSums(exp(actual$log_pmf)), rep(1, 50), tolerance = 1e-10)
  }
})

test_that("cut tilting is invariant to observation and candidate order", {
  x <- c(0.12, 0.12, 0.28, 0.7, 0.99)
  cuts <- c(0.1, 0.35, 0.7, 0.35)
  left <- c(0.01, 0.2, 0.5, 0.2)
  right <- 1 - left
  actual <- poistree:::SPPT_exact_cut_scores(x, cuts, 40, 1, left, right, 0.5, 0.1)
  idx <- c(4, 2, 1, 3)
  permuted <- poistree:::SPPT_exact_cut_scores(rev(x), cuts[idx], 40, 1, left[idx], right[idx], 0.5, 0.1)
  expect_equal(permuted$log_scores, actual$log_scores[idx], tolerance = 1e-12)
  expect_equal(permuted$log_pmf, actual$log_pmf[idx, ], tolerance = 1e-12)
  expect_equal(actual$log_scores[2], actual$log_scores[4], tolerance = 1e-12)
  single <- poistree:::SPPT_exact_cut_scores(x, cuts[1], 40, 1, left[1], right[1], 0.5, 0.1)
  expect_equal(single$log_scores[1], actual$log_scores[1], tolerance = 1e-12)
})

test_that("exact allocation replay retains rare routing probabilities", {
  cases <- list(
    list(logits = numeric(), bits = integer()),
    list(logits = c(50, 1, 2), bits = c(0L, 1L, 0L)),
    list(logits = c(1000, 1, 2), bits = c(0L, 1L, 0L)),
    list(logits = c(-1000, -1, -2), bits = c(1L, 0L, 1L)),
    list(logits = c(1000, 1000), bits = c(0L, 1L)),
    list(logits = c(1000, -1000, 50, -50), bits = rep(0L, 4)),
    list(logits = c(1000, -1000, 50, -50), bits = rep(1L, 4))
  )
  for (case in cases) {
    result <- poistree:::SPPT_exact_allocation_logprob(case$logits, case$bits)
    expect_true(all(is.finite(unlist(result))))
    expect_lte(abs(result$log_conditional - result$direct_log_conditional), 1e-9)
    expect_lte(result$log_conditional, 1e-12)
    expected_joint <- sum(ifelse(
      case$bits == 1L,
      plogis(case$logits, log.p = TRUE),
      plogis(case$logits, lower.tail = FALSE, log.p = TRUE)
    ))
    expect_equal(result$log_joint, expected_joint, tolerance = 1e-12)
    if (all(case$bits == 0L) || all(case$bits == 1L))
      expect_equal(result$log_conditional, 0, tolerance = 1e-12)
  }
})

test_that("conditional allocation replay normalizes and is invariant to cut shifts", {
  logits <- c(-12, -3, 0.2, 0.2, 5, 14)
  bits <- as.matrix(expand.grid(rep(list(0:1), length(logits))))
  count <- rowSums(bits)
  probabilities <- vapply(seq_len(nrow(bits)), function(i) {
    result <- poistree:::SPPT_exact_allocation_logprob(logits, as.integer(bits[i, ]))
    shifted <- poistree:::SPPT_exact_allocation_logprob(logits + 400, as.integer(bits[i, ]))
    expect_lte(abs(result$log_conditional - result$direct_log_conditional), 1e-10)
    expect_lte(abs(shifted$log_conditional - result$log_conditional), 1e-9)
    result$log_conditional
  }, numeric(1))
  for (k in 0:length(logits)) {
    expect_equal(sum(exp(probabilities[count == k])), 1, tolerance = 1e-12)
  }
})

test_that("gate-table right probabilities retain subnormal log-left tails", {
  logits <- c(-1000, -50, 0, 50, 708, 720, 744, 745, 746, 1000)
  result <- poistree:::SPPT_exact_allocation_logprob(logits, integer(length(logits)))
  expected <- plogis(logits, lower.tail = FALSE, log.p = TRUE)
  expect_true(all(is.finite(result$gate_log_right)))
  expect_equal(result$gate_log_right, expected, tolerance = 1e-12)
})

test_that("shared allocation tables draw the independently enumerated collapsed law", {
  logits <- c(-1.7, -0.3, 0.8, 2)
  shifts <- c(-1.2, 0, 0.9)
  H_left <- c(0.4, 0.6, 0.9)
  H_right <- c(1, 0.6, 0.2)
  a <- 0.7
  b <- 0.3
  draws <- 30000L
  set.seed(68302)
  actual <- poistree:::SPPT_exact_reuse_draws(
    logits, shifts, H_left, H_right, a, b, draws
  )
  expect_equal(actual$prefix_builds, c(1L, 0L, 0L))
  expect_equal(actual$prefix_cells, (length(logits) + 1) * (length(logits) + 2) / 2)
  logQ <- function(k, H) a * log(b) - lgamma(a) + lgamma(a + k) - (a + k) * log(b + H)
  for (c in seq_along(shifts)) {
    oracle <- cut_tilt_enumerate(-logits, shifts[c], 1, 1, H_left[c], H_right[c], a, b)
    counts <- rowSums(oracle$allocations)
    probability <- exp(oracle$log_joint + logQ(counts, H_left[c]) +
      logQ(length(logits) - counts, H_right[c]) - oracle$log_score)
    encoded <- as.vector(actual$allocations[[c]] %*% 2^(seq_along(logits) - 1)) + 1L
    observed <- tabulate(encoded, nbins = 2^length(logits)) / draws
    tolerance <- 6 * sqrt(probability * (1 - probability) / draws) + 6 / draws
    expect_true(all(abs(observed - probability) < tolerance))
    expected_count <- vapply(0:length(logits), function(k) sum(probability[counts == k]), numeric(1))
    expect_equal(exp(actual$log_count[c, ]), expected_count, tolerance = 1e-12)
  }
})

test_that("conditional draws share one triangular table at fixed count", {
  logits <- c(-2, -0.1, 0.5, 1.6)
  shifts <- c(-1000, 0, 1000)
  draws <- 20000L
  set.seed(40384)
  actual <- poistree:::SPPT_exact_reuse_draws(
    logits, shifts, rep(0.3, 3), rep(0.7, 3), 0.5, 0.1, draws, fixed_count = 2L
  )
  oracle <- cut_tilt_enumerate(-logits, 0, 1, 1, 0.3, 0.7, 0.5, 0.1)
  count <- rowSums(oracle$allocations)
  probability <- ifelse(count == 2L, exp(oracle$log_joint - oracle$log_pmf[3]), 0)
  expect_equal(actual$prefix_builds, c(1L, 0L, 0L))
  expect_equal(actual$prefix_cells, 15L)
  for (B in actual$allocations) {
    expect_true(all(rowSums(B) == 2L))
    encoded <- as.vector(B %*% 2^(seq_along(logits) - 1)) + 1L
    observed <- tabulate(encoded, nbins = 16) / draws
    tolerance <- 6 * sqrt(probability * (1 - probability) / draws) + 6 / draws
    expect_true(all(abs(observed - probability) < tolerance))
  }
})

test_that("boundary allocation counts bypass the prefix table", {
  for (logits in list(numeric(), c(-1000, -2, 2, 1000))) {
    for (k in unique(c(0L, length(logits)))) {
      actual <- poistree:::SPPT_exact_reuse_draws(
        logits, c(-50, 50), c(0, 1), c(1, 0), 0.5, 0.1, 3L, fixed_count = k
      )
      expect_equal(actual$prefix_builds, c(0L, 0L))
      expect_equal(actual$prefix_cells, 0L)
      for (B in actual$allocations) {
        expect_equal(rowSums(B), rep(k, 3))
        expect_true(all(B == as.integer(k == length(logits))))
      }
    }
  }
})

test_that("exact allocation factors cancel with defensive action mixtures", {
  x <- c(0.03, 0.3, 0.8, 0.97)
  cuts <- c(0.2, 0.5, 0.8)
  H_left <- c(0.05, 0.4, 0.9)
  H_right <- c(0.95, 0.6, 0.1)
  a <- 0.5
  b <- 0.2
  logQ <- function(k, H) a * log(b) - lgamma(a) + lgamma(a + k) - (a + k) * log(b + H)
  score <- poistree:::SPPT_exact_cut_scores(x, cuts, 12, 1, H_left, H_right, a, b)$log_scores
  log_parent <- logQ(length(x), 1)
  prior <- c(0.3, 0.1, 0.2, 0.4)
  informed <- log(prior) + c(log_parent, score)
  informed <- exp(informed - cut_tilt_lse(informed))
  for (epsilon in c(0, 0.05, 0.8)) {
    action <- (1 - epsilon) * informed + epsilon * prior
    for (c in seq_along(cuts)) {
      oracle <- cut_tilt_enumerate(x, cuts[c], 12, 1, H_left[c], H_right[c], a, b)
      k <- rowSums(oracle$allocations)
      log_target <- oracle$log_joint + logQ(k, H_left[c]) + logQ(length(x) - k, H_right[c])
      log_allocation_proposal <- log_target - oracle$log_score
      full_ratio <- log(prior[c + 1L]) + log_target - log_parent -
        log(action[c + 1L]) - log_allocation_proposal
      simplified <- log(prior[c + 1L]) + score[c] - log_parent - log(action[c + 1L])
      expect_equal(full_ratio, rep(simplified, length(full_ratio)), tolerance = 1e-12)
    }
  }
})
