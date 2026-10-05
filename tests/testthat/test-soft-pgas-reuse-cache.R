test_that("root payloads invalidate only changed coordinates and retain lazy tables", {
  gates <- rbind(c(4, 7), c(4, 7), c(5, 7), c(5, 8), c(4, 8), c(4, 8))
  for (exact_max in c(0L, 150L)) {
    cached <- poistree:::SPPT_root_reuse_probe(gates, exact_max, TRUE)
    fresh <- poistree:::SPPT_root_reuse_probe(gates, exact_max, FALSE)
    expect_identical(cached$output, fresh$output)
    expect_equal(cached$builds, c(2, 2, 3, 4, 5, 5))
    expect_equal(cached$hits, c(0, 2, 3, 4, 5, 7))
    expect_equal(fresh$builds, 2 * seq_len(nrow(gates)))
    expect_equal(fresh$hits, rep(0, nrow(gates)))
    expect_equal(cached$prefixes, if (exact_max > 0) cached$builds else rep(0, nrow(gates)))
    expect_equal(fresh$prefixes, if (exact_max > 0) fresh$builds else rep(0, nrow(gates)))
  }
})

test_that("root caching preserves full PGAS paths with fixed and changing gates", {
  for (exact_max in c(0L, 150L)) for (config in list(c(FALSE, FALSE), c(TRUE, FALSE), c(TRUE, TRUE))) {
    set.seed(503)
    cached <- poistree:::SPPT_root_reuse_chain(TRUE, config[1], config[2], exact_max)
    set.seed(503)
    fresh <- poistree:::SPPT_root_reuse_chain(FALSE, config[1], config[2], exact_max)
    expect_identical(cached, fresh)
  }
})

test_that("forced exact allocations route without replay or random draws", {
  for (defensive in c(0, 0.2, 0.9)) for (bits in list(rep(0L, 6), rep(1L, 6), c(0L, 1L, 0L, 1L, 0L, 1L))) {
    set.seed(123)
    rng <- .Random.seed
    result <- poistree:::SPPT_forced_exact_probe(bits, defensive)
    expect_identical(.Random.seed, rng)
    expect_equal(result$left, which(bits == 1L) - 1L)
    expect_equal(result$right, which(bits == 0L) - 1L)
    expect_equal(result$prefix_builds, 0)
    expect_equal(result$forced_routes, 1)
    # Root split at coordinate 1/cut 0.5 is action index 2 in C++, 3 in R.
    expect_equal(result$logw, result$score[3] - result$logQA - result$logq[3], tolerance = 1e-12)
    if (defensive == 0) {
      logZ <- max(result$score) + log(sum(exp(result$score - max(result$score))))
      expect_equal(result$logw, logZ - result$logQA, tolerance = 1e-12)
    }
  }
})
