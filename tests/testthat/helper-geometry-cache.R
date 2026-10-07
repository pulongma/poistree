geometry_cache_args <- function() {
  list(
    x = cbind(c(.05, .09, .12, .19, .24, .32, .56, .66, .73, .85, .92, .96),
              c(.08, .15, .67, .22, .75, .29, .55, .89, .65, .42, .87, .94)),
    region = matrix(c(0, 1, 0, 1), 2L, byrow = TRUE),
    gating = "soft", sampler = "pcg", a = .5, b = .04,
    gate = 7, a_gate = 2, b_gate = .3, sd_gate = .5,
    alpha = .85, eta = .5, max_depth = 3L, cut_candidates = 4L,
    chains = 1L, iter = 60L, burn = 20L, thin = 2L,
    tree_moves = 2L, change_moves = 2L, seed = 219L, verbose = FALSE,
    predict_at = cbind(c(.03, .29, .63, .98), c(.21, .72, .37, .93)),
    test = cbind(c(.17, .44, .81), c(.32, .67, .54))
  )
}

expect_cache_same_fit <- function(cached, reference, cached_rng, reference_rng) {
  expect_identical(cached_rng, reference_rng)
  expect_equal(cached$prediction, reference$prediction, tolerance = 1e-10)
  expect_equal(cached$posterior, reference$posterior, tolerance = 1e-10)
  expect_equal(cached$diagnostics, reference$diagnostics, tolerance = 1e-10)
  expect_identical(cached$diagnostics$leaf_count_trace,
                   reference$diagnostics$leaf_count_trace)
  expect_equal(poistree:::ppt_integral(cached, type = "draws"),
               poistree:::ppt_integral(reference, type = "draws"), tolerance = 1e-10)
  at <- cbind(c(.01, .27, .51, .76, .99), c(.15, .44, .85, .23, .62))
  expect_equal(poistree:::.ppt_state_eval(cached, at),
               poistree:::.ppt_state_eval(reference, at), tolerance = 1e-10)
}

