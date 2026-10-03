test_that("hard-leaf candidate cuts default to no aspect-ratio cap", {
  region <- matrix(c(0, 0.1, 0, 1), ncol = 2, byrow = TRUE)
  points <- cbind(
    c(0.020, 0.021, 0.022, 0.023, 0.070, 0.080, 0.090, 0.095),
    seq(0.1, 0.8, length.out = 8)
  )

  unrestricted <- poistree:::PPT_valid_cuts(
    points, region, min_leaf_n = 2L, cut_grid_n = 30L,
    max_aspect_ratio = Inf, force_mid_cut = FALSE
  )
  compact <- poistree:::PPT_valid_cuts(
    points, region, min_leaf_n = 2L, cut_grid_n = 30L,
    max_aspect_ratio = 6, force_mid_cut = FALSE
  )

  expect_gt(length(unrestricted[[1L]]), 0L)
  expect_length(compact[[1L]], 0L)
})

test_that("finite aspect caps also govern fallback and midpoint cuts", {
  region <- matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE)
  points <- cbind(
    rep(c(0.49, 0.51), each = 4L),
    seq(0.1, 0.8, length.out = 8L)
  )

  fallback_allowed <- poistree:::PPT_valid_cuts(
    points, region, min_leaf_n = 2L, cut_grid_n = 30L,
    max_aspect_ratio = 2, force_mid_cut = FALSE
  )
  fallback_blocked <- poistree:::PPT_valid_cuts(
    points, region, min_leaf_n = 2L, cut_grid_n = 30L,
    max_aspect_ratio = 1.99, force_mid_cut = FALSE
  )
  midpoint_allowed <- poistree:::PPT_valid_cuts(
    points, region, min_leaf_n = 2L, cut_grid_n = 30L,
    max_aspect_ratio = 2, force_mid_cut = TRUE
  )
  midpoint_blocked <- poistree:::PPT_valid_cuts(
    points, region, min_leaf_n = 2L, cut_grid_n = 30L,
    max_aspect_ratio = 1.99, force_mid_cut = TRUE
  )

  expect_equal(fallback_allowed[[1L]], 0.5)
  expect_length(fallback_blocked[[1L]], 0L)
  expect_equal(midpoint_allowed[[1L]], 0.5)
  expect_length(midpoint_blocked[[1L]], 0L)
})

test_that("max_aspect_ratio is exposed by hard-leaf SMC and PGAS", {
  set.seed(920)
  x <- matrix(runif(40), ncol = 2)
  region <- matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE)

  smc <- ppt_fit(
    x, region,
    gating = "hard", scales = "leaf", sampler = "smc",
    max_depth = 1L, min_leaf_n = 2L, particles = 10L,
    max_aspect_ratio = Inf, seed = 921L
  )
  expect_identical(smc$control$max_aspect_ratio, Inf)

  pgas <- ppt_fit(
    x, region,
    gating = "hard", scales = "leaf", sampler = "pgas",
    max_depth = 1L, min_leaf_n = 2L, particles = 4L,
    chains = 1L, iter = 2L, burn = 1L,
    max_aspect_ratio = 12, seed = 922L, verbose = FALSE
  )
  expect_identical(pgas$control$max_aspect_ratio, 12)

  expect_error(
    ppt_fit(
      x, region,
      gating = "hard", scales = "leaf", sampler = "smc",
      particles = 10L, max_aspect_ratio = 0.5
    ),
    "max_aspect_ratio"
  )
})
