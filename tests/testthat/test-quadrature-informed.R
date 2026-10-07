quadrature_informed_args <- function() {
  background <- cbind(
    c(0, .04, .16, .23, .31, .42, .53, .64, .76, .88, 1),
    c(1, .12, .73, .36, .84, .28, .55, .06, .66, .42, 0)
  )
  list(
    x = cbind(c(.05, .09, .12, .19, .24, .32, .56, .66, .73, .85, .92, .96),
              c(.08, .15, .67, .22, .75, .29, .55, .89, .65, .42, .87, .94)),
    region = matrix(c(0, 1, 0, 1), 2L, byrow = TRUE),
    background = background,
    weights = c(.01, .04, .02, .15, .07, .13, .09, .18, .11, .08, .12) * 2.7,
    gating = "soft", sampler = "pcg", a = .5,
    gate = 7, a_gate = 2, b_gate = .3, sd_gate = .5,
    alpha = .85, eta = .5, max_depth = 3L, cut_candidates = 4L,
    chains = 1L, iter = 60L, burn = 20L, thin = 2L,
    tree_moves = 2L, change_moves = 2L, seed = 219L, verbose = FALSE,
    predict_at = background,
    test = cbind(c(.17, .44, .81), c(.32, .67, .54))
  )
}

quadrature_informed_soft_draws <- function(fit, at) {
  evaluate <- function(s) {
    nodes <- fit$posterior$state$nodes[[s]]
    gate <- fit$posterior$state$gate[s, ]
    value <- numeric(nrow(at))
    visit <- function(id, mass, box) {
      node <- nodes[match(id, nodes[, 1L]), ]
      if (node[2L] < 0) {
        value <<- value + mass * node[4L]
        return(invisible(NULL))
      }
      axis <- as.integer(node[2L]) + 1L
      cut <- node[3L]
      width <- if (fit$model$gate_scale == "node") diff(box[axis, ]) else
        diff(fit$data$region[axis, ])
      right <- plogis(gate[axis] * (at[, axis] - cut) / width)
      left_box <- right_box <- box
      left_box[axis, 2L] <- right_box[axis, 1L] <- cut
      visit(2L * id, mass * (1 - right), left_box)
      visit(2L * id + 1L, mass * right, right_box)
    }
    visit(1L, rep(1, nrow(at)), fit$data$region)
    value
  }
  vapply(seq_len(fit$posterior$draws), evaluate, numeric(nrow(at)))
}

quadrature_informed_inside <- function(at, box, region) {
  inside <- rep(TRUE, nrow(at))
  for (j in seq_len(ncol(at))) {
    upper <- if (box[j, 2L] == region[j, 2L]) at[, j] <= box[j, 2L] else
      at[, j] < box[j, 2L]
    inside <- inside & at[, j] >= box[j, 1L] & upper
  }
  inside
}

quadrature_informed_hard_draws <- function(fit, at) {
  evaluate <- function(tree) {
    value <- numeric(nrow(at))
    membership <- integer(nrow(at))
    for (leaf in Filter(function(node) isTRUE(node$is_leaf), tree)) {
      inside <- quadrature_informed_inside(at, leaf$region, fit$data$region)
      value[inside] <- leaf$lambda
      membership <- membership + inside
    }
    expect_identical(membership, rep(1L, nrow(at)))
    value
  }
  vapply(fit$posterior$tree_draws, evaluate, numeric(nrow(at)))
}

quadrature_informed_lppd <- function(draws, integral, weights) {
  ll <- colSums(log(draws)) - integral
  center <- max(ll)
  center + log(sum(weights * exp(ll - center)))
}

quadrature_informed_legacy_fit <- function(args) {
  poistree:::qpp_set_quadrature(args$background, args$weights, args$region)
  on.exit(poistree:::qpp_clear_quadrature(), add = TRUE)
  args$b <- args$a * sum(args$weights) / nrow(args$x)
  args$background <- args$weights <- NULL
  do.call(ppt_fit, args)
}

test_that("quadrature informed PCG preserves weighted targets for all gate scales", {
  for (scale in c("root", "node")) {
    for (structure in c("shared", "dimension")) {
      args <- quadrature_informed_args()
      args$gate_scale <- scale
      args$gate_structure <- structure
      fit <- do.call(ppt_fit_quadrature, c(args, list(
        informed = TRUE, proposal_temperature = .73, proposal_defensive = .23
      )))
      expect_s3_class(fit, "ppt")
      expect_true(fit$control$informed)
      expect_identical(fit$control$proposal_temperature, .73)
      expect_identical(fit$control$proposal_defensive, .23)
      expect_identical(fit$model$algorithm,
                       "Partially collapsed Gibbs (RAM, informed tree proposals)")
      expect_identical(fit$model$gate_scale, scale)
      expect_identical(fit$prior$gate$structure, structure)
      expect_equal(fit$prior$intensity$rate,
                   args$a * sum(args$weights) / nrow(args$x))
      expect_equal(fit$data$quadrature$total_exposure, sum(args$weights))
      expect_true(any(fit$diagnostics$leaf_count_trace > 1L))
      if (structure == "shared") {
        expect_equal(fit$posterior$state$gate[, 1L],
                     fit$posterior$state$gate[, 2L])
      }
      background <- quadrature_informed_soft_draws(fit, args$background)
      integral <- colSums(background * args$weights)
      expect_equal(background, fit$prediction$draws, tolerance = 1e-10)
      expect_equal(integral, fit$posterior$integrated_intensity_draws,
                   tolerance = 1e-10)
      expect_equal(poistree:::ppt_integral(fit, type = "draws"), integral, tolerance = 1e-10)
      train <- quadrature_informed_soft_draws(fit, args$x)
      expect_equal(mean(colSums(log(train)) - integral),
                   fit$posterior$mean_log_likelihood, tolerance = 1e-10)
      test <- quadrature_informed_soft_draws(fit, args$test)
      expected <- quadrature_informed_lppd(test, integral,
                                           rep(1 / ncol(test), ncol(test)))
      expect_equal(fit$posterior$lppd, expected, tolerance = 1e-10)
      expect_equal(as.numeric(ppt_lppd(fit, test = args$test)), expected,
                   tolerance = 1e-10)
    }
  }
})

test_that("quadrature defaults preserve the original uninformed PCG path", {
  args <- quadrature_informed_args()
  default <- do.call(ppt_fit_quadrature, args)
  default_rng <- .Random.seed
  explicit <- do.call(ppt_fit_quadrature, c(args, list(informed = FALSE)))
  expect_identical(default_rng, .Random.seed)
  default$call <- explicit$call <- NULL
  expect_identical(default, explicit)
  ignored <- do.call(ppt_fit_quadrature, c(args, list(
    informed = FALSE, proposal_temperature = .73, proposal_defensive = .23
  )))
  ignored$call <- NULL
  expect_identical(default, ignored)
  legacy <- quadrature_informed_legacy_fit(args)
  expect_identical(default_rng, .Random.seed)
  for (component in c("prediction", "posterior", "diagnostics", "prior")) {
    expect_identical(default[[component]], legacy[[component]])
  }
  expect_false(default$control$informed)
  expect_null(default$control$proposal_temperature)
  expect_null(default$control$proposal_defensive)

  formal_names <- names(formals(ppt_fit_quadrature))
  expect_identical(formal_names[seq_len(9L)],
                   c("x", "region", "background", "weights", "gating",
                     "sampler", "a", "b", "..."))
  expect_true(all(match(c("informed", "proposal_temperature", "proposal_defensive"),
                        formal_names) > match("...", formal_names)))
})

test_that("quadrature errors clear the active integration rule", {
  args <- quadrature_informed_args()
  plain_args <- args
  plain_args$background <- plain_args$weights <- NULL
  plain_args$b <- .04
  before <- do.call(ppt_fit, plain_args)
  before_rng <- .Random.seed

  expect_error(do.call(ppt_fit_quadrature, c(args, list(
    informed = TRUE, cache_geometry = NA
  ))), "cache_geometry")
  after <- do.call(ppt_fit, plain_args)
  expect_identical(before_rng, .Random.seed)
  before$call <- after$call <- NULL
  expect_identical(before, after)

  again <- do.call(ppt_fit_quadrature, c(args, list(informed = TRUE)))
  expect_true(all(is.finite(again$posterior$integrated_intensity_draws)))
})

test_that("quadrature validates informed controls and supported backends", {
  args <- quadrature_informed_args()
  expect_error(do.call(ppt_fit_quadrature, c(args, list(informed = NA))), "informed")
  expect_error(do.call(ppt_fit_quadrature, c(args, list(
    informed = TRUE, proposal_temperature = 0
  ))), "proposal_temperature")
  expect_error(do.call(ppt_fit_quadrature, c(args, list(
    informed = TRUE, proposal_defensive = 1
  ))), "proposal_defensive")
  minimal <- args[c("x", "region", "background", "weights")]
  expect_error(do.call(ppt_fit_quadrature, c(minimal, list(
    gating = "hard", sampler = "smc", informed = TRUE
  ))), "[Ii]nformed")
  expect_error(do.call(ppt_fit_quadrature, c(minimal, list(
    gating = "hard", sampler = "smc", engine = "dense"
  ))), "shared")
  for (gating in c("soft", "hard")) {
    for (sampler in c("pgas", "rjmcmc", "irjmcmc")) {
      expect_error(do.call(ppt_fit_quadrature, c(minimal, list(
        gating = gating, sampler = sampler
      ))), "soft PCG or hard SMC")
    }
  }
})

test_that("hard shared SMC integrates weights once at cut ties and outer edges", {
  x <- matrix(c(rep(.1, 8L), rep(.2, 8L), rep(.8, 2L), .9), ncol = 1L)
  region <- matrix(c(0, 1), nrow = 1L)
  background <- matrix(c(0, .1, .2, .2, .5, .8, .9, 1), ncol = 1L)
  weights <- c(.05, .1, .1, .2, 1, .6, .25, .4)
  test <- matrix(c(0, .2, .8, 1), ncol = 1L)
  fit <- ppt_fit_quadrature(
    x, region, background, weights, "hard", "smc", .5, .07,
    predict_at = background, test = test, max_depth = 2L,
    particles = 40L, cut_candidates = 5L, seed = 29L
  )
  expect_identical(fit$backend, "PPT_fit_SMC_shared")
  expect_identical(fit$control$engine, "shared")
  n <- nrow(x)
  log_root <- lgamma(n + .5) - lgamma(.5) + .5 * log(.07) -
    (n + .5) * log(.07 + sum(weights))
  expect_equal(fit$posterior$log_target_normalizer,
               fit$posterior$log_relative_normalizer + log_root,
               tolerance = 1e-12)
  expect_true(any(fit$diagnostics$leaf_count_trace > 1L))
  for (tree in fit$posterior$tree_draws) {
    leaves <- Filter(function(node) isTRUE(node$is_leaf), tree)
    exposures <- vapply(leaves, function(leaf) {
      inside <- quadrature_informed_inside(background, leaf$region, region)
      expected <- sum(weights[inside])
      expect_equal(leaf$exposure, expected, tolerance = 1e-12)
      expected
    }, numeric(1))
    expect_equal(sum(exposures), sum(weights), tolerance = 1e-12)
    for (node in Filter(function(node) !isTRUE(node$is_leaf), tree)) {

      expect_true(any(background[, 1L] == node$L))
    }
  }
  draws <- quadrature_informed_hard_draws(fit, background)
  integral <- colSums(draws * weights)
  expect_equal(draws, fit$prediction$draws, tolerance = 1e-10)
  expect_equal(integral, fit$posterior$integrated_intensity_draws, tolerance = 1e-10)
  expect_equal(poistree:::ppt_integral(fit, type = "draws"), integral, tolerance = 1e-10)
  expected <- quadrature_informed_lppd(
    quadrature_informed_hard_draws(fit, test), integral, fit$posterior$particle_weights
  )
  expect_equal(fit$posterior$lppd, expected, tolerance = 1e-10)
  expect_equal(as.numeric(ppt_lppd(fit, test = test)), expected, tolerance = 1e-10)
})

test_that("compact informed quadrature agrees with uncached exposure", {
  args <- quadrature_informed_args()
  args$gate_family <- "compact"
  args$gate_scale <- "node"
  for (structure in c("shared", "dimension")) {
    args$gate_structure <- structure
    cached <- do.call(ppt_fit_quadrature, c(args, list(
      informed = TRUE, cache_geometry = TRUE)))
    rng <- .Random.seed
    direct <- do.call(ppt_fit_quadrature, c(args, list(
      informed = TRUE, cache_geometry = FALSE)))
    expect_identical(rng, .Random.seed)
    expect_equal(cached$posterior, direct$posterior, tolerance = 1e-10)
    expect_equal(colSums(cached$prediction$draws * args$weights),
                 cached$posterior$integrated_intensity_draws, tolerance = 1e-10)
    expect_true(cached$control$informed)
  }
})
