test_that("hard leaf RJ-MCMC uses the unified ppt API", {
  set.seed(10)
  x <- matrix(runif(80), ncol = 2)
  colnames(x) <- c("x", "y")
  region <- matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE)
  grid <- matrix(runif(30), ncol = 2)
  test <- matrix(runif(10), ncol = 2)

  fit <- ppt_fit(
    x, region,
    gating = "hard", scales = "leaf", sampler = "rjmcmc",
    predict_at = grid, test = test,
    max_depth = 3, min_leaf_n = 3,
    chains = 1, iter = 50, burn = 10,
    cut_candidates = 5, prediction_draws = 10,
    seed = 11, verbose = FALSE
  )

  expect_s3_class(fit, "ppt")
  expect_identical(fit$model$label, "PPT")
  expect_identical(fit$model$gating, "hard")
  expect_identical(fit$model$scales, "leaf")
  expect_identical(fit$model$sampler, "rjmcmc")
  expect_length(ppt_predict(fit), nrow(grid))
  expect_equal(nrow(ppt_predict(fit, type = "interval")), nrow(grid))
  expect_true(all(is.finite(ppt_predict(fit))))
  expect_true(is.finite(as.numeric(poistree:::ppt_logLik(fit))))
  expect_true(is.finite(as.numeric(ppt_lppd(fit))))
  expect_true(is.na(fit$posterior$log_evidence))
  expect_true(is.finite(fit$posterior$mean_leaves))
  expect_true(is.finite(fit$posterior$mean_max_depth))
  expect_true(is.finite(poistree:::ppt_integral(fit)))
  expect_length(poistree:::ppt_integral(fit, "draws"), fit$posterior$draws)
  expect_equal(nrow(fit$prediction$draws), nrow(grid))
  expect_equal(ncol(fit$prediction$draws), fit$posterior$draws)
  expect_length(
    ppt_diagnostics(fit)$leaf_count_trace,
    fit$control$chains * (fit$control$iter - fit$control$burn)
  )
  expect_named(
    fit$diagnostics$tree_acceptance,
    c("grow", "prune", "change")
  )
  expect_length(fit$posterior$tree_draws, fit$posterior$draws)
  marginal_draws <- ppt_marginal(
    fit, "x", grid = grid[, 1L], type = "draws"
  )
  expect_true(all(is.finite(marginal_draws)))
  x_edges <- sort(unique(c(
    region[1L, ],
    unlist(lapply(fit$posterior$tree_draws, function(tree) {
      unlist(lapply(Filter(function(node) isTRUE(node$is_leaf), tree),
                    function(node) as.matrix(node$region)[1L, ]))
    }))
  )))
  x_width <- diff(x_edges)
  exact_marginal <- ppt_marginal(
    fit, "x", grid = head(x_edges, -1L) + x_width / 2,
    average = FALSE, type = "draws"
  )
  expect_equal(
    unname(colSums(exact_marginal * x_width)),
    poistree:::ppt_integral(fit, "draws"), tolerance = 1e-10
  )
})

test_that("hard leaf PGAS uses the unified ppt API", {
  set.seed(20)
  x <- matrix(runif(40), ncol = 2)
  colnames(x) <- c("x", "y")
  region <- matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE)
  grid <- matrix(runif(12), ncol = 2)
  test <- matrix(runif(6), ncol = 2)

  fit <- ppt_fit(
    x, region,
    gating = "hard", scales = "leaf", sampler = "pgas",
    predict_at = grid, test = test,
    max_depth = 2, min_leaf_n = 2,
    particles = 10, chains = 1, iter = 6, burn = 2,
    seed = 21, verbose = FALSE
  )

  expect_s3_class(fit, "ppt")
  expect_identical(fit$model$label, "PPT")
  expect_identical(fit$model$gating, "hard")
  expect_identical(fit$model$scales, "leaf")
  expect_identical(fit$model$sampler, "pgas")
  expect_length(ppt_predict(fit), nrow(grid))
  expect_true(all(is.finite(ppt_predict(fit))))
  expect_true(is.finite(as.numeric(poistree:::ppt_logLik(fit))))
  expect_true(is.finite(as.numeric(ppt_lppd(fit))))
  expect_true(is.na(fit$posterior$log_evidence))
  expect_true(is.finite(fit$posterior$mean_leaves))
  expect_true(is.finite(fit$posterior$mean_max_depth))
  expect_true(is.finite(poistree:::ppt_integral(fit)))
  expect_length(poistree:::ppt_integral(fit, "draws"), fit$posterior$draws)
  expect_true(is.finite(fit$diagnostics$particle_ess))
  expect_true(isTRUE(fit$control$conditional_smc))
  expect_false(isTRUE(fit$control$ancestor_sampling))
  expect_identical(fit$control$resampling_schedule, "tree_level")
  expect_equal(dim(fit$prediction$draws),
               c(nrow(grid), fit$posterior$draws))
  expect_length(fit$posterior$tree_draws, fit$posterior$draws)
  marginal_draws <- ppt_marginal(
    fit, "x", grid = grid[, 1L], type = "draws"
  )
  expect_true(all(is.finite(marginal_draws)))
  x_edges <- sort(unique(c(
    region[1L, ],
    unlist(lapply(fit$posterior$tree_draws, function(tree) {
      unlist(lapply(Filter(function(node) isTRUE(node$is_leaf), tree),
                    function(node) as.matrix(node$region)[1L, ]))
    }))
  )))
  x_width <- diff(x_edges)
  exact_marginal <- ppt_marginal(
    fit, "x", grid = head(x_edges, -1L) + x_width / 2,
    average = FALSE, type = "draws"
  )
  expect_equal(
    unname(colSums(exact_marginal * x_width)),
    poistree:::ppt_integral(fit, "draws"), tolerance = 1e-10
  )
  expect_length(
    ppt_diagnostics(fit)$leaf_count_trace,
    fit$posterior$draws
  )
})

test_that("transition replay law matches the one-step-ahead proposal", {
  cases <- list(
    list(
      x = matrix(seq(0.05, 0.95, length.out = 20L), ncol = 1L),
      region = matrix(c(0, 1), nrow = 1L)
    ),
    list(
      x = matrix(
        ((seq_len(90L) * c(17, 31, 43)) %% 97) / 97,
        ncol = 3L, byrow = TRUE
      ),
      region = matrix(rep(c(0, 1), 3L), ncol = 2L, byrow = TRUE)
    )
  )

  for (case in cases) {
    law <- poistree:::PPT_transition_probabilities(
      case$x, case$region,
      min_leaf_n = 1L, cut_grid_n = 5L,
      max_aspect_ratio = Inf
    )

    expect_gt(sum(law$action == "split"), 1L)
    expect_true(all(is.finite(law$log_probability)))
    expect_equal(sum(law$probability), 1, tolerance = 1e-12)
    expect_equal(
      law$replay_split_probability,
      law$proposal_split_probability,
      tolerance = 1e-12
    )
  }
})

test_that("serialized split probability is a probability rather than S", {
  x <- matrix(seq(0.03, 0.97, length.out = 30L), ncol = 1L)
  region <- matrix(c(0, 1), nrow = 1L)
  raw <- poistree:::PPT_fit_SMC(
    x, matrix(0.5, ncol = 1L), region,
    max_depth = 1L, P = 12L, min_leaf_n = 1L,
    resample_thresh = 0.5, a = 0.5, b = 0,
    max_aspect_ratio = Inf
  )
  probability <- vapply(
    raw$particle, function(tree) as.numeric(tree[[1L]]$prob_split),
    numeric(1)
  )
  decision <- vapply(
    raw$particle, function(tree) as.numeric(tree[[1L]]$S), numeric(1)
  )

  expect_true(all(is.finite(probability)))
  expect_true(all(probability > 0 & probability < 1))
  expect_false(isTRUE(all.equal(probability, decision)))
})

test_that("conditional SMC preserves the complete reference when P equals one", {
  x <- matrix(
    c(seq(0.02, 0.16, length.out = 20L),
      seq(0.82, 0.98, length.out = 20L)),
    ncol = 1L
  )
  region <- matrix(c(0, 1), nrow = 1L)

  set.seed(1)
  raw <- poistree:::PPT_fit_PG(
    x, matrix(0.5, ncol = 1L), region,
    max_depth = 2L, niter = 8L, P = 1L,
    min_leaf_n = 2L, resample_thresh = 0.5,
    a = 0.5, b = 0, verbose = FALSE,
    max_aspect_ratio = Inf
  )

  topology_signature <- function(tree) {
    paste(vapply(tree, function(node) {
      if (is.null(node)) return("_")
      paste(
        as.integer(node$S), as.integer(node$J),
        format(as.numeric(node$L), digits = 16), sep = ":"
      )
    }, character(1)), collapse = "|")
  }
  signatures <- vapply(
    raw$particles, topology_signature, character(1)
  )

  expect_identical(as.integer(raw$particles[[1L]][[1L]]$S), 1L)
  expect_length(unique(signatures), 1L)
  expect_equal(as.numeric(raw$weights), rep(1, ncol(raw$weights)))
})

test_that("Particle Gibbs honors a proper nondefault Gamma leaf prior", {
  x <- matrix(seq(2.1, 3.9, length.out = 6L), ncol = 1L)
  region <- matrix(c(2, 4), nrow = 1L)
  a <- 2.3
  b <- 4.0

  set.seed(730)
  raw <- poistree:::PPT_fit_PG(
    x, matrix(3, ncol = 1L), region,
    max_depth = 0L, niter = 8000L, P = 1L,
    min_leaf_n = 1L, resample_thresh = 0.5,
    a = a, b = b, verbose = FALSE,
    max_aspect_ratio = Inf
  )
  lambda <- as.numeric(raw$lambda$draws)

  expected_mean <- (a + nrow(x)) / (b + 2)
  expected_variance <- (a + nrow(x)) / (b + 2)^2
  expect_equal(mean(lambda), expected_mean, tolerance = 0.03)
  expect_equal(stats::var(lambda), expected_variance, tolerance = 0.03)
})

test_that("depth-one Particle Gibbs has the enumerated posterior law", {
  x <- matrix(
    c(seq(0.03, 0.23, length.out = 15L),
      seq(0.75, 0.97, length.out = 15L)),
    ncol = 1L
  )
  region <- matrix(c(0, 1), nrow = 1L)
  law <- poistree:::PPT_transition_probabilities(
    x, region, min_leaf_n = 1L, cut_grid_n = 30L,
    max_aspect_ratio = Inf
  )

  set.seed(11)
  raw <- poistree:::PPT_fit_PG(
    x, matrix(0.5, ncol = 1L), region,
    max_depth = 1L, niter = 12000L, P = 4L,
    min_leaf_n = 1L, resample_thresh = 0.5,
    a = 0.5, b = 0, verbose = FALSE,
    max_aspect_ratio = Inf, cut_grid_n = 30L
  )
  roots <- lapply(raw$particles, `[[`, 1L)

  observed <- numeric(length(law$probability))
  observed[law$action == "stop"] <- mean(vapply(
    roots, function(node) node$S == 0L, logical(1)
  ))
  split_rows <- which(law$action == "split")
  for (k in split_rows) {
    observed[k] <- mean(vapply(roots, function(node) {
      node$S == 1L && node$J == law$axis[k] - 1L &&
        abs(node$L - law$cut[k]) < 1e-10
    }, logical(1)))
  }

  expect_equal(sum(observed), 1, tolerance = 1e-12)
  expect_lt(max(abs(observed - law$probability)), 0.02)
})

test_that("conditional-SMC Particle Gibbs returns valid tree partitions", {
  set.seed(912)
  x <- cbind(runif(50), runif(50))
  region <- matrix(rep(c(0, 1), 2L), ncol = 2L, byrow = TRUE)
  fit <- ppt_fit(
    x, region, gating = "hard", scales = "leaf", sampler = "pgas",
    max_depth = 3L, min_leaf_n = 1L,
    particles = 12L, chains = 1L, iter = 14L, burn = 2L,
    seed = 913L, verbose = FALSE
  )

  check_tree <- function(tree) {
    leaves <- integer()
    for (i in seq_along(tree)) {
      node <- tree[[i]]
      if (is.null(node)) next

      if (i > 1L) {
        parent <- tree[[i %/% 2L]]
        expect_false(is.null(parent))
        expect_false(isTRUE(parent$is_leaf))
        expect_equal(node$depth, parent$depth + 1L)
      }

      box <- as.matrix(node$region)
      expect_true(all(is.finite(box)))
      expect_true(all(box[, 2L] > box[, 1L]))

      if (isTRUE(node$is_leaf)) {
        leaves <- c(leaves, i)
        next
      }

      expect_identical(as.integer(node$S), 1L)
      left_id <- 2L * i
      right_id <- left_id + 1L
      expect_lte(right_id, length(tree))
      left <- tree[[left_id]]
      right <- tree[[right_id]]
      expect_false(is.null(left))
      expect_false(is.null(right))
      expect_equal(
        sort(c(as.integer(left$idx), as.integer(right$idx))),
        sort(as.integer(node$idx))
      )
      expect_length(intersect(left$idx, right$idx), 0L)

      axis <- as.integer(node$J) + 1L
      expect_gte(axis, 1L)
      expect_lte(axis, nrow(box))
      expect_equal(as.matrix(left$region)[axis, 2L], node$L)
      expect_equal(as.matrix(right$region)[axis, 1L], node$L)
    }

    leaf_nodes <- tree[leaves]
    expect_equal(
      sort(unlist(lapply(leaf_nodes, `[[`, "idx"), use.names = FALSE)),
      seq_len(nrow(x))
    )
    leaf_volume <- sum(vapply(leaf_nodes, function(node) {
      box <- as.matrix(node$region)
      prod(box[, 2L] - box[, 1L])
    }, numeric(1)))
    expect_equal(leaf_volume, 1, tolerance = 1e-12)
  }

  invisible(lapply(fit$posterior$tree_draws, check_tree))
})

test_that("hard and soft samplers return the same unified schema", {
  x <- matrix(seq(0.05, 0.95, length.out = 20), ncol = 1)
  region <- matrix(c(0, 1), nrow = 1)
  grid <- matrix(seq(0.1, 0.9, length.out = 5), ncol = 1)

  smc <- ppt_fit(
    x, region, gating = "hard", scales = "leaf", sampler = "smc",
    predict_at = grid, max_depth = 2, min_leaf_n = 2,
    particles = 10, seed = 12
  )
  rjmcmc <- ppt_fit(
    x, region, gating = "hard", scales = "leaf", sampler = "rjmcmc",
    predict_at = grid, max_depth = 2, min_leaf_n = 2,
    chains = 1, iter = 30, burn = 10,
    cut_candidates = 5, prediction_draws = 5,
    seed = 13, verbose = FALSE
  )
  pgas <- ppt_fit(
    x, region, gating = "hard", scales = "leaf", sampler = "pgas",
    predict_at = grid, max_depth = 2, min_leaf_n = 2,
    particles = 8, chains = 1, iter = 5, burn = 2,
    seed = 14, verbose = FALSE
  )

  soft <- ppt_fit(
    x, region, gating = "soft", scales = "leaf", sampler = "rjmcmc",
    predict_at = grid, gate = 10, update_gate = FALSE,
    max_depth = 2, min_leaf_n = 3,
    chains = 1, iter = 20, burn = 5, thin = 2,
    tree_moves = 1, change_moves = 1, cut_candidates = 4,
    seed = 14, verbose = FALSE
  )

  for (fit in list(rjmcmc, pgas, soft)) {
    expect_identical(names(fit), names(smc))
    expect_identical(names(fit$posterior), names(smc$posterior))
    expect_identical(names(fit$diagnostics), names(smc$diagnostics))
  }
  for (fit in list(smc, soft)) {
    expect_true(is.finite(fit$posterior$mean_max_depth))
    expect_true(is.finite(poistree:::ppt_integral(fit)))
    expect_length(poistree:::ppt_integral(fit, "draws"), fit$posterior$draws)
  }
  expect_true(is.na(smc$posterior$mean_gate))
})

test_that("only terminal-leaf scales are accepted", {
  x <- matrix(seq(0.1, 0.9, length.out = 20), ncol = 1)
  region <- matrix(c(0, 1), nrow = 1)

  expect_error(
    ppt_fit(x, region, gating = "hard", scales = "other"),
    "'arg' should be"
  )
  expect_error(
    ppt_fit(x, region, gating = "soft", sampler = "smc"),
    "not available"
  )
})

test_that("hard RJ-MCMC snaps boundary roundoff safely", {
  x <- matrix(c(0.2, 0.8, 1 + 5e-14), ncol = 1)
  region <- matrix(c(0, 1), nrow = 1)
  fit <- ppt_fit(
    x, region, gating = "hard", scales = "leaf", sampler = "rjmcmc",
    max_depth = 1, min_leaf_n = 1,
    chains = 1, iter = 20, burn = 5,
    cut_candidates = 3, prediction_draws = 5,
    verbose = FALSE
  )
  expect_equal(max(fit$data$x), 1)
})
