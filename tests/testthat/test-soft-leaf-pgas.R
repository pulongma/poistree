# Exact posterior over (tree, labels) for d = 1, max_depth = 2 and a fixed grid,
# used to check that the soft Particle-Gibbs kernel is invariant.
soft_exact_posterior <- function(x, grid, region, a, b, gate, alpha, eta) {
  n <- nrow(x)
  rho <- function(depth) alpha * (1 + depth)^(-eta)          # split probability
  exposure <- function(path) {
    if (!nrow(path)) return(region[1, 2] - region[1, 1])
    poistree:::ppstree_geometry(
      rep(1L, nrow(path)), path$cut, rep(region[1, 2] - region[1, 1], nrow(path)),
      path$side, matrix(0.5, 1, 1), region, gate
    )$H
  }
  log_phi <- function(path, xi) {
    if (!nrow(path)) return(0)
    z <- gate * (xi - path$cut) / (region[1, 2] - region[1, 1])
    sum(ifelse(path$side < 0, -log1p(exp(z)), -log1p(exp(-z))))
  }
  log_Q <- function(m, H) a * log(b) - lgamma(a) + lgamma(a + m) - (a + m) * log(b + H)
  node_options <- c(list(list(split = FALSE, cut = NA)),
                    lapply(grid, function(c) list(split = TRUE, cut = c)))
  trees <- list(list(root = list(split = FALSE, cut = NA)))
  for (c1 in grid) for (o2 in node_options) for (o3 in node_options) {
    trees[[length(trees) + 1L]] <- list(
      root = list(split = TRUE, cut = c1), left = o2, right = o3
    )
  }
  summaries <- list()
  for (tree in trees) {
    leaves <- list(); lp <- 0
    if (!tree$root$split) {
      leaves[["1"]] <- data.frame(cut = numeric(), side = numeric())
      lp <- log(1 - rho(0))
    } else {
      lp <- log(rho(0)) - log(length(grid))
      for (child in c("left", "right")) {
        side <- if (child == "left") -1 else 1
        id <- if (child == "left") 2 else 3
        base <- data.frame(cut = tree$root$cut, side = side)
        if (!tree[[child]]$split) {
          leaves[[as.character(id)]] <- base; lp <- lp + log(1 - rho(1))
        } else {
          lp <- lp + log(rho(1)) - log(length(grid))
          leaves[[as.character(2 * id)]] <- rbind(base, data.frame(cut = tree[[child]]$cut, side = -1))
          leaves[[as.character(2 * id + 1)]] <- rbind(base, data.frame(cut = tree[[child]]$cut, side = 1))
        }
      }
    }
    ids <- names(leaves)
    H <- vapply(leaves, exposure, numeric(1))
    phi <- sapply(leaves, function(path) vapply(x[, 1], function(xi) log_phi(path, xi), numeric(1)))
    phi <- matrix(phi, nrow = n)
    labelings <- as.matrix(expand.grid(rep(list(seq_along(ids)), n)))
    lw <- apply(labelings, 1L, function(z) {
      m <- tabulate(z, nbins = length(ids))
      lp + sum(log_Q(m, H)) + sum(phi[cbind(seq_len(n), z)])
    })
    summaries[[length(summaries) + 1L]] <- data.frame(
      split = tree$root$split, cut = tree$root$cut, nleaf = length(ids), lw = lw
    )
  }
  out <- do.call(rbind, summaries)
  out$p <- exp(out$lw - max(out$lw)); out$p <- out$p / sum(out$p)
  out
}

test_that("soft Particle Gibbs with ancestor sampling is invariant for the exact posterior", {
  skip_on_cran()
  x <- matrix(c(0.12, 0.33, 0.41, 0.58, 0.77, 0.9), ncol = 1)
  region <- matrix(c(0, 1), 1)
  grid <- c(0.3, 0.5, 0.7)
  a <- 1; b <- 0.2; gate <- 6; alpha <- 0.95; eta <- 2

  exact <- soft_exact_posterior(x, grid, region, a, b, gate, alpha, eta)
  p_split <- sum(exact$p[exact$split])
  p_cut <- vapply(grid, function(c) sum(exact$p[exact$split & exact$cut == c]), numeric(1))
  p_nleaf <- vapply(1:4, function(k) sum(exact$p[exact$nleaf == k]), numeric(1))
  expect_equal(sum(exact$p), 1)

  tolerance <- 0.03                   # ~4 batch-means standard errors at 12k correlated draws
  # (resampling after every node | after every level) x (exact | approximate proposal)
  # x (sequential-imputation | auxiliary-rate allocation above exact_max)
  configs <- list(
    list(TRUE, 500L, FALSE, "laplace", 0.5, 0.1),
    list(FALSE, 500L, FALSE, "laplace", 0.5, 0.1),
    list(TRUE, 0L, FALSE, "laplace", 0.5, 0.1),
    list(TRUE, 0L, TRUE, "laplace", 0.5, 0.1),
    list(TRUE, 0L, FALSE, "hard", 0.5, 0.1),
    list(FALSE, 0L, FALSE, "hard", 0.5, 0.1),
    list(TRUE, 0L, TRUE, "hard", 0.5, 0.1),
    list(TRUE, 3L, FALSE, "hard", 0.5, 0.1),
    list(TRUE, 0L, FALSE, "hard", 1, 0.02),
    list(TRUE, 500L, FALSE, "hard", 0.5, 0.1)
  )
  for (config in configs) {
    resample_node <- config[[1]]; exact_max <- config[[2]]
    allocation_rates <- config[[3]]
    set.seed(11)
    raw <- poistree:::SPPT_fit_PGAS(
      x, matrix(0.5, 1, 1), matrix(numeric(), 0, 1), region, list(grid),
      a, b, gate, 1, 1, 0.1, 0, TRUE, alpha, eta, 2L,
      3L, 12000L, 0L, 1L, 0L, FALSE, TRUE, exact_max, 0, resample_node, 1,
      allocation_rates, FALSE, config[[4]], config[[5]], config[[6]]
    )
    states <- raw$state_nodes
    root <- t(vapply(states, function(st) st[st[, 1] == 1, 2:3], numeric(2)))
    chain_split <- mean(root[, 1] >= 0)
    chain_cut <- vapply(grid, function(c) mean(root[, 1] >= 0 & abs(root[, 2] - c) < 1e-12), numeric(1))
    chain_nleaf <- vapply(1:4, function(k) mean(raw$nleaf == k), numeric(1))

    expect_lt(abs(chain_split - p_split), tolerance)
    expect_true(all(abs(chain_cut - p_cut) < tolerance))
    expect_true(all(abs(chain_nleaf - p_nleaf) < tolerance))
    expect_gt(mean(raw$as_rate, na.rm = TRUE), 0.2)
    # Dmax = 2: positions 1, 2, 3 and no event after the last one
    expect_true(all(raw$resampled <= if (resample_node) 2 else 1))
  }
})

test_that("soft Particle Gibbs is reachable through ppt_fit and returns a usable ppt object", {
  set.seed(5)
  x <- matrix(runif(80), ncol = 2)
  colnames(x) <- c("u", "v")
  region <- matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE)
  grid <- matrix(runif(20), ncol = 2)
  test <- matrix(runif(10), ncol = 2)

  fit <- ppt_fit(
    x, region, gating = "soft", scales = "leaf", sampler = "pgas",
    predict_at = grid, test = test,
    particles = 6, iter = 30, burn = 10, max_depth = 3, cut_candidates = 5,
    exact_max = 20, verbose = FALSE
  )
  expect_s3_class(fit, "ppt")
  expect_identical(fit$model$label, "S-PPT")
  expect_identical(fit$model$sampler, "pgas")
  expect_identical(fit$model$algorithm, "PGAS")
  expect_true(fit$control$ancestor_sampling)
  expect_equal(fit$posterior$draws, 20L)
  expect_length(fit$posterior$state$nodes, 20L)
  expect_equal(dim(fit$posterior$state$gate), c(20L, 2L))
  expect_true(all(is.finite(ppt_predict(fit))))
  expect_length(ppt_predict(fit, newdata = grid[1:3, , drop = FALSE]), 3L)
  expect_true(is.finite(ppt_lppd(fit)))
  expect_true(is.finite(fit$diagnostics$particle_ess))
  expect_true(is.finite(fit$diagnostics$ancestor_move_rate))
  # shared-path store: distinct expansions per sweep never exceed the number of
  # particle-node pairs, and at least the root is expanded
  expect_gte(fit$diagnostics$expanded_nodes, 1)
  expect_lte(fit$diagnostics$expanded_nodes, 6 * (2^3 - 1))
  expect_identical(fit$control$resampling_schedule, "tree_node")
  expect_lte(fit$diagnostics$resampling_events, 2^3 - 2)
  expect_output(print(ppt_summary(fit)), "S-PPT")
  marginal <- ppt_marginal(fit, "u", grid = c(0.2, 0.8))
  expect_s3_class(marginal, "ppt_marginal")
  expect_true(all(is.finite(marginal$mean)))
  surface <- ppt_lambda(fit, n = 4)
  expect_equal(nrow(surface), 16L)

  fit_noas <- ppt_fit(
    x, region, gating = "soft", scales = "leaf", sampler = "pgas",
    particles = 4, iter = 12, burn = 2, max_depth = 2, cut_candidates = 4,
    ancestor_sampling = FALSE, resampling = "level", verbose = FALSE
  )
  expect_identical(fit_noas$model$algorithm, "PG")
  expect_identical(fit_noas$control$resampling_schedule, "tree_level")
})

test_that("the forward soft step reproduces the fully adapted increment", {
  # With exact proposals and defensive = 0 the log increment is action and
  # colouring independent, so every particle carries the same weight after the
  # root step and the final ESS of a plain SMC sweep equals the particle count.
  set.seed(3)
  x <- matrix(runif(40), ncol = 2)
  region <- matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE)
  raw <- poistree:::SPPT_fit_PGAS(
    x, matrix(0.5, 1, 2), matrix(numeric(), 0, 2), region, list(c(0.4, 0.6), c(0.5)),
    0.5, 0.5, c(10, 10), c(1, 1), c(1, 1), c(0.1, 0.1), c(0, 0), TRUE, 0.9, 1, 1L,
    8L, 3L, 0L, 1L, 0L, FALSE, TRUE, 500L, 0, TRUE, 1, FALSE, FALSE
  )
  expect_true(all(abs(raw$ess - 8) < 1e-8))
  expect_true(all(raw$resampled == 0))      # one position only: no resampling event
})
