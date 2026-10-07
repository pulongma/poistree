region_marginal <- matrix(c(0, 2, -1, 2), ncol = 2, byrow = TRUE)

small_heap_state <- function(rates) {
  stopifnot(length(rates) == 5L)
  cbind(
    heap_id = c(1, 2, 3, 4, 5),
    axis = c(0, 1, -1, -1, -1),
    cut = c(0.8, 0.2, NA, NA, NA),
    lambda = rates,
    xi = rep(NA_real_, 5),
    m = rep(0, 5)
  )
}

make_heap_marginal_fit <- function(gating, scales, rates, gate_mode,
                                   gate = matrix(numeric(), 0, 0),
                                   gate_depth = 0) {
  nodes <- if (is.list(rates)) {
    lapply(rates, small_heap_state)
  } else {
    list(small_heap_state(rates))
  }
  x <- matrix(c(0.1, -0.8, 1.8, 1.8), ncol = 2, byrow = TRUE)
  colnames(x) <- c("x", "y")
  structure(
    list(
      model = list(gating = gating, scales = scales),
      data = list(x = x, dimension = 2L, region = region_marginal),
      prediction = list(
        locations = matrix(numeric(), 0L, 2L),
        draws = matrix(numeric(), 0L, length(nodes))
      ),
      posterior = list(
        tree_draws = list(), particle_weights = numeric(),
        state = list(
          mode = "heap", nodes = nodes, gate = gate,
          gate_mode = gate_mode, gate_depth = gate_depth
        )
      )
    ),
    class = "ppt"
  )
}

make_hard_leaf_marginal_fit <- function() {
  leaf <- function(xlo, xhi, ylo, yhi, lambda) {
    list(
      is_leaf = TRUE,
      region = matrix(c(xlo, xhi, ylo, yhi), ncol = 2, byrow = TRUE),
      lambda = lambda
    )
  }
  tree <- list(
    leaf(0, 0.8, -1, 0.2, 2),
    leaf(0, 0.8, 0.2, 2, 5),
    leaf(0.8, 2, -1, 2, 7)
  )
  x <- matrix(c(0.1, -0.8, 1.8, 1.8), ncol = 2, byrow = TRUE)
  colnames(x) <- c("x", "y")
  structure(
    list(
      model = list(gating = "hard", scales = "leaf"),
      data = list(x = x, dimension = 2L, region = region_marginal),
      posterior = list(
        tree_draws = list(tree), particle_weights = 1,
        state = list(mode = "leafbox")
      )
    ),
    class = "ppt"
  )
}

integrate_state_over_y <- function(x, fit) {
  stats::integrate(
    function(y) {
      at <- cbind(x = rep(x, length(y)), y = y)
      as.numeric(poistree:::.ppt_state_eval(fit, at)[, 1L])
    },
    lower = region_marginal[2, 1],
    upper = region_marginal[2, 2],
    rel.tol = 1e-11, abs.tol = 1e-12,
    subdivisions = 2000L
  )$value / diff(region_marginal[2, ])
}

compact_right <- function(x, cut, width, gate) {
  h <- width / gate
  t <- (x - (cut - h)) / (2 * h)
  ifelse(t <= 0, 0, ifelse(t >= 1, 1, t^2 * (3 - 2 * t)))
}

log1pexp <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

logistic_right_integral <- function(cut, lower, upper, gate) {
  width <- upper - lower
  width / gate * (
    log1pexp(gate * (upper - cut) / width) -
      log1pexp(gate * (lower - cut) / width)
  )
}

test_that("hard terminal-leaf PPT keeps exact box marginalization", {
  fit <- make_hard_leaf_marginal_fit()
  grid <- c(0, 0.4, 0.8, 1.5, 2)
  draws <- ppt_marginal(fit, "x", grid = grid, type = "draws")

  left_average <- (2 * 1.2 + 5 * 1.8) / 3
  expect_equal(as.numeric(draws[, 1]), c(left_average, left_average, 7, 7, 7))
  expect_identical(attr(draws, "method"), "exact tree integration")

  projected <- ppt_marginal(
    fit, "x", grid = grid, average = FALSE, type = "draws"
  )
  expect_equal(as.numeric(projected), as.numeric(draws) * 3)
})

test_that("logistic S-PPT uses draw-specific analytic path integrals", {
  rates <- c(0, 0, 7, 2, 5)
  gates <- rbind(c(8, 5), c(3, 11))
  fit <- make_heap_marginal_fit(
    "soft", "leaf", rates = list(rates, rates), gate_mode = 2L,
    gate = gates
  )
  grid <- c(0, 0.35, 0.8, 1.2, 2)
  draws <- ppt_marginal(fit, "x", grid = grid, type = "draws")

  expected <- vapply(seq_len(nrow(gates)), function(s) {
    right_x <- stats::plogis(gates[s, 1] * (grid - 0.8) / 2)
    right_y <- logistic_right_integral(0.2, -1, 2, gates[s, 2])
    left_y <- 3 - right_y
    (1 - right_x) * (2 * left_y + 5 * right_y) / 3 + 7 * right_x
  }, numeric(length(grid)))

  expect_equal(dim(draws), dim(expected))
  expect_equal(as.numeric(draws), as.numeric(expected), tolerance = 2e-10)
  expect_false(isTRUE(all.equal(draws[, 1], draws[, 2])))
  numerical <- vapply(grid, integrate_state_over_y, numeric(1), fit = fit)
  expect_equal(as.numeric(draws[, 1]), numerical, tolerance = 2e-9)
})

test_that("compact S-PPT agrees with numerical domain integration", {
  fit <- make_heap_marginal_fit(
    "soft", "leaf", rates = c(0, 0, 7, 2, 5), gate_mode = 1L,
    gate = matrix(c(5, 7), nrow = 1L)
  )
  grid <- c(0, 0.35, 0.8, 1.2, 2)
  draws <- ppt_marginal(fit, "x", grid = grid, type = "draws")
  numerical <- vapply(grid, integrate_state_over_y, numeric(1), fit = fit)
  expect_equal(as.numeric(draws[, 1]), numerical, tolerance = 2e-9)
  expect_identical(
    attr(draws, "method"), "exact posterior-state integration"
  )
})
