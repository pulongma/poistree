pcg_test_args <- function() {
  list(
    x = cbind(c(.08, .21, .34, .62, .79, .93),
              c(.12, .71, .25, .86, .43, .61)),
    region = matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE),
    sampler = "pcg", a = .5, b = .2, max_depth = 2L,
    cut_candidates = 3L, chains = 1L, iter = 18L, burn = 8L,
    thin = 2L, tree_moves = 1L, change_moves = 1L,
    seed = 112L, verbose = FALSE
  )
}

pcg_test_fit <- function(...) {
  do.call(ppt_fit, utils::modifyList(pcg_test_args(), list(...)))
}

test_that("PCG has a soft-only interface and validates RAM controls", {
  args <- pcg_test_args()
  expect_error(do.call(ppt_fit, c(args, list(gating = "hard"))),
               "not available")
  validate <- poistree:::.ppt_validate_ram_controls
  expect_identical(validate(.234, .7, NULL, 11L)$adapt, 8L)
  expect_identical(validate(.234, .7, 0, 0L)$adapt, 0L)
  expect_identical(validate(.234, 1, 8, 8L)$adapt, 8L)
  for (value in list(0, 1, -1, Inf, NA_real_, c(.2, .3), "0.2", 1i)) {
    expect_error(validate(value, .7, 2, 8), "ram_target")
  }
  for (value in list(.5, 1.1, Inf, NA_real_, c(.7, .8), "0.7", 1i)) {
    expect_error(validate(.234, value, 2, 8), "ram_decay")
  }
  for (value in list(-1, 1.5, 9, Inf, NA_real_, c(1, 2), "2", 1i)) {
    expect_error(validate(.234, .7, value, 8), "ram_adapt")
  }
  expect_error(pcg_test_fit(ram_adapt = 9), "ram_adapt")
  expect_error(pcg_test_fit(sd_gate = 0), "sd_gate")
  expect_error(pcg_test_fit(gate_structure = "shared", sd_gate = c(.1, .2)),
               "Shared gates require scalar")
})

test_that("PCG fits logistic and compact gates with both gate structures", {
  for (family in c("logistic", "compact")) {
    for (structure in c("dimension", "shared")) {
      fit <- pcg_test_fit(gate_family = family, gate_structure = structure,
                          chains = 2L, ram_adapt = 5L)
      expect_s3_class(fit, "ppt")
      expect_identical(fit$model$sampler, "pcg")
      expect_identical(fit$model$algorithm, "Partially collapsed Gibbs (RAM)")
      expect_identical(fit$prior$gate$family, family)
      expect_identical(fit$prior$gate$structure, structure)
      expect_identical(fit$control$ram_adapt, 5L)
      expect_identical(fit$control$cut_candidates, 3L)
      expect_length(fit$posterior$particle_weights, 0L)
      expect_equal(fit$posterior$draws, 10L)
      expect_true(all(is.finite(fit$prediction$draws)))
      expect_true(all(fit$prediction$draws > 0))
      ram <- fit$diagnostics$ram
      dimension <- if (structure == "shared") 1L else 2L
      expect_length(ram$covariance, 2L)
      expect_length(ram$factor, 2L)
      expect_equal(ram$updates, rep(5L, 2L))
      expect_equal(ram$failures, rep(0L, 2L))
      for (chain in seq_len(2L)) {
        expect_equal(dim(ram$covariance[[chain]]), rep(dimension, 2L))
        expect_equal(ram$covariance[[chain]],
                     tcrossprod(ram$factor[[chain]]), tolerance = 1e-12)
        expect_true(all(eigen(ram$covariance[[chain]], symmetric = TRUE,
                              only.values = TRUE)$values > 0))
      }
      if (structure == "shared") {
        expect_equal(fit$posterior$state$gate[, 1L],
                     fit$posterior$state$gate[, 2L])
      }
      expect_true(fit$diagnostics$gate_joint_acceptance >= 0)
      expect_true(fit$diagnostics$gate_joint_acceptance <= 1)
      expect_equal(unname(fit$diagnostics$gate_acceptance),
                   rep(fit$diagnostics$gate_joint_acceptance, 2L))
      expect_length(fit$diagnostics$chain_gate_joint_acceptance, 2L)
      evaluated <- ppt_lambda(fit, at = fit$data$x, type = "draws")
      expect_equal(evaluated$draws, fit$prediction$draws, tolerance = 1e-11)
      expect_equal(ppt_predict(fit, type = "mean"), rowMeans(fit$prediction$draws))
      expect_true(all(is.finite(poistree:::ppt_integral(fit, "draws"))))
      expect_s3_class(poistree:::ppt_logLik(fit), "logLik")
      expect_s3_class(summary(fit), "summary.ppt")
      diagnostics <- ppt_diagnostics(fit)
      expect_identical(diagnostics$sampler, "pcg")
      expect_identical(diagnostics$ram, ram)
    }
  }
})

test_that("PCG respects fixed gates and nonadaptive proposal scales", {
  for (family in c("logistic", "compact")) {
    fixed <- pcg_test_fit(gate_family = family, update_gate = FALSE,
                          gate = c(7, 11), sd_gate = c(.04, .09))
    expect_equal(fixed$posterior$state$gate,
                 matrix(rep(c(7, 11), each = fixed$posterior$draws), ncol = 2L))
    expect_identical(fixed$diagnostics$ram$updates, 0L)
    expect_equal(fixed$diagnostics$ram$covariance[[1L]], diag(c(.04, .09)^2))
    nonadaptive <- pcg_test_fit(gate_family = family, ram_adapt = 0L,
                               sd_gate = c(.04, .09))
    expect_identical(nonadaptive$diagnostics$ram$updates, 0L)
    expect_equal(nonadaptive$diagnostics$ram$covariance[[1L]],
                 diag(c(.04, .09)^2))
  }
})

test_that("PCG is reproducible and freezes RAM before retained draws", {
  first <- pcg_test_fit()
  repeat_fit <- pcg_test_fit()
  longer <- pcg_test_fit(iter = 24L)
  expect_identical(first$prediction, repeat_fit$prediction)
  expect_identical(first$posterior, repeat_fit$posterior)
  expect_identical(first$diagnostics, repeat_fit$diagnostics)
  expect_identical(first$control$ram_adapt, 6L)
  expect_identical(first$diagnostics$ram, longer$diagnostics$ram)
  root <- pcg_test_fit(max_depth = 0L, burn = 0L, iter = 4L,
                       ram_adapt = 0L, tree_moves = 0L, change_moves = 0L)
  expect_true(all(root$diagnostics$leaf_count_trace == 1))
  expect_true(all(root$diagnostics$max_depth_trace == 0))
  expect_identical(root$diagnostics$ram$updates, 0L)
})

test_that("RAM rank-one updates match their matrix formula", {
  for (dimension in c(1L, 3L)) {
    L <- diag(seq(.07, .13, length.out = dimension), nrow = dimension)
    if (dimension == 3L) L[3, 1] <- .02
    u <- seq_len(dimension) * c(-1, 1, -1)[seq_len(dimension)]
    for (acceptance in c(0, .234, 1)) {
      target <- .234
      eta <- min(1, dimension * 12^(-.7))
      expected <- L %*% (diag(dimension) + eta * (acceptance - target) *
                          tcrossprod(u) / sum(u^2)) %*% t(L)
      result <- poistree:::ppstree_ram_inspect(
        L, u, acceptance, target, .7, 12L
      )
      expect_true(result$success)
      expect_equal(result$eta, eta, tolerance = 1e-14)
      expect_equal(result$covariance, expected, tolerance = 1e-12)
      expect_equal(tcrossprod(result$factor), expected, tolerance = 1e-12)
      expect_true(all(eigen(result$covariance, symmetric = TRUE,
                            only.values = TRUE)$values > 0))
    }
  }
})

test_that("PCG gate targets equal the finite sum over labeled likelihoods", {
  x <- cbind(c(.1, .3, .65, .9), c(.2, .9, .5, .95))
  region <- matrix(c(0, 1, 0, 1), ncol = 2L, byrow = TRUE)
  splits <- rbind(c(1, 0, .5), c(2, 0, .25), c(3, 1, .7))
  lambda <- c(.3, .7, 1.1, 1.8)
  assignments <- as.matrix(expand.grid(rep(list(seq_len(4L)), nrow(x))))
  for (family in c("logistic", "compact")) {
    for (shared in c(FALSE, TRUE)) {
      gate <- if (shared) c(6, 6) else c(5, 8)
      shape <- if (shared) c(2, 2) else c(2, 4)
      rate <- if (shared) c(.2, .2) else c(.2, .4)
      right <- function(value, cut, width, slope) {
        if (family == "logistic") return(plogis(slope * (value - cut)))
        t <- pmax(0, pmin(1, (value - cut + width / slope) /
                           (2 * width / slope)))
        3 * t^2 - 2 * t^3
      }
      root_right <- function(value) right(value, .5, 1, gate[1L])
      left_right <- function(value) right(value, .25, .5, gate[1L])
      right_right <- function(value) right(value, .7, 1, gate[2L])
      phi <- cbind(
        (1 - root_right(x[, 1L])) * (1 - left_right(x[, 1L])),
        (1 - root_right(x[, 1L])) * left_right(x[, 1L]),
        root_right(x[, 1L]) * (1 - right_right(x[, 2L])),
        root_right(x[, 1L]) * right_right(x[, 2L])
      )
      integral <- function(f) integrate(f, 0, 1, rel.tol = 1e-11,
                                        subdivisions = 1000L)$value
      exposure <- c(
        integral(function(t) (1 - root_right(t)) * (1 - left_right(t))),
        integral(function(t) (1 - root_right(t)) * left_right(t)),
        integral(root_right) * integral(function(t) 1 - right_right(t)),
        integral(root_right) * integral(right_right)
      )
      prior <- if (shared) {
        (shape[1L] - 1) * log(gate[1L]) - rate[1L] * gate[1L]
      } else sum((shape - 1) * log(gate) - rate * gate)
      weighted <- sweep(phi, 2L, lambda, "*")
      labeled_weights <- apply(assignments, 1L, function(z) {
        prod(weighted[cbind(seq_len(nrow(x)), z)])
      })
      expected <- prior - sum(lambda * exposure) + log(sum(labeled_weights))
      result <- poistree:::ppstree_pcg_inspect(
        x, region, splits, gate, lambda, shape, rate, c(0, 0),
        as.integer(shared), as.integer(family == "compact")
      )
      expect_equal(as.integer(result$leaf_ids), 4:7)
      expect_equal(result$phi, phi, tolerance = 1e-12)
      expect_equal(result$exposure, exposure, tolerance = 1e-8,
                   ignore_attr = TRUE)
      expect_equal(result$log_target, expected, tolerance = 1e-8)
      jacobian <- if (shared) log(gate[1L]) else sum(log(gate))
      expect_equal(result$log_scale_target, expected + jacobian, tolerance = 1e-8)
      expected_prob <- weighted / rowSums(weighted)
      expect_equal(result$allocation_prob, expected_prob, tolerance = 1e-12)
      independent_prob <- apply(assignments, 1L, function(z) {
        prod(result$allocation_prob[cbind(seq_len(nrow(x)), z)])
      })
      expect_equal(independent_prob, labeled_weights / sum(labeled_weights),
                   tolerance = 1e-12)
      expect_equal(sum(independent_prob), 1, tolerance = 1e-12)
    }
  }
})
