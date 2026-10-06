test_that("cached candidate cuts have the expected tied-data support", {
  x <- cbind(c(0, .2, .2, .4, .4, .6, .6, .8, .8, 1), rep(.5, 10L))
  region <- matrix(c(0, 1, 0, 1), 2L, byrow = TRUE)
  splits <- matrix(c(1, 0, .5), ncol = 3L)
  # These expectations are hand-calculated for the tied observations above.
  # The constant second coordinate has no admissible split in any mode.
  expected <- list(
    list(c(.3, .5, .7), .3, .7),
    list(c(.4, .6, .8), .4, .8),
    list((7:24) / 31, (13:24) / 62, .5 + (7:18) / 62)
  )
  for (mode in 0:2) {
    args <- list(X = x, region = region, splits = splits, Dmax = 2L,
                 nmin = 2L, cut_mode = mode, ncand = 5L)
    cached <- do.call(poistree:::ppstree_cuts_inspect,
                      c(args, list(cache_cuts = TRUE)))
    direct <- do.call(poistree:::ppstree_cuts_inspect,
                      c(args, list(cache_cuts = FALSE)))
    expect_identical(cached[c("node_ids", "cuts", "can_split")],
                     direct[c("node_ids", "cuts", "can_split")])
    expect_equal(cached$cache_entries, 3)
    expect_equal(direct$cache_entries, 0)
    expect_equal(as.integer(cached$node_ids), 1:3)
    expect_identical(as.logical(cached$can_split), c(FALSE, TRUE, TRUE))
    for (node in 1:3) {
      expect_equal(as.numeric(cached$cuts[[node]][[1L]]),
                   expected[[mode + 1L]][[node]], tolerance = 1e-14)
      expect_length(cached$cuts[[node]][[2L]], 0L)
    }
    args$Dmax <- 1L
    shallow <- do.call(poistree:::ppstree_cuts_inspect, args)
    expect_false(any(shallow$can_split))
    args$Dmax <- 2L
    args$nmin <- 3L
    too_small <- do.call(poistree:::ppstree_cuts_inspect, args)
    expect_false(any(too_small$can_split))
    expect_length(too_small$cuts[[2L]][[1L]], 0L)
    expect_length(too_small$cuts[[3L]][[1L]], 0L)
  }

  # Reusing the same inspection entry point with different controls/data
  # must not retain a previous fit's support cache.
  args <- list(X = x, region = region, splits = matrix(numeric(), 0L, 3L),
               Dmax = 1L, nmin = 4L, cut_mode = 0L, ncand = 2L)
  root <- do.call(poistree:::ppstree_cuts_inspect, args)
  expect_equal(as.numeric(root$cuts[[1L]][[1L]]), .5)
  args$X[,] <- .5
  constant <- do.call(poistree:::ppstree_cuts_inspect, args)
  expect_false(any(constant$can_split))
  expect_length(constant$cuts[[1L]][[1L]], 0L)
})

test_that("cut caching preserves complete seeded sampler trajectories", {
  accepted <- rejected <- rep(FALSE, 3L)
  gate_both <- changed_rule <- FALSE
  for (sampler in c("pcg", "rjmcmc", "irjmcmc")) {
    for (family in 0:2) {
      for (mode in 0:2) {
        args <- cut_cache_native_args()
        args$pcg <- sampler == "pcg"
        args$informed <- sampler == "irjmcmc"
        args$gate_family <- family
        args$cut_mode <- mode
        args$gate_shared <- as.integer((family + mode) %% 2L == 0L)
        args$nmin <- if (mode == 2L) 2L else 1L
        # Ties exercise stable candidate ordering and exact boundary routing.
        if (mode == 0L) args$X[2L, ] <- args$X[1L, ]
        cached <- cut_cache_run(args, TRUE)
        direct <- cut_cache_run(args, FALSE)
        expect_identical(cached$rng, direct$rng)
        expect_identical(cached$value, direct$value)
        acceptance <- as.numeric(cached$value$accept)
        accepted <- accepted | acceptance > 0
        rejected <- rejected | acceptance < 1
        gate_both <- gate_both || any(cached$value$gate_accept > 0 &
                                      cached$value$gate_accept < 1)
        changed_rule <- changed_rule ||
          cut_cache_changed_rules(cached$value$state_nodes)
      }
    }
  }
  # Establish that agreement includes accepted/rejected topology proposals,
  # accepted/rejected gate proposals, and changed rules at existing node IDs.
  expect_true(all(accepted))
  expect_true(all(rejected))
  expect_true(gate_both)
  expect_true(changed_rule)
})

test_that("cut caches remain isolated across chains, fits, and direct diagnostics", {
  for (family in 0:2) {
    args <- cut_cache_native_args()
    args$gate_family <- family
    args$chains <- 2L
    args$gate_shared <- as.integer(family == 1L)
    first <- cut_cache_run(args, TRUE)
    expect_identical(first, cut_cache_run(args, FALSE))

    other <- args
    other$X <- 1 - args$X
    other$cut_mode <- 2L
    other$ncand <- 37L
    other$nmin <- 3L
    other$Dmax <- 2L
    expect_identical(cut_cache_run(other, TRUE), cut_cache_run(other, FALSE))
    expect_identical(first, cut_cache_run(args, TRUE))

    for (sampler in c("pcg", "rjmcmc", "irjmcmc")) {
      args$pcg <- sampler == "pcg"
      args$informed <- sampler == "irjmcmc"
      expect_identical(cut_cache_run(args, TRUE, diagnostic = TRUE),
                       cut_cache_run(args, FALSE, diagnostic = TRUE))
    }
    args$Dmax <- 0L
    args$nmove <- args$ncc <- 0L
    args$update_gate <- 0L
    args$grid <- matrix(numeric(), 0L, 2L)
    expect_identical(cut_cache_run(args, TRUE), cut_cache_run(args, FALSE))
  }
})

test_that("reused PCG normalizers equal independent mixture intensities", {
  x <- cbind(c(-2, -1.2, -.3, .9, 2.1, 3), c(4, 6.8, 4.9, 5.6, 6.3, 7))
  region <- matrix(c(-2, 3, 4, 7), 2L, byrow = TRUE)
  splits <- rbind(c(1, 0, .5), c(2, 0, -.75), c(3, 1, 6))
  # A very wide range of rates also checks that normalization stays on the
  # log scale. Compact gates additionally exercise exactly zero memberships.
  rates <- c(1e-100, .2, 3, 1e80)
  for (family in 0:2) {
    for (shared in c(FALSE, TRUE)) {
      gate <- if (shared) c(5, 5) else c(3, 8)
      gate_logs <- function(value, cut, axis, local_width) {
        width <- if (family == 0L) diff(region[axis, ]) else local_width
        if (family == 1L) {
          t <- pmax(0, pmin(1, (value - cut + width / gate[axis]) /
                            (2 * width / gate[axis])))
          right <- t^2 * (3 - 2 * t)
          return(cbind(log1p(-right), log(right)))
        }
        z <- gate[axis] * (value - cut) / width
        cbind(plogis(z, lower.tail = FALSE, log.p = TRUE),
              plogis(z, log.p = TRUE))
      }
      root <- gate_logs(x[, 1L], .5, 1L, 5)
      left <- gate_logs(x[, 1L], -.75, 1L, 2.5)
      right <- gate_logs(x[, 2L], 6, 2L, 3)
      log_phi <- cbind(root[, 1L] + left[, 1L],
                       root[, 1L] + left[, 2L],
                       root[, 2L] + right[, 1L],
                       root[, 2L] + right[, 2L])
      log_weight <- sweep(log_phi, 2L, log(rates), "+")
      expected <- apply(log_weight, 1L, function(row) {
        pivot <- max(row)
        pivot + log(sum(exp(row - pivot)))
      })
      result <- poistree:::ppstree_pcg_inspect(
        x, region, splits, gate, rates, c(2, 2), c(.2, .2), c(0, 0),
        as.integer(shared), family
      )
      expect_equal(as.numeric(result$log_normalizers), expected, tolerance = 1e-12)
      expect_equal(result$allocation_prob, exp(log_weight - expected),
                   tolerance = 1e-12)
      expect_equal(rowSums(result$allocation_prob), rep(1, nrow(x)),
                   tolerance = 1e-12)
      expect_true(all(is.finite(result$log_normalizers)))
    }
  }
})
