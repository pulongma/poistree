test_that("geometry cache control is a logical backend option", {
  for (sampler in c("rjmcmc", "irjmcmc", "pcg")) {
    args <- geometry_cache_args()
    args$sampler <- sampler
    for (invalid in list(NULL, NA, 0, 1, "TRUE", logical(), c(TRUE, FALSE))) {

      args["cache_geometry"] <- list(invalid)
      expect_error(do.call(ppt_fit, args), "cache_geometry")
    }
    args$cache_geometry <- TRUE
    enabled <- do.call(ppt_fit, args)
    args$cache_geometry <- NULL
    default <- do.call(ppt_fit, args)
    expect_identical(enabled$prediction, default$prediction)
    expect_identical(enabled$posterior, default$posterior)
    expect_identical(default$control$cache_geometry, TRUE)
  }
})

test_that("cached and reference samplers agree through accepted and rejected moves", {
  topology_changed <- FALSE
  tree_rejection_seen <- FALSE
  gate_changed <- FALSE
  for (sampler in c("rjmcmc", "irjmcmc", "pcg")) {
    for (scale in c("root", "node")) {
      for (structure in c("dimension", "shared")) {
        for (update in c(FALSE, TRUE)) {
          args <- geometry_cache_args()
          args$sampler <- sampler
          args$gate_scale <- scale
          args$gate_structure <- structure
          args$update_gate <- update
          cached <- do.call(ppt_fit, c(args, list(cache_geometry = TRUE)))
          cached_rng <- .Random.seed
          reference <- do.call(ppt_fit, c(args, list(cache_geometry = FALSE)))
          expect_cache_same_fit(cached, reference, cached_rng, .Random.seed)
          expect_identical(reference$control$cache_geometry, FALSE)
          topology_changed <- topology_changed ||
            length(unique(cached$diagnostics$leaf_count_trace)) > 1L
          acceptance <- cached$diagnostics$tree_acceptance
          tree_rejection_seen <- tree_rejection_seen ||
            any(is.finite(acceptance) & acceptance > 0 & acceptance < 1)
          gate_changed <- gate_changed ||
            any(cached$posterior$state$gate != args$gate)
        }
      }
    }
  }
  expect_true(topology_changed)
  expect_true(tree_rejection_seen)
  expect_true(gate_changed)
})

test_that("compact gates and root-only empty prediction cases remain compatible", {
  for (sampler in c("rjmcmc", "irjmcmc", "pcg")) {
    for (root_only in c(FALSE, TRUE)) {
      args <- geometry_cache_args()
      args$sampler <- sampler
      args$gate_family <- "compact"
      if (root_only) {
        args$max_depth <- 0L
        args$tree_moves <- args$change_moves <- 0L
        args$chains <- 2L
        args$predict_at <- matrix(numeric(), 0L, 2L)
      }
      cached <- do.call(ppt_fit, c(args, list(cache_geometry = TRUE)))
      cached_rng <- .Random.seed
      reference <- do.call(ppt_fit, c(args, list(cache_geometry = FALSE)))
      expect_cache_same_fit(cached, reference, cached_rng, .Random.seed)
      expect_identical(cached$model$gate_scale, "node")
    }
  }
})

test_that("shared ancestor calculations agree on repeated dimensions and steep gates", {
  region <- matrix(c(-2, 3, 4, 7), 2L, byrow = TRUE)
  x <- cbind(c(-2, -1.1, -.2, .8, 2.4, 3), c(4, 5.7, 6.1, 4.9, 6.8, 7))
  splits <- rbind(c(1, 0, 1), c(2, 0, -.4), c(4, 1, 5.9),
                  c(8, 0, -1.2), c(3, 1, 6.2))
  for (family in 0:2) {
    for (gate in list(c(1.5, 2.1), c(400, 700))) {
      for (shared in c(FALSE, TRUE)) {
        if (shared) gate[] <- gate[1L]
        args <- list(X = x, region = region, splits = splits, gate = gate,
                     lambda = c(.2, .4, .8, 1.6, 3.2, 6.4),
                     a_gate = c(2, 2), b_gate = c(.1, .1), gate_min = c(0, 0),
                     gate_shared = as.integer(shared), gate_family = family)
        cached <- do.call(poistree:::ppstree_pcg_inspect,
                           c(args, list(cache_geometry = TRUE)))
        reference <- do.call(poistree:::ppstree_pcg_inspect,
                              c(args, list(cache_geometry = FALSE)))
        expect_equal(cached, reference, tolerance = 1e-10)
        expect_equal(rowSums(cached$phi), rep(1, nrow(x)), tolerance = 1e-12)
        expect_equal(sum(cached$exposure), prod(region[, 2] - region[, 1]),
                     tolerance = 1e-7)
      }
    }
  }
})

test_that("diagnostic sampler caches match uncached branching chains", {
  data <- geometry_cache_args()
  topology_changed <- FALSE
  gate_changed <- FALSE
  tree_rejection_seen <- FALSE
  for (family in 0:2) {
    for (sampler in c("rjmcmc", "irjmcmc", "pcg")) {
      args <- list(
        X = data$x, mon = data$predict_at, region = data$region,
        a = data$a, b = data$b, gate = data$gate, a_gate = data$a_gate,
        b_gate = data$b_gate, sd_gate = data$sd_gate, gate_min = 0,
        gate_shared = as.integer(family == 1L), alpha = data$alpha,
        eta = data$eta, Dmax = data$max_depth, nmin = 1L,
        iters = data$iter, burn = data$burn, thin = data$thin,
        nmove = data$tree_moves, ncc = data$change_moves,
        cut_mode = 1L, ncand = data$cut_candidates, update_gate = 1L,
        gate_family = family, informed = sampler == "irjmcmc",
        pcg = sampler == "pcg", ram_adapt = 12L, verbose = FALSE
      )
      set.seed(data$seed)
      cached <- do.call(poistree:::ppstree_diag,
                         c(args, list(cache_geometry = TRUE)))
      cached_rng <- .Random.seed
      set.seed(data$seed)
      reference <- do.call(poistree:::ppstree_diag,
                            c(args, list(cache_geometry = FALSE)))
      expect_identical(cached_rng, .Random.seed)
      expect_equal(cached, reference, tolerance = 1e-10)
      expect_identical(cached$nleaf, reference$nleaf)
      topology_changed <- topology_changed || length(unique(cached$nleaf)) > 1L
      gate_changed <- gate_changed || any(cached$gate != data$gate)
      tree_rejection_seen <- tree_rejection_seen ||
        any(is.finite(cached$accept) & cached$accept > 0 & cached$accept < 1)
    }
  }
  expect_true(topology_changed)
  expect_true(gate_changed)
  expect_true(tree_rejection_seen)
})
