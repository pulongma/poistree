membership_with_quadrature <- function(background, weights, region, expression) {
  poistree:::qpp_set_quadrature(background, weights, region)
  on.exit(poistree:::qpp_clear_quadrature())
  force(expression)
}

membership_weighted_data <- function(kind = "matched") {
  args <- geometry_cache_args()
  q <- rbind(args$x, args$test, c(0, 0), c(1, 1), args$x[c(1L, 7L), ])
  args$x <- args$x[c(8, 1, 4, 1, 7, 3, 10, 12, 2, 6, 9, 11), ]
  if (kind == "partial") args$x[c(2L, 9L), ] <- rbind(c(.17, .38), c(.79, .58))
  if (kind == "unmatched") args$x <- .001 + .997 * args$x
  args$background <- q
  args$weights <- seq_len(nrow(q)) / sum(seq_len(nrow(q))) * 2.7
  args$iter <- 44L
  args$burn <- 14L
  # Deliberately distinct from both the training rows and integration rule.
  args$predict_at <- rbind(c(.11, .23), c(.38, .62), c(.57, .43), c(.89, .71))
  args$test <- rbind(c(.13, .28), c(.47, .52), c(.84, .81))
  args
}

test_that("event to background mapping is exact and keeps repeated events", {
  region <- rbind(c(0, 1), c(0, 1))
  q <- rbind(c(0, .25), c(.2, .4), c(.2, .4), c(.5, .75), c(1, 1))
  full <- q[c(4, 2, 3, 1, 5, 2), ]
  near <- c(.2, .4 + .Machine$double.eps)
  signed_zero <- c(-0, .25)
  x <- rbind(full, near, signed_zero, c(.1, .9))
  expected <- c(4L, 2L, 2L, 1L, 5L, 2L, NA_integer_, NA_integer_, NA_integer_)
  for (family in 0:2) {
    mapping <- membership_with_quadrature(q, rep(1, nrow(q)), region,
      poistree:::ppstree_training_background_inspect(x, region, c(5, 8), family))
    expect_identical(mapping, expected)
    reversed <- q[nrow(q):1L, ]
    other <- membership_with_quadrature(reversed, rep(1, nrow(q)), region,
      poistree:::ppstree_training_background_inspect(x, region, c(5, 8), family))
    expect_identical(other, c(2L, 3L, 3L, 5L, 1L, 3L,
                              NA_integer_, NA_integer_, NA_integer_))
    no_rule <- poistree:::ppstree_training_background_inspect(x, region, c(5, 8), family)
    expect_identical(no_rule, rep(NA_integer_, nrow(x)))
  }
})

test_that("reused event bases agree with independent gates and weighted exposure", {
  splits <- rbind(c(1, 0, .5), c(2, 0, .25), c(4, 1, .45),
                  c(8, 0, .1), c(3, 1, .65))
  rates <- c(.2, .4, .8, 1.6, 3.2, 6.4)
  for (kind in c("matched", "partial", "unmatched")) {
    args <- membership_weighted_data(kind)
    for (family in 0:2) {
      for (shared in 0:1) {
        gate <- if (shared) c(5, 5) else c(4, 8)
        native <- list(X = args$x, region = args$region, splits = splits,
                       gate = gate, lambda = rates, a_gate = c(2, 2),
                       b_gate = c(.2, .2), gate_min = c(0, 0),
                       gate_shared = shared, gate_family = family)
        result <- membership_with_quadrature(args$background, args$weights, args$region, {
          cached <- do.call(poistree:::ppstree_pcg_inspect,
                            c(native, list(cache_geometry = TRUE)))
          direct <- do.call(poistree:::ppstree_pcg_inspect,
                            c(native, list(cache_geometry = FALSE)))
          # The older uncached exposure loop already differs at rounding
          # precision; event memberships themselves remain bit-identical.
          expect_identical(cached$phi, direct$phi)
          expect_identical(cached$log_normalizers, direct$log_normalizers)
          expect_identical(cached$allocation_prob, direct$allocation_prob)
          expect_equal(cached, direct, tolerance = 1e-12)
          cached
        })
        train <- membership_tree_log_phi(args$x, args$region, splits, gate, family)
        background <- membership_tree_log_phi(args$background, args$region,
                                               splits, gate, family)
        phi <- exp(train$log_phi)
        exposure <- colSums(exp(background$log_phi) * args$weights)
        weighted <- sweep(phi, 2L, rates, "*")
        expect_equal(as.integer(result$leaf_ids), train$ids)
        expect_equal(result$phi, phi, tolerance = 1e-12)
        expect_equal(as.numeric(result$exposure), exposure, tolerance = 1e-12)
        expect_equal(result$allocation_prob, weighted / rowSums(weighted),
                     tolerance = 1e-12)
        expect_equal(as.numeric(result$log_normalizers), log(rowSums(weighted)),
                     tolerance = 1e-12)
        expect_equal(sum(result$exposure), sum(args$weights), tolerance = 1e-12)
      }
    }
  }
})

test_that("matched event chains preserve RNG and predictions away from the rule", {
  for (kind in c("matched", "partial", "unmatched")) {
    for (family in c("root", "node", "compact")) {
      for (structure in c("shared", "dimension")) {
        args <- membership_weighted_data(kind)
        args$gate_family <- if (family == "compact") "compact" else "logistic"
        args$gate_scale <- if (family == "compact") "node" else family
        args$gate_structure <- structure
        cached <- do.call(ppt_fit_quadrature, c(args, list(cache_geometry = TRUE)))
        rng <- .Random.seed
        direct <- do.call(ppt_fit_quadrature, c(args, list(cache_geometry = FALSE)))
        expect_identical(rng, .Random.seed)
        expect_equal(cached$posterior, direct$posterior, tolerance = 1e-10)
        expect_equal(cached$prediction, direct$prediction, tolerance = 1e-10)
        expect_equal(cached$diagnostics, direct$diagnostics, tolerance = 1e-10)
        expect_identical(cached$diagnostics$leaf_count_trace,
                         direct$diagnostics$leaf_count_trace)
        topology <- function(fit) lapply(fit$posterior$state$nodes,
                                        function(nodes) nodes[, 1:3, drop = FALSE])
        expect_identical(topology(cached), topology(direct))
        at <- rbind(c(.07, .34), c(.67, .47))
        expect_equal(ppt_lambda(cached, at = at, type = "draws"),
                     ppt_lambda(direct, at = at, type = "draws"), tolerance = 1e-10)
      }
    }
  }
})

test_that("matched events also preserve RJMCMC and native diagnostic chains", {
  compare <- function(args, diagnostic = FALSE) {
    cached <- membership_native_run(args, TRUE, diagnostic)
    direct <- membership_native_run(args, FALSE, diagnostic)
    expect_identical(cached$rng, direct$rng)
    expect_equal(cached$value, direct$value, tolerance = 1e-10)
    expect_identical(cached$value$accept, direct$value$accept)
    if (diagnostic) expect_identical(cached$value$nleaf, direct$value$nleaf)
    else {
      topology <- function(result) lapply(result$value$state_nodes,
                                          function(nodes) nodes[, 1:3, drop = FALSE])
      expect_identical(topology(cached), topology(direct))
    }
  }
  for (family in 0:2) {
    for (shared in 0:1) {
      args <- cut_cache_native_args()
      data <- membership_weighted_data("matched")
      args$X <- data$x
      args$gate_family <- family
      args$gate_shared <- shared
      args$pcg <- FALSE
      membership_with_quadrature(data$background, data$weights, data$region, {
        compare(args)
        compare(args, diagnostic = TRUE)
        args$pcg <- TRUE
        compare(args, diagnostic = TRUE)
      })
    }
  }
})
