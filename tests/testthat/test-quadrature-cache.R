# Independent R evaluation on the retained tree. This deliberately traverses
# each serialized draw directly rather than calling the C++ state evaluator.
quadrature_cache_r_intensity <- function(fit, at) {
  evaluate_draw <- function(draw) {
    nodes <- fit$posterior$state$nodes[[draw]]
    gate <- fit$posterior$state$gate[draw, ]
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
  vapply(seq_len(fit$posterior$draws), evaluate_draw, numeric(nrow(at)))
}

quadrature_cache_args <- function() {
  args <- geometry_cache_args()
  args$background <- cbind(c(0, .04, .16, .23, .31, .42, .53, .64, .76, .88, 1),
                           c(1, .12, .73, .36, .84, .28, .55, .06, .66, .42, 0))
  args$weights <- c(.01, .04, .02, .15, .07, .13, .09, .18, .11, .08, .12) * 2.7
  # Prediction and integration use exactly the same locations in California.
  args$predict_at <- args$background
  args
}

test_that("weighted cached fitting matches the reference and independent R exposure", {
  for (scale in c("root", "node")) {
    for (structure in c("dimension", "shared")) {
      for (update in c(FALSE, TRUE)) {
        args <- quadrature_cache_args()
        args$gate_scale <- scale
        args$gate_structure <- structure
        args$update_gate <- update
        cached <- do.call(ppt_fit_quadrature, c(args, list(cache_geometry = TRUE)))
        cached_rng <- .Random.seed
        reference <- do.call(ppt_fit_quadrature, c(args, list(cache_geometry = FALSE)))
        expect_cache_same_fit(cached, reference, cached_rng, .Random.seed)
        at_background <- quadrature_cache_r_intensity(cached, args$background)
        integral <- colSums(at_background * args$weights)
        expect_equal(at_background, cached$prediction$draws, tolerance = 1e-10)
        expect_equal(integral, cached$posterior$integrated_intensity_draws,
                     tolerance = 1e-10)
        train_intensity <- quadrature_cache_r_intensity(cached, args$x)
        expected_ll <- colSums(log(train_intensity)) - integral
        expect_equal(mean(expected_ll), cached$posterior$mean_log_likelihood,
                     tolerance = 1e-10)
        test_intensity <- quadrature_cache_r_intensity(cached, args$test)
        test_ll <- colSums(log(test_intensity)) - integral
        expected_lppd <- max(test_ll) + log(mean(exp(test_ll - max(test_ll))))
        expect_equal(expected_lppd, cached$posterior$lppd, tolerance = 1e-10)
      }
    }
  }
})

test_that("exposure caches cannot leak across quadrature rules, ordinary fits, or errors", {
  for (scale in c("root", "node")) {
    plain_args <- geometry_cache_args()
    plain_args$gate_scale <- scale
    plain_before <- do.call(ppt_fit, plain_args)
    args <- quadrature_cache_args()
    args$gate_scale <- scale
    first <- do.call(ppt_fit_quadrature, args)
    rng_before <- .Random.seed
    other_args <- args
    other_args$weights <- rev(args$weights) * 3
    other_args$background <- 1 - args$background
    other <- do.call(ppt_fit_quadrature, other_args)
    expect_false(isTRUE(all.equal(first$posterior$integrated_intensity_draws,
                                 other$posterior$integrated_intensity_draws)))
    expect_error(do.call(ppt_fit_quadrature, c(args, list(cache_geometry = NA))),
                 "cache_geometry")
    again <- do.call(ppt_fit_quadrature, args)
    expect_identical(first$prediction, again$prediction)
    expect_identical(first$posterior, again$posterior)
    expect_identical(rng_before, .Random.seed)
    plain_after <- do.call(ppt_fit, plain_args)
    expect_identical(plain_before$posterior, plain_after$posterior)
    expect_identical(plain_before$prediction, plain_after$prediction)
  }
})
