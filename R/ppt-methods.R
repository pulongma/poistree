#' Predict intensity from a fitted PPT
#'
#' Rows of `newdata` that were part of `predict_at` during fitting are served
#' from the stored posterior summaries. Any other rows are evaluated
#' post hoc from the retained posterior state draws via the same machinery as
#' [ppt_lambda()], so prediction no longer requires listing every location
#' in `predict_at` before fitting. (The two paths can differ very slightly in
#' the interval convention: stored summaries follow the fitting backend,
#' post-hoc summaries use the weighted quantile of [ppt_lambda()].)
#' Mean-only predictions avoid calculating posterior quantiles.
#'
#' @param object A fitted `ppt` object.
#' @param newdata Optional matrix of prediction locations. If `NULL`, all
#'   stored predictions are returned.
#' @param type Posterior mean, median, or interval.
#' @param ... Reserved for future prediction backends.
#'
#' @return A numeric vector for `mean` or `median`, or a data frame for
#'   `interval`.
#' @export
ppt_predict <- function(object, newdata = NULL,
                        type = c("median", "mean", "interval"), ...) {
  if (!inherits(object, "ppt")) {
    stop("`object` must inherit from class \"ppt\".", call. = FALSE)
  }
  type <- match.arg(type)
  index <- seq_len(nrow(object$prediction$locations))
  if (!is.null(newdata)) {
    newdata <- .ppt_validate_points(
      newdata, object$data$dimension, object$data$region, "newdata",
      allow_empty = TRUE
    )
    index <- .ppt_match_prediction_rows(
      newdata, object$prediction$locations
    )
    if (anyNA(index)) {
      if (!.ppt_has_state(object)) {
        stop(
          "Some `newdata` rows were not evaluated during fitting and this ",
          "fit carries no posterior state draws. Refit with the current ",
          "poistree version, or include those rows in `predict_at`.",
          call. = FALSE
        )
      }
      if (identical(type, "mean")) {
        lambda <- ppt_lambda(object, at = newdata, type = "draws")
        return(as.numeric(lambda$draws %*% lambda$weights))
      }
      lambda <- ppt_lambda(object, at = newdata, type = "summary")
      if (identical(type, "interval")) {
        return(data.frame(
          median = lambda$median,
          lower = lambda$lower,
          upper = lambda$upper
        ))
      }
      return(lambda[[type]])
    }
  }

  if (identical(type, "interval")) {
    return(data.frame(
      median = object$prediction$median[index],
      lower = object$prediction$lower[index],
      upper = object$prediction$upper[index]
    ))
  }
  object$prediction[[type]][index]
}

#' Posterior marginal intensity for one input
#'
#' For input coordinate \eqn{x_j}, this function computes, for every posterior
#' tree draw,
#' \deqn{m_j(u)=|\mathcal D_{-j}|^{-1}
#'   \int_{\mathcal D_{-j}}\lambda(u,x_{-j})\,d x_{-j}.}
#' By default, the integral is evaluated from every retained posterior state:
#' exactly from axis-aligned boxes for hard models, from exact
#' piecewise-polynomial path integrals for compact soft gates, and from the
#' stable analytic root-scaled logistic path integrals used during fitting
#' (with adaptive quadrature for numerically ill-conditioned paths).
#' Node-scaled logistic gates use deterministic adaptive quadrature for
#' products with different node widths. Supplying
#' `reference` explicitly instead uses a uniform reference-design
#' approximation. Its one-coordinate substitution designs are evaluated from
#' the stored posterior state draws, so `reference` need not have appeared in
#' `predict_at` (rows that were pre-evaluated during fitting are reused when
#' available). Set `average = FALSE` to omit the volume normalization and
#' obtain the projected point-process intensity.
#'
#' @param object A fitted `ppt` object.
#' @param variable One input name or one-based column index.
#' @param grid Optional numeric vector of values for the selected input. The
#'   default is an equally spaced grid over its fitted region.
#' @param n Number of default grid values when `grid` is `NULL`.
#' @param level Pointwise posterior credible level.
#' @param average Divide by the volume of the other input dimensions. The
#'   default preserves the units and overall scale of the fitted intensity.
#' @param reference Optional numeric reference-design matrix with one column
#'   per input. When supplied, its rows approximate the uniform measure on the
#'   fitted region and explicitly select reference-design approximation instead
#'   of state-based integration.
#' @param type Return posterior summaries or the grid-by-draw matrix.
#'
#' @return For `type = "summary"`, a data frame containing the input value,
#'   posterior mean and median, and pointwise credible limits. For
#'   `type = "draws"`, a matrix with grid values in rows and posterior draws
#'   in columns. Draw weights are stored in the `weights` attribute.
#' @export
ppt_marginal <- function(object, variable, grid = NULL, n = 100L,
                         level = 0.95, average = TRUE,
                         reference = NULL,
                         type = c("summary", "draws")) {
  if (!inherits(object, "ppt")) {
    stop("`object` must inherit from class \"ppt\".", call. = FALSE)
  }
  type <- match.arg(type)
  d <- object$data$dimension
  region <- object$data$region
  input_names <- colnames(object$data$x)
  if (is.null(input_names)) input_names <- paste0("x", seq_len(d))

  if (is.character(variable)) {
    if (length(variable) != 1L || is.na(variable) ||
        !variable %in% input_names) {
      stop("`variable` must be one input name or column index.",
           call. = FALSE)
    }
    j <- match(variable, input_names)
  } else {
    if (length(variable) != 1L || !is.finite(variable) ||
        variable != floor(variable) || variable < 1L || variable > d) {
      stop("`variable` must be one input name or column index.",
           call. = FALSE)
    }
    j <- as.integer(variable)
  }

  if (is.null(grid)) {
    if (length(n) != 1L || !is.finite(n) || n != floor(n) || n < 2L) {
      stop("`n` must be an integer of at least 2.", call. = FALSE)
    }
    grid <- seq(region[j, 1L], region[j, 2L], length.out = as.integer(n))
  } else {
    grid <- as.numeric(grid)
    if (!length(grid) || any(!is.finite(grid))) {
      stop("`grid` must be a non-empty finite numeric vector.",
           call. = FALSE)
    }
    tolerance <- sqrt(.Machine$double.eps) *
      max(1, abs(region[j, ]))
    if (any(grid < region[j, 1L] - tolerance |
            grid > region[j, 2L] + tolerance)) {
      stop("Every value in `grid` must lie inside the fitted region.",
           call. = FALSE)
    }
    grid <- pmin(region[j, 2L], pmax(region[j, 1L], grid))
  }
  if (length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.logical(average) || length(average) != 1L || is.na(average)) {
    stop("`average` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  trees <- object$posterior$tree_draws
  exact_hard_leaf <- identical(object$model$gating, "hard") &&
    identical(object$model$scales, "leaf") && length(trees)
  exact_heap_state <- identical(object$posterior$state$mode, "heap") &&
    length(object$posterior$state$nodes)
  marginal_method <- "exact tree integration"

  if (is.null(reference) && exact_hard_leaf) {
    draws <- vapply(
      trees, .ppt_tree_marginal,
      numeric(length(grid)), grid = grid, variable = j,
      region = region, average = isTRUE(average)
    )
    # A scalar grid makes vapply() return one value per tree as a vector.
    # Retain every posterior draw when restoring the grid-by-draw layout.
    draws <- matrix(draws, nrow = length(grid), ncol = length(trees))
  } else if (is.null(reference) && exact_heap_state) {
    draws <- .ppt_state_marginal(
      object, grid = grid, variable = j, average = isTRUE(average)
    )
    marginal_method <- if (identical(object$posterior$state$gate_mode, 3L)) {
      "adaptive posterior-state integration"
    } else {
      "exact posterior-state integration"
    }
  } else {
    if (is.null(reference)) {
      stop(
        "This fit predates the posterior state required for exact marginal ",
        "integration. Supply a uniform `reference` design or refit with the ",
        "current poistree version.",
        call. = FALSE
      )
    }
    reference <- .ppt_validate_points(
      reference, d, region, "reference", allow_empty = FALSE
    )
    evaluation <- do.call(rbind, lapply(grid, function(value) {
      out <- reference
      out[, j] <- value
      out
    }))
    index <- if (is.null(object$prediction$locations)) NA_integer_ else
      .ppt_match_prediction_rows(evaluation, object$prediction$locations)
    if (!anyNA(index) && !is.null(object$prediction$draws) &&
        ncol(as.matrix(object$prediction$draws))) {
      evaluated_draws <- as.matrix(object$prediction$draws)[
        index, , drop = FALSE
      ]
    } else if (.ppt_has_state(object)) {
      evaluated_draws <- .ppt_state_eval(object, evaluation)
    } else {
      stop(
        "The required marginal-design rows were not evaluated during ",
        "fitting and this fit predates the posterior state store. Refit ",
        "with the current poistree version (or include the rows in ",
        "`predict_at`).",
        call. = FALSE
      )
    }
    n_reference <- nrow(reference)
    draws <- vapply(seq_along(grid), function(k) {
      rows <- (k - 1L) * n_reference + seq_len(n_reference)
      colMeans(evaluated_draws[rows, , drop = FALSE])
    }, numeric(ncol(evaluated_draws)))
    # `vapply()` drops its matrix dimension when there is only one posterior
    # draw. Restore the draws-by-grid layout before transposing so single-draw
    # fits follow the same grid-by-draw contract as all other fits.
    draws <- t(matrix(
      draws, nrow = ncol(evaluated_draws), ncol = length(grid)
    ))
    if (!isTRUE(average) && d > 1L) {
      other_axes <- setdiff(seq_len(d), j)
      draws <- draws * prod(
        region[other_axes, 2L] - region[other_axes, 1L]
      )
    }
    marginal_method <- "uniform reference-design approximation"
  }
  weights <- .ppt_posterior_weights(object, ncol(draws))
  colnames(draws) <- paste0("draw", seq_len(ncol(draws)))
  rownames(draws) <- format(grid, digits = 10L, trim = TRUE)
  attr(draws, "grid") <- grid
  attr(draws, "variable") <- input_names[j]
  attr(draws, "weights") <- weights
  attr(draws, "average") <- isTRUE(average)
  attr(draws, "method") <- marginal_method
  if (identical(type, "draws")) return(draws)

  alpha <- (1 - level) / 2
  posterior_quantile <- function(probability) {
    apply(draws, 1L, .ppt_weighted_quantile,
          weights = weights, probability = probability)
  }
  out <- data.frame(
    variable = rep(input_names[j], length(grid)),
    value = grid,
    mean = as.numeric(draws %*% weights),
    median = posterior_quantile(0.5),
    lower = posterior_quantile(alpha),
    upper = posterior_quantile(1 - alpha),
    stringsAsFactors = FALSE
  )
  attr(out, "level") <- level
  attr(out, "average") <- isTRUE(average)
  attr(out, "method") <- marginal_method
  class(out) <- c("ppt_marginal", "data.frame")
  out
}

#' @rdname ppt_predict
#' @param newdata Optional prediction matrix.
#' @method predict ppt
#' @export
predict.ppt <- function(object, newdata = NULL, ...) {
  ppt_predict(object, newdata = newdata, ...)
}

#' Summarize a fitted PPT
#'
#' SMC fits report a target normalizer rather than an absolute Bayesian
#' evidence estimate: the tree-prior normalizing constant is not computed,
#' so `posterior$log_evidence` is `NA`. With a proper intensity prior (`b > 0`),
#' `posterior$log_target_normalizer` adds the root Gamma--Poisson log marginal
#' to `posterior$log_relative_normalizer`. With the improper `b = 0` prior,
#' only the relative normalizer is available, under the backend's formal
#' improper-prior leaf-score convention. The summary labels these quantities
#' separately. RJ-MCMC and Particle-Gibbs do not estimate these normalizers.
#'
#' @param object A fitted `ppt` object.
#' @param ... Reserved for future methods.
#' @return An object of class `summary.ppt`.
#' @export
ppt_summary <- function(object, ...) {
  if (!inherits(object, "ppt")) {
    stop("`object` must inherit from class \"ppt\".", call. = FALSE)
  }
  structure(
    list(
      call = object$call,
      model = object$model$label,
      components = .ppt_model_components(object$model),
      n = object$data$n,
      dimension = object$data$dimension,
      prediction_locations = nrow(object$prediction$locations),
      posterior_draws = object$posterior$draws,
      mean_leaves = object$posterior$mean_leaves,
      mean_max_depth = object$posterior$mean_max_depth %||% NA_real_,
      mean_gate = object$posterior$mean_gate,
      mean_log_likelihood = object$posterior$mean_log_likelihood,
      log_evidence = object$posterior$log_evidence %||% NA_real_,
      log_target_normalizer =
        object$posterior$log_target_normalizer %||% NA_real_,
      log_relative_normalizer =
        object$posterior$log_relative_normalizer %||% NA_real_,
      mean_integrated_intensity =
        object$posterior$mean_integrated_intensity %||% NA_real_,
      lppd = object$posterior$lppd,
      acceptance = object$diagnostics$acceptance
    ),
    class = "summary.ppt"
  )
}

#' @rdname ppt_summary
#' @method summary ppt
#' @export
summary.ppt <- function(object, ...) ppt_summary(object, ...)

#' @method print ppt
#' @export
print.ppt <- function(x, ...) {
  cat("<poistree fit> ", x$model$label, "\n", sep = "")
  cat("  configuration : ",
      paste(.ppt_model_components(x$model), collapse = " + "),
      "\n", sep = "")
  if (!is.null(x$model$algorithm)) {
    cat("  algorithm     : ", x$model$algorithm, "\n", sep = "")
  }
  cat("  observations  : ", x$data$n, " in ", x$data$dimension,
      " dimension(s)\n", sep = "")
  cat("  posterior     : ", x$posterior$draws, " draws; mean leaves ",
      format(x$posterior$mean_leaves, digits = 4L), "\n", sep = "")
  invisible(x)
}

#' @method print summary.ppt
#' @export
print.summary.ppt <- function(x, ...) {
  cat("Bayesian Poisson point-process tree\n")
  cat("  Model          :", x$model, "\n")
  cat("  Configuration  :", paste(x$components, collapse = " + "), "\n")
  cat("  Observations   :", x$n, "\n")
  cat("  Dimension      :", x$dimension, "\n")
  cat("  Prediction rows:", x$prediction_locations, "\n")
  cat("  Posterior draws:", x$posterior_draws, "\n")
  cat("  Mean leaves    :", format(x$mean_leaves, digits = 5L), "\n")
  if (is.finite(x$mean_max_depth)) {
    cat("  Mean max depth :", format(x$mean_max_depth, digits = 5L), "\n")
  }
  if (is.finite(x$mean_gate)) {
    cat("  Mean gate      :", format(x$mean_gate, digits = 5L), "\n")
  }
  cat("  Mean logLik    :", format(x$mean_log_likelihood, digits = 7L), "\n")
  if (is.finite(x$log_evidence %||% NA_real_)) {
    cat("  Log evidence   :", format(x$log_evidence, digits = 7L), "\n")
  } else if (is.finite(x$log_target_normalizer %||% NA_real_)) {
    cat("  Log target normalizer:",
        format(x$log_target_normalizer, digits = 7L), "\n")
  } else if (is.finite(x$log_relative_normalizer %||% NA_real_)) {
    cat("  Log relative normalizer:",
        format(x$log_relative_normalizer, digits = 7L), "\n")
  }
  if (is.finite(x$mean_integrated_intensity)) {
    cat("  Mean integral  :",
        format(x$mean_integrated_intensity, digits = 7L), "\n")
  }
  if (is.finite(x$lppd)) {
    cat("  Test lppd      :", format(x$lppd, digits = 7L), "\n")
  }
  if (length(x$acceptance)) {
    cat("  Acceptance:\n")
    print(round(x$acceptance, 3L))
  } else {
    cat("  Acceptance     : not returned by this summary backend\n")
  }
  invisible(x)
}

#' Extract the posterior mean point-process log likelihood
#'
#' This is the posterior mean in-sample log likelihood, not model evidence.
#'
#' @param object A fitted `ppt` object.
#' @param ... Reserved for compatibility.
#' @return An object of class `logLik`.
#' @export
ppt_logLik <- function(object, ...) {
  if (!inherits(object, "ppt")) {
    stop("`object` must inherit from class \"ppt\".", call. = FALSE)
  }
  structure(
    object$posterior$mean_log_likelihood,
    class = "logLik",
    nobs = object$data$n,
    df = NA_integer_
  )
}

#' @rdname ppt_logLik
#' @method logLik ppt
#' @export
logLik.ppt <- function(object, ...) ppt_logLik(object, ...)

#' Test-pattern log predictive density
#'
#' The joint posterior log predictive density of a test point pattern over
#' the fitted observation region,
#' \deqn{\mathrm{lppd}
#'   = \log \sum_s w_s \exp\Big\{\sum_i \log\{r\lambda^{(s)}(t_i)\}
#'     - r\textstyle\int_{\mathcal D}\lambda^{(s)}\Big\},}
#' with uniform draw weights for MCMC fits and particle weights for SMC.
#' The factor \eqn{r} (`scale`) rescales the fitted intensity to the
#' intensity of the test process, e.g. \eqn{r = (1-p)/p} when training and
#' test patterns are obtained by independent \eqn{p}-thinning of one pattern.
#' With `type = "plugin"` the score is the log-likelihood of the test pattern
#' under the posterior mean intensity \eqn{\bar\lambda=\sum_s w_s\lambda^{(s)}},
#' \deqn{\sum_i \log\{r\bar\lambda(t_i)\} - r\textstyle\int_{\mathcal D}\bar\lambda.}
#' Both scores are log densities, so larger values are better.
#' Stored prediction draws are reused at matching test locations; only
#' missing locations are evaluated from the retained posterior states.
#' When `test` is supplied here, the intensities \eqn{\lambda^{(s)}(t_i)}
#' are evaluated post hoc from the retained posterior state draws and the
#' stored per-draw intensity integrals, so the test pattern does NOT have to
#' be supplied at fitting time. When `test` is `NULL`, the stored value from
#' `ppt_fit(..., test = )` is returned if `scale = 1` and
#' `type = "posterior"`; otherwise the stored test pattern is rescored.
#'
#' @param object A fitted `ppt` object.
#' @param test Optional numeric matrix of test points (one column per
#'   input). If `NULL`, the test pattern supplied at fitting time is used.
#' @param scale Positive factor \eqn{r} multiplying the fitted intensity.
#' @param type `"posterior"` (joint posterior predictive density) or
#'   `"plugin"` (log-likelihood under the posterior mean intensity).
#' @return A numeric scalar with attributes `n_test`, `joint`, `scale`,
#'   and `type`.
#' @export
ppt_lppd <- function(object, test = NULL, scale = 1,
                     type = c("posterior", "plugin")) {
  if (!inherits(object, "ppt")) {
    stop("`object` must inherit from class \"ppt\".", call. = FALSE)
  }
  type <- match.arg(type)
  if (length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    stop("`scale` must be a positive number.", call. = FALSE)
  }
  out <- function(value, n_test) {
    structure(value, n_test = n_test, joint = identical(type, "posterior"),
              scale = scale, type = type)
  }
  if (is.null(test)) {
    if (scale == 1 && identical(type, "posterior") &&
        is.finite(object$posterior$lppd)) {
      return(out(object$posterior$lppd, nrow(object$data$test)))
    }
    test <- object$data$test
    if (is.null(test)) {
      stop(
        "No test pattern is stored. Supply `test` here, or refit with a ",
        "non-empty `test` matrix.",
        call. = FALSE
      )
    }
  }
  test <- .ppt_validate_points(
    test, object$data$dimension, object$data$region, "test",
    allow_empty = FALSE
  )
  integral <- object$posterior$integrated_intensity_draws
  if (is.null(integral) || !length(integral)) {
    stop(
      "This fit predates exact-integral storage; refit it with the current ",
      "poistree version.",
      call. = FALSE
    )
  }
  draws <- .ppt_prediction_draws(object, test)
  if (ncol(draws) != length(integral)) {
    stop("Posterior state draws and integral draws are inconsistent.",
         call. = FALSE)
  }
  draws <- pmax(draws, .Machine$double.xmin)
  weights <- .ppt_posterior_weights(object, ncol(draws))
  if (identical(type, "plugin")) {
    lambda_bar <- as.numeric(draws %*% weights)
    value <- sum(log(scale * lambda_bar)) - scale * sum(weights * integral)
    return(out(value, nrow(test)))
  }
  log_predictive_draw <- colSums(log(scale * draws)) - scale * as.numeric(integral)
  log_weighted <- log(weights) + log_predictive_draw
  center <- max(log_weighted)
  out(center + log(sum(exp(log_weighted - center))), nrow(test))
}

#' Extract the exact posterior intensity integral
#'
#' Every posterior tree draw has an analytically evaluated total intensity 
#' \eqn{A=\int_{\mathcal D}\lambda(x)\,dx}. The hard-tree model uses exact box
#' volumes; no numerical grid is used by this accessor.
#'
#' @param object A fitted `ppt` object.
#' @param type Return the posterior mean or all retained draw-specific
#'   integrals.
#' @return A numeric scalar for `mean` or a numeric vector for `draws`.
#' @export
ppt_integral <- function(object, type = c("mean", "draws")) {
  if (!inherits(object, "ppt")) {
    stop("`object` must inherit from class \"ppt\".", call. = FALSE)
  }
  type <- match.arg(type)
  draws <- object$posterior$integrated_intensity_draws
  mean_integral <- object$posterior$mean_integrated_intensity
  if (is.null(draws) || !length(draws) ||
      is.null(mean_integral) || !is.finite(mean_integral)) {
    stop(
      "This fit predates exact-integral storage; refit it with the current ",
      "poistree version.", call. = FALSE
    )
  }
  if (identical(type, "draws")) as.numeric(draws) else as.numeric(mean_integral)
}

#' Plot fitted PPT intensity summaries
#'
#' Plots the posterior intensity from the stored posterior state draws
#' evaluated over an automatic grid (see [ppt_lambda()]), so it works for any
#' fit whether or not `predict_at` was supplied. One selected dimension gives
#' a curve with a shaded pointwise credible band; two selected dimensions
#' give a raster image of the intensity surface. When the fit has more input
#' dimensions than are shown, the remaining coordinates are fixed at the
#' midpoints of their fitted ranges, so the display is a slice of the
#' intensity. Fits from earlier package versions without stored state fall
#' back to displaying the stored `predict_at` summaries. For
#' publication-quality graphics use [ppt_lambda()] with \pkg{ggplot2}.
#'
#' @param x A fitted `ppt` object.
#' @param type Posterior median or mean.
#' @param dims One dimension for a curve or two dimensions for an intensity
#'   image.
#' @param n Grid resolution per shown dimension.
#' @param level Pointwise credible level for the one-dimensional band.
#' @param points Overlay the observed points: a rug in one dimension, white
#'   points on the surface in two dimensions.
#' @param main,xlab,ylab Optional labels.
#' @param ... Additional graphical arguments passed to [graphics::plot()] or
#'   [graphics::image()].
#' @return The fitted object, invisibly.
#' @method plot ppt
#' @export
plot.ppt <- function(x, type = c("median", "mean"), dims = NULL,
                     n = 60L, level = 0.95, points = FALSE,
                     main = NULL, xlab = NULL, ylab = NULL, ...) {
  type <- match.arg(type)
  d <- x$data$dimension
  region <- x$data$region
  input_names <- colnames(x$data$x)
  if (is.null(input_names)) input_names <- paste0("x", seq_len(d))
  if (is.null(dims)) dims <- if (d == 1L) 1L else c(1L, 2L)
  dims <- as.integer(dims)
  if (!length(dims) || length(dims) > 2L || anyNA(dims) ||
      any(dims < 1L | dims > d) || anyDuplicated(dims)) {
    stop("`dims` must select one or two distinct input dimensions.",
         call. = FALSE)
  }
  if (length(n) != 1L || !is.finite(n) || n != floor(n) || n < 2L) {
    stop("`n` must be an integer of at least 2.", call. = FALSE)
  }
  n <- as.integer(n)
  if (!.ppt_has_state(x)) {
    return(.ppt_plot_stored(x, type, dims, main, xlab, ylab, ...))
  }
  if (is.null(main)) {
    main <- paste(x$model$label, type, "intensity")
    if (d > length(dims)) main <- paste0(main, " (midpoint slice)")
  }
  midpoints <- (region[, 1L] + region[, 2L]) / 2

  if (length(dims) == 1L) {
    j <- dims[1L]
    grid <- seq(region[j, 1L], region[j, 2L], length.out = n)
    at <- matrix(rep(midpoints, each = n), nrow = n)
    at[, j] <- grid
    colnames(at) <- input_names
    lambda <- ppt_lambda(x, at = at, level = level)
    if (is.null(xlab)) xlab <- input_names[j]
    if (is.null(ylab)) ylab <- "Intensity"
    graphics::plot(
      grid, lambda[[type]], type = "n",
      ylim = range(0, lambda$upper, na.rm = TRUE),
      main = main, xlab = xlab, ylab = ylab, ...
    )
    graphics::polygon(
      c(grid, rev(grid)), c(lambda$lower, rev(lambda$upper)),
      col = grDevices::adjustcolor("steelblue", 0.25), border = NA
    )
    graphics::lines(grid, lambda[[type]], lwd = 2, col = "steelblue4")
    if (isTRUE(points)) graphics::rug(x$data$x[, j])
  } else {
    j <- dims[1L]
    k <- dims[2L]
    grid_j <- seq(region[j, 1L], region[j, 2L], length.out = n)
    grid_k <- seq(region[k, 1L], region[k, 2L], length.out = n)
    pairs <- as.matrix(expand.grid(grid_j, grid_k))
    at <- matrix(rep(midpoints, each = nrow(pairs)), nrow = nrow(pairs))
    at[, j] <- pairs[, 1L]
    at[, k] <- pairs[, 2L]
    colnames(at) <- input_names
    lambda <- ppt_lambda(x, at = at, level = level)
    z <- matrix(lambda[[type]], n, n)
    if (is.null(xlab)) xlab <- input_names[j]
    if (is.null(ylab)) ylab <- input_names[k]
    graphics::image(
      grid_j, grid_k, z, col = grDevices::hcl.colors(64L, "Viridis"),
      useRaster = TRUE, main = main, xlab = xlab, ylab = ylab, ...
    )
    if (isTRUE(points)) {
      graphics::points(
        x$data$x[, j], x$data$x[, k], pch = 16, cex = 0.4,
        col = grDevices::adjustcolor("white", 0.7)
      )
    }
  }
  invisible(x)
}

# Legacy display for fits from package versions without stored posterior
# state: show the summaries stored at the `predict_at` locations.
#' @keywords internal
.ppt_plot_stored <- function(x, type, dims, main, xlab, ylab, ...) {
  locations <- x$prediction$locations
  if (is.null(locations) || !nrow(locations)) {
    stop(
      "This fit has neither posterior state draws nor stored predictions; ",
      "refit with the current poistree version.", call. = FALSE
    )
  }
  z <- x$prediction[[type]]
  if (is.null(main)) main <- paste(x$model$label, type, "intensity")
  input_names <- colnames(locations)
  if (is.null(input_names)) {
    input_names <- paste0("x", seq_len(ncol(locations)))
  }
  if (length(dims) == 1L) {
    j <- dims[1L]
    ord <- order(locations[, j])
    if (is.null(xlab)) xlab <- input_names[j]
    if (is.null(ylab)) ylab <- "Intensity"
    graphics::plot(
      locations[ord, j], z[ord], type = "l",
      main = main, xlab = xlab, ylab = ylab, ...
    )
  } else {
    if (is.null(xlab)) xlab <- input_names[dims[1L]]
    if (is.null(ylab)) ylab <- input_names[dims[2L]]
    palette <- grDevices::hcl.colors(64L, "Viridis")
    breaks <- seq(min(z), max(z), length.out = length(palette) + 1L)
    if (diff(range(z)) == 0) {
      color <- rep(palette[ceiling(length(palette) / 2)], length(z))
    } else {
      color <- palette[pmax(1L, pmin(length(palette),
                                    findInterval(z, breaks,
                                                 all.inside = TRUE)))]
    }
    graphics::plot(
      locations[, dims[1L]], locations[, dims[2L]],
      col = color, pch = 16, main = main, xlab = xlab, ylab = ylab, ...
    )
  }
  invisible(x)
}
