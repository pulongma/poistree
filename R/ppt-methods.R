#' Predict intensity from a fitted PPT
#'
#' When every requested row was part of `predict_at` during fitting, stored
#' posterior summaries are returned. Otherwise the complete request is
#' summarized from retained draws using [ppt_lambda()], reusing available
#' stored draws and evaluating missing locations from the posterior states.
#' Prediction therefore does not require listing every location in
#' `predict_at` before fitting. The two paths can use different quantile
#' conventions: stored summaries follow the fitting backend, whereas
#' post-hoc summaries use the weighted quantile of [ppt_lambda()]. Their
#' medians and interval endpoints can therefore differ.
#' Mean-only predictions avoid calculating posterior quantiles.
#'
#' @param object An object inheriting from class `ppt`, usually returned by
#'   [ppt_fit()]. Evaluation at locations not stored during fitting requires
#'   retained posterior states.
#' @param newdata A finite numeric matrix with one row per prediction location
#'   and one column per fitted input, in the same column order as the fitting
#'   data. All locations must lie in the fitted observation region. The
#'   default, `NULL`, returns summaries at the stored `predict_at` locations;
#'   it does not generate a grid. A matrix with zero rows is allowed. Supplied
#'   row order and repeated locations are preserved.
#' @param type Character scalar selecting the returned summary:
#'   \describe{
#'     \item{`"median"`}{Posterior median intensity; the default.}
#'     \item{`"mean"`}{Posterior mean intensity.}
#'     \item{`"interval"`}{Posterior median and lower and upper limits of the
#'       pointwise 95 percent credible interval.}
#'   }
#' @param ... For `predict()`, arguments forwarded to `ppt_predict()`, notably
#'   `type`. Additional arguments received directly by `ppt_predict()` are
#'   currently ignored.
#'
#' @return A numeric vector for `mean` or `median`, or a data frame for
#'   `interval`.
#' @md
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
#' `predict_at`. Stored draws are reused when every row of the complete
#' substitution design was evaluated during fitting; otherwise the full
#' design is evaluated from the retained states. Set `average = FALSE` to
#' omit the volume normalization and
#' obtain the projected point-process intensity.
#'
#' @param object An object inheriting from class `ppt`, usually returned by
#'   [ppt_fit()], with retained posterior states or the prediction draws
#'   needed for the supplied `reference` design.
#' @param variable A single input name or one-based integer column index in
#'   `1, ..., d`, where `d` is the fitted input dimension. If the fitting data
#'   have no column names, use the generated names `"x1"`, `"x2"`, and so on.
#' @param grid A nonempty, finite numeric vector of values for the selected
#'   input, all inside its fitted region. The supplied order is preserved.
#'   The default, `NULL`, creates `n` equally spaced values including both
#'   region endpoints.
#' @param n Integer scalar of at least 2, default `100`, specifying the grid
#'   length when `grid = NULL`. Ignored when `grid` is supplied.
#' @param level Finite numeric scalar strictly between 0 and 1, default `0.95`.
#'   The summary contains the weighted posterior quantiles at
#'   `(1 - level) / 2` and `(1 + level) / 2`. This is a pointwise credible
#'   level, not a simultaneous-band level; it is validated but does not affect
#'   `type = "draws"`.
#' @param average Logical scalar, default `TRUE`. If `TRUE`, divide each
#'   integrated intensity by the volume of all input dimensions except
#'   `variable`, preserving the fitted intensity's units. If `FALSE`, return
#'   the projected point-process intensity without this division.
#' @param reference A finite numeric matrix with at least one row and one
#'   column per fitted input, in fitting-data column order and within the
#'   observation region. Its rows should approximate the uniform distribution
#'   over that region. At each grid value, its `variable` column is replaced
#'   and the resulting intensities are averaged over rows. This explicitly
#'   selects a reference-design approximation. The default, `NULL`, instead
#'   integrates the retained tree states using the model-specific method
#'   described above.
#' @param type Character scalar selecting the output:
#'   \describe{
#'     \item{`"summary"`}{Pointwise posterior summaries; the default.}
#'     \item{`"draws"`}{The grid-by-posterior-draw intensity matrix, with
#'       posterior weights in its `weights` attribute.}
#'   }
#'
#' @return For `type = "summary"`, a data frame containing the input value,
#'   posterior mean and median, and pointwise credible limits. For
#'   `type = "draws"`, a matrix with grid values in rows and posterior draws
#'   in columns. Draw weights are stored in the `weights` attribute.
#' @md
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
#' @method predict ppt
#' @export
predict.ppt <- function(object, newdata = NULL, ...) {
  ppt_predict(object, newdata = newdata, ...)
}

#' Assemble a fitted PPT summary
#'
#' Internal helper used by the public `summary()` method.
#'
#' @param object An object inheriting from class `ppt`, usually returned by
#'   [ppt_fit()]. Stored posterior and diagnostic fields supply the summary.
#' @param ... Additional arguments passed from `summary()`; currently ignored.
#' @return A `summary.ppt` list containing model, data, posterior, and sampler
#'   summaries. Missing backend summaries are represented by `NA`.
#' @md
#' @keywords internal
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
#' @param object An object inheriting from class `ppt`, usually returned by
#'   [ppt_fit()].
#' @param ... Additional arguments accepted for compatibility with the S3
#'   `summary()` and `print()` generics; currently ignored.
#' @return `summary()` returns an object of class
#'   `summary.ppt` containing model, data, posterior, and available sampler
#'   summaries. The print method returns its input invisibly.
#' @md
#' @method summary ppt
#' @export
summary.ppt <- function(object, ...) ppt_summary(object, ...)

#' Print a fitted PPT
#'
#' Display the model configuration, sample size, and posterior draw count.
#'
#' @param x An object inheriting from class `ppt`, usually returned by
#'   [ppt_fit()].
#' @param ... Additional arguments accepted for compatibility with [print()];
#'   currently ignored.
#' @return `x`, invisibly.
#' @md
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

#' @rdname summary.ppt
#' @param x An object of class `summary.ppt`, returned by `summary()` applied
#'   to a `ppt` fit.
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
#' @param object An object inheriting from class `ppt`, usually returned by
#'   [ppt_fit()]. The stored posterior mean in-sample log likelihood is used;
#'   no likelihood is reevaluated by this accessor.
#' @param ... Additional arguments; currently ignored.
#' @return A numeric scalar of class `logLik`, with `nobs` equal to the number
#'   of observed points and `df = NA` because no effective parameter count is
#'   estimated.
#' @md
#' @keywords internal
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
#' @param object An object inheriting from class `ppt`, usually returned by
#'   [ppt_fit()]. Rescoring requires retained intensity-integral draws and
#'   either retained states or stored prediction draws at the test locations.
#' @param test A finite numeric matrix with at least one row and one column
#'   per fitted input, in fitting-data column order. All test points must lie
#'   within the fitted observation region. The default, `NULL`, uses the test
#'   pattern supplied to [ppt_fit()]; an error is raised if none is stored.
#' @param scale Positive, finite numeric scalar \eqn{r} multiplying the fitted
#'   intensity and its integral. The default is `1`. For independent
#'   training/test thinning with training probability \eqn{p}, use
#'   `(1 - p) / p` to predict the test-process intensity.
#' @param type Character scalar selecting the log score:
#'   \describe{
#'     \item{`"posterior"`}{Log of the posterior average joint point-pattern
#'       likelihood; the default.}
#'     \item{`"plugin"`}{Point-pattern log likelihood evaluated at the
#'       posterior mean intensity.}
#'   }
#' @return A numeric scalar with attributes `n_test`, `joint`, `scale`,
#'   and `type`.
#' @md
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

#' Extract the posterior intensity integral
#'
#' Every posterior tree draw has a stored total intensity
#' \eqn{A=\int_{\mathcal D}\lambda(x)\,dx}. Ordinary hard-tree fits use exact
#' box volumes. Fits from [ppt_fit_quadrature()] instead store the weighted
#' background sum over the physical domain. This accessor returns the stored
#' integrals without performing additional numerical integration.
#'
#' @param object An object inheriting from class `ppt`, usually returned by
#'   [ppt_fit()], with stored posterior intensity integrals. Fits created
#'   before integral storage was introduced must be refitted.
#' @param type Character scalar selecting the returned integrals:
#'   \describe{
#'     \item{`"mean"`}{Stored posterior mean total intensity; the default.
#'       SMC fits use posterior particle weights.}
#'     \item{`"draws"`}{One total intensity per retained posterior draw, in
#'       the same order as the fitted posterior states.}
#'   }
#' @return A numeric scalar for `mean` or a numeric vector for `draws`.
#' @md
#' @keywords internal
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
#' @param x An object inheriting from class `ppt`, usually returned by
#'   [ppt_fit()].
#' @param type Character scalar selecting the plotted intensity:
#'   \describe{
#'     \item{`"median"`}{Pointwise posterior median; the default.}
#'     \item{`"mean"`}{Pointwise posterior mean.}
#'   }
#' @param dims One or two distinct one-based input column indices, each in
#'   `1, ..., d`, where `d` is the fitted input dimension. One index gives a
#'   curve; two give an intensity image. The default, `NULL`, selects input 1
#'   for a one-dimensional fit and inputs 1 and 2 otherwise. Unselected
#'   coordinates are fixed at their region midpoints when states are available.
#' @param n Integer scalar of at least 2, default `60`, giving the number of
#'   equally spaced grid values per selected dimension. A two-dimensional
#'   image evaluates `n^2` locations. This is validated but not used when
#'   displaying a legacy fit without posterior states.
#' @param level Finite numeric scalar strictly between 0 and 1, default `0.95`.
#'   Sets the equal-tail pointwise credible band for a one-dimensional plot.
#'   With retained states it is validated for both plot types; a
#'   two-dimensional image does not display a band. Ignored for legacy fits
#'   without retained states.
#' @param points Logical scalar, default `FALSE`. If `TRUE`, overlay the
#'   observed coordinates as a rug in one dimension or white points in two
#'   dimensions. Ignored for legacy fits without retained states.
#' @param main Plot title as a character string or expression, or `NULL` (the
#'   default) to construct a title from the model label and summary type,
#'   with a midpoint-slice label when needed.
#' @param xlab Horizontal-axis label as a character string or expression, or
#'   `NULL` (the default) to use the name of the first selected input.
#' @param ylab Vertical-axis label as a character string or expression, or
#'   `NULL` (the default) to use `"Intensity"` for a curve or the second
#'   selected input name for an image.
#' @param ... Additional graphical arguments passed to [graphics::plot()] for
#'   curves and legacy displays, or [graphics::image()] for intensity images.
#'   Do not repeat arguments already set by the method, including `type` and
#'   `ylim` for curves or `col` and `useRaster` for images.
#' @return The fitted object, invisibly.
#' @md
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

# Display stored predictions for fits without retained posterior states.
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
