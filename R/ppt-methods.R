#' Predict intensity from a fitted PPT
#'
#' All current backends store posterior intensity summaries at the locations
#' supplied through `predict_at` during fitting. New prediction rows must
#' therefore be among those stored locations.
#'
#' @param object A fitted `ppt` object.
#' @param newdata Optional matrix selecting rows from the stored prediction
#'   locations. If `NULL`, all stored predictions are returned.
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
      stop(
        "Some `newdata` rows were not evaluated during fitting. Refit with ",
        "those rows included in `predict_at`.",
        call. = FALSE
      )
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

#' @rdname ppt_predict
#' @param newdata Optional prediction matrix.
#' @method predict ppt
#' @export
predict.ppt <- function(object, newdata = NULL, ...) {
  ppt_predict(object, newdata = newdata, ...)
}

#' Summarize a fitted PPT
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
      kappa = object$posterior$kappa %||% NA_real_,
      tau = object$posterior$tau %||% NA_real_,
      mean_gate = object$posterior$mean_gate,
      mean_log_likelihood = object$posterior$mean_log_likelihood,
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
  if (is.finite(x$kappa)) {
    cat("  Mean kappa     :", format(x$kappa, digits = 5L), "\n")
  }
  if (is.finite(x$tau)) {
    cat("  Mean tau       :", format(x$tau, digits = 5L), "\n")
  }
  if (is.finite(x$mean_gate)) {
    cat("  Mean gate      :", format(x$mean_gate, digits = 5L), "\n")
  }
  cat("  Mean logLik    :", format(x$mean_log_likelihood, digits = 7L), "\n")
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

#' Extract model evidence
#'
#' SMC fits return their log marginal-likelihood estimate. RJ-MCMC and PGAS do
#' not currently estimate the marginal likelihood and return `NA` with an
#' explanatory attribute.
#'
#' @param object A fitted `ppt` object.
#' @param warn Emit a warning when evidence is unavailable.
#' @return A numeric scalar.
#' @export
ppt_evidence <- function(object, warn = TRUE) {
  if (!inherits(object, "ppt")) {
    stop("`object` must inherit from class \"ppt\".", call. = FALSE)
  }
  log_evidence <- object$posterior$log_evidence %||% NA_real_
  if (length(log_evidence) == 1L && is.finite(log_evidence)) {
    return(structure(
      as.numeric(log_evidence),
      sampler = object$model$sampler,
      estimate = "log marginal likelihood"
    ))
  }
  reason <- paste(
    "Marginal likelihood is not estimated by the current",
    "RJ-MCMC or PGAS backend."
  )
  if (isTRUE(warn)) warning(reason, call. = FALSE)
  structure(NA_real_, reason = reason)
}

#' Extract test-pattern log predictive density
#'
#' `test` must be supplied to `ppt_fit()`. The returned value is the joint
#' posterior log predictive density of that test point pattern over the fitted
#' observation region.
#'
#' @param object A fitted `ppt` object.
#' @return A numeric scalar.
#' @export
ppt_lppd <- function(object) {
  if (!inherits(object, "ppt")) {
    stop("`object` must inherit from class \"ppt\".", call. = FALSE)
  }
  if (!is.finite(object$posterior$lppd)) {
    stop("No lppd is stored. Refit with a non-empty `test` matrix.",
         call. = FALSE)
  }
  structure(
    object$posterior$lppd,
    n_test = nrow(object$data$test),
    joint = TRUE
  )
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
#' @param x A fitted `ppt` object.
#' @param type Posterior median or mean.
#' @param dims One dimension for a line plot or two dimensions for a colored
#'   point plot.
#' @param main,xlab,ylab Optional labels.
#' @param ... Additional graphical arguments passed to [graphics::plot()].
#' @return The fitted object, invisibly.
#' @method plot ppt
#' @export
plot.ppt <- function(x, type = c("median", "mean"), dims = NULL,
                     main = NULL, xlab = NULL, ylab = NULL, ...) {
  type <- match.arg(type)
  locations <- x$prediction$locations
  d <- ncol(locations)
  if (is.null(dims)) dims <- if (d == 1L) 1L else c(1L, 2L)
  dims <- as.integer(dims)
  if (!length(dims) || length(dims) > 2L ||
      anyNA(dims) || any(dims < 1L | dims > d)) {
    stop("`dims` must select one or two valid input dimensions.",
         call. = FALSE)
  }
  z <- x$prediction[[type]]
  if (is.null(main)) main <- paste(x$model$label, type, "intensity")
  input_names <- colnames(locations)
  if (is.null(input_names)) input_names <- paste0("x", seq_len(d))

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
