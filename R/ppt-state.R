#' @keywords internal
.ppt_has_state <- function(object) {
  identical(object$posterior$state$mode, "heap") ||
    length(object$posterior$tree_draws) > 0L
}

#' @keywords internal
.ppt_state_eval <- function(object, at) {
  at <- .ppt_validate_points(
    at, object$data$dimension, object$data$region, "at",
    allow_empty = FALSE
  )
  state <- object$posterior$state
  if (identical(state$mode, "heap")) {
    return(ppt_eval_state(
      state$nodes, state$gate, object$data$region, at,
      as.integer(state$gate_mode), as.numeric(state$gate_depth %||% 0)
    ))
  }
  trees <- object$posterior$tree_draws
  if (length(trees)) {
    return(.ppt_eval_leafbox(trees, at, object$data$region))
  }
  stop(
    "This fit does not carry posterior state draws; refit with the ",
    "current poistree version to enable evaluation at new locations.",
    call. = FALSE
  )
}

# Reuse retained prediction draws; evaluate only locations missing from them.
# .ppt_state_eval remains the independent, unconditional state evaluator.
#' @keywords internal
.ppt_prediction_draws <- function(object, at) {
  stored <- object$prediction$draws
  locations <- object$prediction$locations
  if (!is.matrix(stored) || is.null(locations) ||
      nrow(stored) != nrow(locations) || !ncol(stored)) {
    return(.ppt_state_eval(object, at))
  }
  state <- object$posterior$state
  state_draws <- if (identical(state$mode, "heap")) length(state$nodes) else {
    length(object$posterior$tree_draws)
  }
  if (state_draws > 0L && state_draws != ncol(stored)) {
    return(.ppt_state_eval(object, at))
  }
  # Full-grid requests need neither row keys nor a copy of the draw matrix.
  if (identical(dim(at), dim(locations)) && all(at == locations)) return(stored)
  index <- .ppt_match_prediction_rows(at, locations)
  known <- !is.na(index)
  if (!any(known)) return(.ppt_state_eval(object, at))
  if (all(known)) {
    if (identical(index, seq_len(nrow(stored)))) return(stored)
    return(stored[index, , drop = FALSE])
  }
  draws <- matrix(0, nrow(at), ncol(stored))
  draws[known, ] <- stored[index[known], , drop = FALSE]
  draws[!known, ] <- .ppt_state_eval(object, at[!known, , drop = FALSE])
  draws
}

# Integrate every retained heap-state intensity draw over all coordinates
# except `variable`.  The native evaluator reconstructs node boxes and soft
# ancestor paths from the serialized topology, then uses exact box widths,
# piecewise-polynomial compact-gate integrals, or the stable logistic path
# integral used by the fitting backends (adaptive for node-scaled logistic).
#' @keywords internal
.ppt_state_marginal <- function(object, grid, variable, average = TRUE) {
  state <- object$posterior$state
  if (!identical(state$mode, "heap") || !length(state$nodes)) {
    stop(
      "This fit does not carry heap-state posterior draws for exact ",
      "marginal integration.",
      call. = FALSE
    )
  }
  ppt_marginal_state(
    state$nodes,
    as.matrix(state$gate),
    object$data$region,
    as.numeric(grid),
    as.integer(variable - 1L),
    as.integer(state$gate_mode),
    as.numeric(state$gate_depth %||% 0),
    isTRUE(average)
  )
}

# Evaluate hard terminal-leaf tree draws at arbitrary locations by leaf-box
# membership.  The half-open convention (lower closed, upper open, except at
# the upper region boundary) reproduces the sampler's "left iff x < cut"
# routing exactly, so boundary points land in the same leaf as prediction
# during fitting.
#' @keywords internal
.ppt_eval_leafbox <- function(trees, at, region) {
  d <- nrow(region)
  upper <- region[, 2L]
  tolerance <- sqrt(.Machine$double.eps) * pmax(1, abs(upper))
  n <- nrow(at)
  draws <- matrix(0, n, length(trees))
  for (s in seq_along(trees)) {
    leaves <- Filter(
      function(node) !is.null(node) && isTRUE(node$is_leaf), trees[[s]]
    )
    if (!length(leaves)) {
      stop("A posterior tree draw contains no terminal nodes.",
           call. = FALSE)
    }
    value <- numeric(n)
    for (leaf in leaves) {
      box <- as.matrix(leaf$region)
      lambda <- as.numeric(leaf$lambda)
      if (!identical(dim(box), c(d, 2L)) || length(lambda) != 1L ||
          !is.finite(lambda) || lambda < 0) {
        stop("A posterior tree draw contains an invalid leaf.",
             call. = FALSE)
      }
      inside <- rep(TRUE, n)
      for (j in seq_len(d)) {
        is_last <- abs(box[j, 2L] - upper[j]) <= tolerance[j]
        inside <- inside & at[, j] >= box[j, 1L] &
          (at[, j] < box[j, 2L] |
             (is_last & at[, j] <= box[j, 2L] + tolerance[j]))
      }
      value[inside] <- lambda
    }
    draws[, s] <- value
  }
  draws
}

#' Evaluate the posterior intensity at arbitrary locations
#'
#' Every fit stores compact posterior state draws (tree topology, node rates,
#' and gate vectors), so the intensity \eqn{\lambda(s)} of each retained draw
#' can be reproduced exactly at any location after fitting -- no location has
#' to be listed in `predict_at` beforehand. This is the workhorse behind
#' plotting the estimated intensity over a grid and behind post-hoc
#' [ppt_predict()] and [ppt_lppd()] evaluation. Locations already present in
#' the stored prediction draws reuse those draws. Other locations alone are
#' evaluated from the posterior states; input order and duplicate rows are kept.
#'
#' @param object A fitted `ppt` object.
#' @param at Numeric matrix of evaluation locations with one column per input.
#'   When `NULL`, an equally spaced grid over the fitted region is built:
#'   `n` points in one dimension, an `n` by `n` lattice in two dimensions.
#'   Supply `at` explicitly for more than two dimensions.
#' @param n Grid resolution per dimension used when `at` is `NULL`.
#' @param type Return pointwise posterior summaries or the raw
#'   locations-by-draws matrix.
#' @param level Pointwise posterior credible level for `type = "summary"`.
#'
#' @return For `type = "summary"`, a data frame with the input coordinates
#'   and columns `mean`, `median`, `lower`, and `upper` -- ready for
#'   `ggplot2` plotting of the estimated intensity surface. For
#'   `type = "draws"`, a list with elements `at`, `draws` (locations by
#'   draws), and `weights` (posterior draw weights; uniform for MCMC,
#'   particle weights for SMC).
#'
#' @examples
#' \dontrun{
#' fit <- ppt_fit(x, region, chains = 2, iter = 2000, burn = 500)
#' surface <- ppt_lambda(fit, n = 80)
#' ggplot2::ggplot(surface, ggplot2::aes(x1, x2, fill = mean)) +
#'   ggplot2::geom_raster() +
#'   ggplot2::scale_fill_viridis_c()
#' }
#' @export
ppt_lambda <- function(object, at = NULL, n = 50L,
                       type = c("summary", "draws"), level = 0.95) {
  if (!inherits(object, "ppt")) {
    stop("`object` must inherit from class \"ppt\".", call. = FALSE)
  }
  type <- match.arg(type)
  if (length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must lie strictly between 0 and 1.", call. = FALSE)
  }
  d <- object$data$dimension
  region <- object$data$region
  input_names <- colnames(object$data$x)
  if (is.null(input_names)) input_names <- paste0("x", seq_len(d))

  if (is.null(at)) {
    if (length(n) != 1L || !is.finite(n) || n != floor(n) || n < 2L) {
      stop("`n` must be an integer of at least 2.", call. = FALSE)
    }
    n <- as.integer(n)
    if (d == 1L) {
      at <- matrix(seq(region[1L, 1L], region[1L, 2L], length.out = n),
                   ncol = 1L)
    } else if (d == 2L) {
      at <- as.matrix(expand.grid(
        seq(region[1L, 1L], region[1L, 2L], length.out = n),
        seq(region[2L, 1L], region[2L, 2L], length.out = n)
      ))
    } else {
      stop("Supply `at` explicitly for more than two input dimensions.",
           call. = FALSE)
    }
    colnames(at) <- input_names
  }
  at <- .ppt_validate_points(at, d, region, "at", allow_empty = FALSE)

  draws <- .ppt_prediction_draws(object, at)
  weights <- .ppt_posterior_weights(object, ncol(draws))
  if (identical(type, "draws")) {
    return(list(at = at, draws = draws, weights = weights))
  }
  alpha <- (1 - level) / 2
  quantiles <- t(apply(draws, 1L, .ppt_weighted_quantiles,
                      weights = weights, probabilities = c(0.5, alpha, 1 - alpha)))
  out <- data.frame(at, check.names = FALSE)
  names(out) <- input_names
  out$mean <- as.numeric(draws %*% weights)
  out$median <- quantiles[, 1L]
  out$lower <- quantiles[, 2L]
  out$upper <- quantiles[, 3L]
  attr(out, "level") <- level
  class(out) <- c("ppt_lambda", "data.frame")
  out
}
