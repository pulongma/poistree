#' Summarize sampler diagnostics
#'
#' Reports the acceptance information exposed by the selected backend together
#' with the number of retained draws. Full scalar traces will be added when the
#' backends retain them in fitted objects.
#'
#' @param object A fitted `ppt` object.
#' @param plot Draw an acceptance-rate bar plot.
#' @return An object of class `ppt_diagnostics`.
#' @export
ppt_diagnostics <- function(object, plot = FALSE) {
  if (!inherits(object, "ppt")) {
    stop("`object` must inherit from class \"ppt\".", call. = FALSE)
  }
  out <- structure(
    list(
      model = object$model$label,
      sampler = object$model$sampler,
      scale_prior = object$model$scale_prior,
      draws = object$posterior$draws,
      mean_leaves = object$posterior$mean_leaves,
      mean_max_depth = object$posterior$mean_max_depth %||% NA_real_,
      kappa = object$posterior$kappa %||% NA_real_,
      tau = object$posterior$tau %||% NA_real_,
      acceptance = object$diagnostics$acceptance,
      gate_by_dimension = object$posterior$gate_by_dimension,
      particle_ess = object$diagnostics$particle_ess %||% NA_real_,
      ess_history = object$diagnostics$ess_history %||% numeric(),
      unique_trees = object$diagnostics$unique_trees %||% NA_integer_
    ),
    class = "ppt_diagnostics"
  )
  if (isTRUE(plot)) plot(out)
  out
}

#' @method print ppt_diagnostics
#' @export
print.ppt_diagnostics <- function(x, ...) {
  cat("poistree sampler diagnostics\n")
  cat("  Model       :", x$model, "\n")
  cat("  Sampler     :", x$sampler, "\n")
  cat("  Draws       :", x$draws, "\n")
  cat("  Mean leaves :", format(x$mean_leaves, digits = 5L), "\n")
  if (is.finite(x$mean_max_depth)) {
    cat("  Max depth   :", format(x$mean_max_depth, digits = 5L), "\n")
  }
  if (is.finite(x$kappa)) {
    cat("  Kappa       :", format(x$kappa, digits = 5L), "\n")
  }
  if (is.finite(x$tau)) {
    cat("  Tau         :", format(x$tau, digits = 5L), "\n")
  }
  if (is.finite(x$particle_ess)) {
    cat("  Particle ESS:", format(x$particle_ess, digits = 5L), "\n")
  }
  if (is.finite(x$unique_trees)) {
    cat("  Unique trees:", x$unique_trees, "\n")
  }
  if (length(x$acceptance)) {
    cat("  Acceptance:\n")
    print(round(x$acceptance, 3L))
  } else {
    cat("  Acceptance  : not returned by this summary backend\n")
  }
  invisible(x)
}

#' @method plot ppt_diagnostics
#' @export
plot.ppt_diagnostics <- function(x, ...) {
  if (x$sampler %in% c("smc", "pgas")) {
    ess <- as.numeric(x$ess_history)
    ess <- ess[is.finite(ess) & ess > 0]
    if (!length(ess)) {
      stop("This SMC fit did not return an ESS trajectory.",
           call. = FALSE)
    }
    graphics::plot(
      seq_along(ess), ess, type = "b",
      xlab = "SMC step", ylab = "Effective sample size",
      main = paste(x$model, toupper(x$sampler), "diagnostics"), ...
    )
    return(invisible(x))
  }
  if (!length(x$acceptance)) {
    stop("This backend did not return acceptance-rate summaries.",
         call. = FALSE)
  }
  graphics::barplot(
    x$acceptance,
    ylim = c(0, 1),
    ylab = "Acceptance rate",
    main = paste(x$model, "RJ-MCMC acceptance"),
    las = 2,
    ...
  )
  invisible(x)
}

#' Simulate a Poisson point process by thinning
#'
#' @param intensity Function accepting an `n` by `d` matrix and returning
#'   nonnegative intensities.
#' @param region Numeric `d` by 2 observation bounds.
#' @param lambda_max Finite upper bound for `intensity` over `region`.
#' @param nsim Number of independent point patterns.
#' @param seed Optional random seed.
#' @return A matrix when `nsim = 1`, otherwise a list of matrices.
#' @export
ppt_sim <- function(intensity, region, lambda_max, nsim = 1L, seed = NULL) {
  if (!is.function(intensity)) {
    stop("`intensity` must be a function of an n by d matrix.",
         call. = FALSE)
  }
  region <- as.matrix(region)
  storage.mode(region) <- "double"
  if (ncol(region) != 2L || !nrow(region) ||
      any(!is.finite(region)) || any(region[, 2L] <= region[, 1L])) {
    stop("`region` must be a finite d by 2 matrix.", call. = FALSE)
  }
  if (length(lambda_max) != 1L || !is.finite(lambda_max) ||
      lambda_max <= 0) {
    stop("`lambda_max` must be a positive finite scalar.", call. = FALSE)
  }
  nsim <- as.integer(nsim)
  if (length(nsim) != 1L || is.na(nsim) || nsim < 1L) {
    stop("`nsim` must be a positive integer.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  d <- nrow(region)
  volume <- prod(region[, 2L] - region[, 1L])
  one_simulation <- function() {
    n_proposal <- stats::rpois(1L, lambda_max * volume)
    if (!n_proposal) return(matrix(numeric(0), 0L, d))
    proposal <- vapply(seq_len(d), function(j) {
      stats::runif(n_proposal, region[j, 1L], region[j, 2L])
    }, numeric(n_proposal))
    proposal <- matrix(proposal, ncol = d)
    lambda <- as.numeric(intensity(proposal))
    if (length(lambda) != n_proposal || any(!is.finite(lambda)) ||
        any(lambda < 0)) {
      stop("`intensity` must return one finite nonnegative value per row.",
           call. = FALSE)
    }
    tolerance <- sqrt(.Machine$double.eps) * max(1, lambda_max)
    if (any(lambda > lambda_max + tolerance)) {
      stop("`lambda_max` is smaller than an evaluated intensity.",
           call. = FALSE)
    }
    keep <- stats::runif(n_proposal) <
      pmin(1, lambda / lambda_max)
    proposal[keep, , drop = FALSE]
  }

  result <- lapply(seq_len(nsim), function(i) one_simulation())
  if (nsim == 1L) result[[1L]] else result
}
