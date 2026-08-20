#' Fit a Bayesian Poisson point-process tree
#'
#' `ppt_fit()` is the extensible interface for the `poistree` model family.
#' Model components are selected independently: `gating` controls the split
#' mechanism, `scales` controls the intensity representation, and `sampler`
#' controls posterior computation. The current release implements the
#' hard-gated terminal-leaf PPT with SMC, RJ-MCMC, or Particle Gibbs with
#' ancestor sampling (PGAS); the component
#' arguments are retained so additional backends can be registered without
#' changing the public interface.
#'
#' All samplers return the same S3 class and output schema.
#'
#' @param x Numeric `n` by `d` matrix containing the observed point-process
#'   inputs. Spatial coordinates, time, and environmental covariates are all
#'   treated as input dimensions.
#' @param region Numeric `d` by 2 matrix containing lower and upper bounds.
#' @param gating Split mechanism. Only `hard` is currently implemented;
#'   `soft` is reserved for a future backend.
#' @param scales Intensity structure: terminal `leaf` rates or
#'   `multiscale` rates.
#' @param sampler Posterior sampler: `rjmcmc`, `smc`, or `pgas`.
#' @param scale_prior Reserved multiscale-prior selector: `markov` or
#'   `independent`. It is ignored by the currently implemented leaf models.
#' @param ... Backend arguments. Common arguments include `predict_at`, `test`,
#'   `max_depth`, and `min_leaf_n`. SMC uses `particles`, `a`, `b`, and
#'   `resample_thresh`. RJ-MCMC uses `chains`, `iter`, `burn`, `cut_candidates`,
#'   and `prediction_draws`. PGAS uses `particles`, `chains`, `iter`, and `burn`.
#'
#' @return An object of S3 class `ppt`.
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' x <- matrix(runif(200), ncol = 2)
#' region <- matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE)
#' fit <- ppt_fit(
#'   x, region,
#'   gating = "hard", scales = "leaf", sampler = "rjmcmc",
#'   predict_at = x,
#'   chains = 1, iter = 1000, burn = 300
#' )
#' ppt_summary(fit)
#' }
ppt_fit <- function(x, region,
                    gating = c("hard", "soft"),
                    scales = c("leaf", "multiscale"),
                    sampler = c("rjmcmc", "smc", "pgas"),
                    scale_prior = c("markov", "independent"),
                    ...) {
  gating <- match.arg(gating)
  scales <- match.arg(scales)
  sampler <- match.arg(sampler)
  scale_prior <- match.arg(scale_prior)

  if (identical(sampler, "pgas") &&
      (!identical(gating, "hard") || !identical(scales, "leaf"))) {
    stop(
      "`sampler = \"pgas\"` is available only for the hard-gated, ",
      "terminal-leaf PPT model (`gating = \"hard\"`, `scales = \"leaf\"`).",
      call. = FALSE
    )
  }

  backend_key <- .ppt_backend_key(gating, scales, sampler, scale_prior)
  backend_name <- unname(.ppt_backend_registry[backend_key])
  if (!length(backend_name) || is.na(backend_name)) {
    requested <- paste(
      c(gating, scales,
        if (identical(scales, "multiscale")) scale_prior,
        sampler),
      collapse = " + "
    )
    stop(
      "The requested configuration (", requested, ") is not yet available ",
      "through `ppt_fit()`. Currently supported: hard + leaf + smc; ",
      "hard + leaf + rjmcmc; and hard + leaf + pgas. Other component ",
      "combinations are reserved for future model backends.",
      call. = FALSE
    )
  }

  backend <- get(backend_name, mode = "function", inherits = TRUE)
  fit <- backend(x = x, region = region, ...)
  fit$call <- match.call()
  fit
}
