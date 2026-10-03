#' Fit a Bayesian Poisson point-process tree
#'
#' `ppt_fit()` fits terminal-leaf intensity models with hard or soft gating.
#' Hard PPT supports SMC, RJ-MCMC, or Particle Gibbs based on an exact
#' conditional-SMC update. Soft S-PPT uses RJ-MCMC and is the default.
#' All samplers return the same S3 class and output schema.
#' The sampler token `pgas` selects Particle Gibbs for compatibility;
#' ancestor sampling is currently disabled because an exact tree-suffix
#' backward weight has not yet been implemented.
#'
#' @param x Numeric `n` by `d` matrix containing the observed point-process
#'   inputs. Spatial coordinates, time, and environmental covariates are all
#'   treated as input dimensions.
#' @param region Numeric `d` by 2 matrix containing lower and upper bounds.
#' @param gating Split mechanism: `soft` (the default) or `hard`.
#' @param scales Intensity structure: terminal `leaf` rates.
#' @param sampler Posterior sampler: `rjmcmc`, `smc`, or `pgas`.
#'   SMC and Particle Gibbs are available for hard gates only.
#' @param ... Backend arguments. Common arguments include `predict_at`, `test`,
#'   `max_depth`, and `min_leaf_n`; the minimum leaf occupancy defaults to 1.
#'   SMC uses `particles`, `a`, `b`, and `resample_thresh`. Hard-leaf SMC and
#'   Particle Gibbs also accept `max_aspect_ratio`, whose default `Inf`
#'   imposes no shape restriction on otherwise valid child regions.
#'   RJ-MCMC uses `chains`, `iter`, `burn`, and `cut_candidates`; hard
#'   RJ-MCMC also accepts `prediction_draws`. Particle Gibbs uses `particles`,
#'   `chains`, `iter`, and `burn`. Soft models additionally accept
#'   `gate_family`, `gate_structure`, `gate`, and gate-prior controls.
#'
#' @return An object of S3 class `ppt`.
#' @export
#' @examples
#' \dontrun{
#' set.seed(1)
#' x <- matrix(runif(200), ncol = 2)
#' region <- matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE)
#' fit <- ppt_fit(x, region, gating = "hard", sampler = "rjmcmc",
#'                predict_at = x, chains = 1, iter = 1000, burn = 300)
#' ppt_summary(fit)
#' demo("soft-tree-surface", package = "poistree")
#' }
ppt_fit <- function(x, region,
                    gating = c("soft", "hard"),
                    scales = "leaf",
                    sampler = c("rjmcmc", "smc", "pgas"),
                    ...) {
  gating <- match.arg(gating)
  scales <- match.arg(scales, "leaf")
  sampler <- match.arg(sampler)
  if (identical(sampler, "pgas") && !identical(gating, "hard")) {
    stop(
      "`sampler = \"pgas\"` is available only for the hard-gated, ",
      "terminal-leaf PPT model (`gating = \"hard\"`, `scales = \"leaf\"`).",
      call. = FALSE
    )
  }
  backend_key <- .ppt_backend_key(gating, scales, sampler)
  backend_name <- unname(.ppt_backend_registry[backend_key])
  if (!length(backend_name) || is.na(backend_name)) {
    stop(
      "The requested configuration (", paste(gating, scales, sampler, sep = " + "),
      ") is not yet available through `ppt_fit()`. Supported configurations ",
      "are hard + leaf with SMC, RJ-MCMC, or Particle Gibbs (token `pgas`); ",
      "and soft + leaf + RJ-MCMC.", call. = FALSE
    )
  }
  backend <- get(backend_name, mode = "function", inherits = TRUE)
  fit <- backend(x = x, region = region, ...)
  fit$call <- match.call()
  fit
}
