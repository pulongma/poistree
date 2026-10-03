#' Fit a Bayesian Poisson point-process tree
#'
#' `ppt_fit()` is the extensible interface for the `poistree` model family.
#' Model components are selected independently: `gating` controls the split
#' mechanism, `scales` controls the intensity representation, and `sampler`
#' controls posterior computation. The current release implements terminal-
#' leaf PPT and multiscale PPT intensity models with hard or soft gating.
#' Hard terminal-leaf PPT is available with SMC, RJ-MCMC, or Particle Gibbs
#' based on an exact conditional-SMC update; the other combinations use
#' RJ-MCMC or informed RJ-MCMC (iRJ-MCMC).
#' The legacy sampler token `pgas` selects Particle Gibbs for compatibility.
#' Ancestor sampling is currently disabled because an exact tree-suffix
#' backward weight has not yet been implemented.
#'
#' All samplers return the same S3 class and output schema. The default call
#' fits S-MPPT with independent multiscale rates and RJ-MCMC. Supplying
#' `scales = "leaf"` without changing `gating` fits S-PPT. Explicit hard-gated
#' calls retain PPT and MPPT.
#'
#' @param x Numeric `n` by `d` matrix containing the observed point-process
#'   inputs. Spatial coordinates, time, and environmental covariates are all
#'   treated as input dimensions.
#' @param region Numeric `d` by 2 matrix containing lower and upper bounds.
#' @param gating Split mechanism: `soft` or `hard`. Soft gating is the
#'   default; supply `gating = "hard"` explicitly for PPT or MPPT.
#' @param scales Intensity structure: terminal `leaf` rates or
#'   `multiscale` rates.
#' @param sampler Posterior sampler: `rjmcmc`, `irjmcmc`, `smc`, or `pgas`.
#'   The `irjmcmc` option is currently available for independent-scale MPPT
#'   and S-MPPT, with either logistic or compact soft gates. The legacy token
#'   `pgas` currently invokes Particle Gibbs
#'   without ancestor sampling.
#' @param scale_prior Multiscale-prior selector: `independent` or `markov`.
#'   It is ignored by terminal-leaf models.
#' @param ... Backend arguments. Common arguments include `predict_at`, `test`,
#'   `max_depth`, and `min_leaf_n`; the minimum leaf occupancy defaults to 1
#'   for every backend. SMC uses `particles`, `a`, `b`, and
#'   `resample_thresh`; hard-leaf SMC and Particle Gibbs also accept
#'   `max_aspect_ratio`, whose default `Inf` imposes no shape restriction on
#'   otherwise valid child regions. RJ-MCMC uses `chains`, `iter`, `burn`,
#'   `cut_candidates`, and `prediction_draws`. Informed RJ-MCMC always combines
#'   a collapsed-score informed cut proposal (with the fitted soft-gate
#'   exposure when applicable) and sequential conditional allocation
#'   proposals. Normalized forward and reverse proposal probabilities are
#'   included in every acceptance ratio.
#'   The implementation informs the cut conditional on the ordinary selection
#'   of move type, eligible node, and split dimension; it does not use the
#'   full-neighborhood informed-importance-tempering transition.
#'   Particle Gibbs uses `particles`, `chains`,
#'   `iter`, and `burn`. Soft models additionally accept `gate_family`,
#'   `gate_structure`, `gate`, and gate-prior controls; multiscale models
#'   accept scale-prior and hyperparameter controls.
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
#'   gating = "hard", scales = "multiscale", sampler = "rjmcmc",
#'   predict_at = x,
#'   chains = 1, iter = 1000, burn = 300
#' )
#' ppt_summary(fit)
#'
#' # Complete S-MPPT surface example with post-hoc intensity evaluation:
#' demo("soft-tree-surface", package = "poistree")
#' }
ppt_fit <- function(x, region,
                    gating = c("soft", "hard"),
                    scales = c("multiscale", "leaf"),
                    sampler = c("rjmcmc", "irjmcmc", "smc", "pgas"),
                    scale_prior = c("independent", "markov"),
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

  if (identical(sampler, "irjmcmc") &&
      (!identical(scales, "multiscale") ||
       !identical(scale_prior, "independent"))) {
    stop(
      "`sampler = \"irjmcmc\"` is currently available only for ",
      "independent-scale MPPT and S-MPPT (`scales = \"multiscale\"`, ",
      "`scale_prior = \"independent\"`).",
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
      "through `ppt_fit()`. Supported configurations are hard + leaf with ",
      "SMC, RJ-MCMC, or Particle Gibbs (legacy token `pgas`); soft + ",
      "leaf + RJ-MCMC; hard or soft + ",
      "multiscale + independent or Markov + RJ-MCMC; and hard or soft + ",
      "multiscale + independent + informed RJ-MCMC.",
      call. = FALSE
    )
  }

  backend <- get(backend_name, mode = "function", inherits = TRUE)
  fit <- backend(x = x, region = region, ...)
  fit$call <- match.call()
  fit
}
