#' Fit a Bayesian Poisson point-process tree
#'
#' `ppt_fit()` is the extensible interface for the `poistree` model family.
#' Model components are selected independently: `gating` controls the split
#' mechanism and `sampler` controls posterior computation. The current
#' release implements two terminal-leaf intensity models: the hard-gated
#' PPT, available with SMC, RJ-MCMC, informed MH, or Particle Gibbs based on an exact
#' conditional-SMC update, and the soft-gated S-PPT, available with RJ-MCMC,
#' informed MH,
#' or Particle Gibbs with ancestor sampling (PGAS), whose conditional SMC
#' runs over the tree and the latent allocation variables of the observations jointly on an
#' extended complete-tree space. For the hard PPT the `pgas` token runs
#' conditional SMC without ancestor sampling.
#'
#' All samplers return the same S3 class and output schema. The default call
#' fits S-PPT with RJ-MCMC; `gating = "hard"` selects PPT.
#'
#' @param x Numeric `n` by `d` matrix containing the observed point-process
#'   inputs. Spatial coordinates, time, and environmental covariates are all
#'   treated as input dimensions.
#' @param region Numeric `d` by 2 matrix containing lower and upper bounds.
#' @param gating Split mechanism: `soft` (the default, S-PPT) or `hard`
#'   (PPT).
#' @param scales Intensity structure. Only terminal `leaf` rates are
#'   implemented; the argument is reserved so that existing calls with
#'   `scales = "leaf"` keep working.
#' @param sampler Posterior sampler: `rjmcmc`, `irjmcmc` (cached locally
#'   informed Metropolis--Hastings), `smc` (hard gating only), or
#'   `pgas`. With `gating = "soft"` the `pgas` token runs Particle Gibbs
#'   with ancestor sampling; with `gating = "hard"` it runs conditional SMC
#'   without ancestor sampling.
#' @param ... Backend arguments. All backends accept `predict_at`, `test`,
#'   `max_depth`, and `cut_candidates`. The positive integer `min_leaf_n`
#'   defaults to 1 in all hard backends and in soft RJ-MCMC and informed MH;
#'   it screens candidate splits by their hard-routed child counts. Soft
#'   PGAS does not accept `min_leaf_n`.
#'   `cut_candidates` defaults to 50 for every backend;
#'   for a quantile proposal this is the requested number of probabilities
#'   before duplicate and inadmissible cuts are removed. It must be an integer
#'   of at least 1 for hard SMC and hard Particle Gibbs, and at least 2 for
#'   the other backends. Soft RJ-MCMC and
#'   informed MH also accept `cut_proposal = "uniform"`, whose grid retains
#'   the existing minimum of 30 locations before filtering, and
#'   `cut_proposal = "data"`, which uses all admissible data midpoints and
#'   ignores `cut_candidates`. The default quantile and uniform grids thus
#'   request 50 locations. SMC uses `particles`,
#'   `a`, `b`, `resample_thresh`, `cut_candidates` (per input and node), and
#'   `engine` (`"shared"`, the default, runs the sweep on a shared-path store
#'   of boxes so that particles reaching the same box share its candidate
#'   cuts and scores; `"dense"` is the original per-particle implementation,
#'   identical in law); hard-leaf SMC and Particle Gibbs also accept
#'   `max_aspect_ratio`, whose default `Inf` imposes no shape restriction on
#'   otherwise valid child regions. RJ-MCMC and informed MH use `chains`,
#'   `iter`, `burn`, and `cut_candidates`. Their hard backends also accept
#'   `prediction_draws`, a positive integer controlling retained predictions
#'   and tree states. Their soft backends use `thin`, `tree_moves`, and
#'   `change_moves`, and do not accept `prediction_draws`.
#'   Particle Gibbs uses `particles`, `chains`,
#'   `iter`, `burn`, and `cut_candidates`. The soft PGAS backend also accepts `thin`,
#'   `label_sweeps` (label Gibbs sweeps per iteration), `ancestor_sampling`,
#'   `exact_max` (largest node occupancy for which the exact Poisson-binomial
#'   one-step-ahead proposal is used, default 150; above it the same
#'   expectation is evaluated by a Laplace approximation of its Beta-integral
#'   representation, O(m) per candidate),
#'   `defensive` (mixture weight on the prior action proposal),
#'   `resampling` (`"node"`, the default, resamples after every heap node at
#'   which some particle advanced; `"level"` resamples after every tree
#'   level), `ess_threshold` (resample at such an event only when the ESS is
#'   at most this fraction of `particles`; default 1, i.e. always),
#'   `allocation` (`"sequential"`, the default, allocates the points of a node
#'   with more than `exact_max` points to its children by sequential
#'   imputation; `"rates"` draws auxiliary child rates and allocates the
#'   points independently given them), and
#'   `cut_grid` (a list of fixed cut locations per input; by default
#'   `cut_candidates` global quantiles). All soft backends accept
#'   `gate_structure` (`"dimension"`, the default, one gating
#'   parameter per input updated by a systematic Metropolis scan over the
#'   inputs in every iteration; or `"shared"`, one common parameter), `gate`,
#'   and gate-prior controls, which are scalars or vectors of length `d`.
#'   Only soft RJ-MCMC and informed MH accept `gate_family`, with choices
#'   `"logistic"` (the default) and `"compact"`. Soft PGAS uses logistic gates
#'   and does not accept a `gate_family` argument.
#'
#' @details
#' `max_depth` must be one finite integer. All hard backends and soft RJ-MCMC
#' and informed MH accept values from 0 through 20; zero keeps only the root.
#' Soft PGAS requires at least 1 and its current reference-tree storage limit
#' permits at most 19. Dense hard SMC and hard Particle Gibbs also impose a
#' tree-storage limit that depends on both depth and particle count. If this
#' limit is exceeded, reduce `max_depth` or, for those hard backends,
#' `particles`. These guards cover tree indexing and preallocated tree
#' storage; they do not bound all cache and output memory. Defaults are
#' unchanged by these limits.
#'
#' With `sampler = "irjmcmc"`, the current RJ-MCMC target and admissible tree
#' support are retained. Each reversible tree kernel uses square-root
#' balanced neighbor weights. For target density \eqn{\pi}, base proposal
#' \eqn{q_0}, and neighboring states \eqn{s,t}, define
#' \deqn{r(s,t)=\frac{\pi(t)q_0(t,s)}{\pi(s)q_0(s,t)},\qquad
#'   w(s,t)=q_0(s,t)\sqrt{r(s,t)},\qquad Z(s)=\sum_t w(s,t).}
#' A candidate is proposed with probability \eqn{w(s,t)/Z(s)} and accepted
#' with probability \eqn{\min\{1,Z(s)/Z(t)\}}. A rejection retains the current
#' state. Draws therefore have ordinary posterior weights; this sampler does
#' not use importance-tempering weights.
#'
#' The hard sampler scores its grow, prune, and terminal-split change
#' neighborhood using exact Gamma--Poisson collapsed likelihoods. Candidate
#' cuts and local scores are cached, and the current neighborhood is reused
#' after rejection. The soft sampler updates the joint tree-and-allocation
#' state at fixed gates. It sums allocation contributions exactly using
#' Poisson-binomial count recursion, followed by an exact conditional draw
#' of the allocation. For \eqn{m} allocated observations at a candidate node,
#' this count recursion costs \eqn{O(m^2)} per candidate.
#' `tree_moves` and `change_moves` retain their meanings
#' as counts of separate grow/prune and change kernels per iteration.
#' Allocation Gibbs sweeps and gate Metropolis updates remain part of the
#' soft sampler. Exact scoring can be costly for large nodes; informed MH
#' is not guaranteed to improve effective samples per second.
#'
#' The default candidate count is shared by all samplers, but their existing
#' cut-generation rules and tree priors remain different. In particular,
#' soft PGAS uses a global cut grid and RJ-MCMC uses node-specific candidates.
#' The same candidate count does not make different sampler families target
#' the same tree model. Explicit `cut_candidates` values follow the
#' proposal-specific rules above; a supplied soft PGAS `cut_grid` overrides
#' its automatically generated quantile grid.
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
#'   gating = "hard", sampler = "rjmcmc",
#'   predict_at = x,
#'   chains = 1, iter = 1000, burn = 300
#' )
#' ppt_summary(fit)
#'
#' # Complete S-PPT surface example with post-hoc intensity evaluation:
#' demo("soft-tree-surface", package = "poistree")
#' }
ppt_fit <- function(x, region,
                    gating = c("soft", "hard"),
                    scales = "leaf",
                    sampler = c("rjmcmc", "irjmcmc", "smc", "pgas"),
                    ...) {
  gating <- match.arg(gating)
  scales <- match.arg(scales)
  sampler <- match.arg(sampler)

  backend_key <- .ppt_backend_key(gating, scales, sampler)
  backend_name <- unname(.ppt_backend_registry[backend_key])
  if (!length(backend_name) || is.na(backend_name)) {
    stop(
      "The requested configuration (",
      paste(c(gating, scales, sampler), collapse = " + "),
      ") is not available through `ppt_fit()`. Supported configurations ",
      "are hard + leaf with SMC, RJ-MCMC, informed MH (`irjmcmc`), or ",
      "Particle Gibbs (`pgas`), and soft + leaf with RJ-MCMC, informed MH, ",
      "or Particle Gibbs with ancestor sampling ",
      "(`pgas`).",
      call. = FALSE
    )
  }

  backend <- get(backend_name, mode = "function", inherits = TRUE)
  fit <- backend(x = x, region = region, ...)
  fit$call <- match.call()
  fit
}
