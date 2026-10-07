#' Fit a Bayesian Poisson point-process tree
#'
#' `ppt_fit()` is the extensible interface for the `poistree` model family.
#' Model components are selected independently: `gating` controls the split
#' mechanism and `sampler` controls posterior computation. The current
#' release implements two terminal-leaf intensity models: the hard-gated
#' PPT, available with SMC, RJ-MCMC, informed MH, or Particle Gibbs based on an exact
#' conditional-SMC update, and the soft-gated S-PPT, available with RJ-MCMC,
#' informed MH, partially collapsed Gibbs with joint adaptive gate updates,
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
#'   informed Metropolis--Hastings), `pcg` (partially collapsed Gibbs with
#'   robust adaptive Metropolis gate updates, soft gating only),
#'   `smc` (hard gating only), or
#'   `pgas`. With `gating = "soft"` the `pgas` token runs Particle Gibbs
#'   with ancestor sampling; with `gating = "hard"` it runs conditional SMC
#'   without ancestor sampling.
#' @param ... Backend arguments. All backends accept `predict_at`, `test`,
#'   `max_depth`, and `cut_candidates`. The positive integer `min_leaf_n`
#'   defaults to 1 in all hard backends and in soft RJ-MCMC, informed MH, and PCG;
#'   it screens candidate splits by their hard-routed child counts. Soft
#'   PGAS does not accept `min_leaf_n`.
#'   `cut_candidates` defaults to 50 for every backend;
#'   for a quantile proposal this is the requested number of probabilities
#'   before duplicate and inadmissible cuts are removed. It must be an integer
#'   of at least 1 for hard SMC and hard Particle Gibbs, and at least 2 for
#'   the other backends. Soft RJ-MCMC, informed MH, and PCG also accept `cut_proposal = "uniform"`, whose grid retains
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
#'   `change_moves`, and do not accept `prediction_draws`. PCG accepts the same
#'   controls as soft RJ-MCMC, plus `ram_target` (default 0.234),
#'   `ram_decay` (default 0.7), and `ram_adapt` (default `floor(0.8 * burn)`).
#'   Soft RJ-MCMC, informed MH, and PCG also accept `cache_geometry = TRUE`.
#'   This reuses unchanged leaf exposures and computes point memberships
#'   by sharing ancestor-gate calculations across leaves. It supports both
#'   `gate_scale = "root"` and `"node"`, and the compact gate family.
#'   Set `cache_geometry = FALSE` to use the original uncached calculations
#'   for numerical comparisons. The cache is local to each fit, requires no
#'   saved files, and does not change the model, priors, or sampler moves.
#'   Floating-point rounding can differ because products are shared.
#'   Particle Gibbs uses `particles`, `chains`,
#'   `iter`, `burn`, and `cut_candidates`. The soft PGAS backend also accepts `thin`,
#'   `label_sweeps` (label Gibbs sweeps per iteration), `ancestor_sampling`,
#'   `exact_max` (largest node occupancy for which the exact Poisson-binomial
#'   one-step-ahead proposal is used, default 150). Above that threshold,
#'   `proposal_score = "laplace"` (the default) uses the previous Laplace score;
#'   `proposal_score = "hard"` ranks actions using cumulative hard-split counts
#'   and exposures under the current parent path. This is a proposal surrogate
#'   for the soft model: selected splits still use true logistic gates and the
#'   full importance correction. Hard scores use `proposal_temperature`
#'   (default 0.5, in (0, 1]) and a prior-action mixture `proposal_defensive`
#'   (default 0.1, in (0, 1)); these controls do not affect exact-scored nodes.
#'   The hard option evaluates true routing gates lazily and retains all cuts.
#'   Its speed and effective-sample-size tradeoff depends on the data and gates.
#'   `defensive` controls the prior-action mixture for exact and Laplace scoring,
#'   separately from the hard scorer's `proposal_defensive`. Further controls are
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
#'   parameter per input; or `"shared"`, one common parameter), `gate`,
#'   and gate-prior controls, which are scalars or vectors of length `d`.
#'   Soft RJ-MCMC, informed MH, and PCG accept `gate_family`, with choices
#'   `"logistic"` (the default) and `"compact"`. Soft PGAS uses logistic gates
#'   and does not accept a `gate_family` argument.
#'   All soft backends accept `gate_scale = NULL`, independently of
#'   `gate_structure`. For logistic gates, `NULL` retains the historical
#'   `"root"` scaling: the gate on input `j` is
#'   `plogis(gate[j] * (x[j] - cut) / root_width[j])`.
#'   Use `gate_scale = "node"` to divide by the splitting node's local
#'   width instead, so the same gate parameter gives a sharper transition
#'   in a narrower node. `gate_structure = "dimension"` continues to share
#'   one parameter per input across its nodes; `"shared"` uses one parameter
#'   across all inputs. Compact gates keep their historical node scaling
#'   when `gate_scale` is `NULL` or `"node"`; `"root"` is not supported.
#'   Node-relative logistic path integrals use deterministic adaptive
#'   quadrature. Stored posterior states retain the scale choice, and older
#'   fitted objects continue to use their stored gate mode.
#'   Node scaling for logistic gates is available with `sampler = "pcg"`,
#'   `"rjmcmc"`, and `"irjmcmc"`. Soft PGAS currently supports only root
#'   scaling; it rejects `gate_scale = "node"` before fitting because its
#'   fixed-grid ancestor sampler requires a different derivation for
#'   node-dependent widths.
#'   All MCMC backends accept `verbose` (default `TRUE`), which displays an
#'   RcppProgress bar for each chain. The bar counts completed iterations,
#'   including burn-in and iterations discarded by thinning. Use
#'   `verbose = FALSE` to suppress progress bars and chain messages.
#'
#' @details
#' `max_depth` must be one finite integer. All hard backends and soft RJ-MCMC
#' and informed MH, as well as PCG, accept values from 0 through 20;
#' zero keeps only the root.
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
#' With `sampler = "pcg"`, each iteration has the following order:
#' draw temporary leaf intensities conditional on the current tree, labels,
#' and gates; propose all free log gates jointly with the labels integrated
#' out; restore every observation label by an independent categorical draw
#' conditional on the resulting gates and temporary intensities, including
#' after a rejected gate proposal; discard the temporary intensities; run
#' standard collapsed grow/prune and change tree updates; and draw fresh
#' intensities when retaining output. The tree target and candidate rules are
#' those of soft RJ-MCMC.
#'
#' With `sampler = "pcg"` and `informed = TRUE`, the grow, prune, and change
#' proposals are informed by a hard surrogate instead of uniform draws. Every
#' candidate split of a node is scored once, from the node's hard-routed
#' counts and the box exposures of its two children (the one-step lookahead
#' score of the hard SMC, with the same tree prior), and the proposal mixes a
#' tempered draw \eqn{\propto \exp(	au\,	ilde s)} with a uniform draw over
#' all candidates: `proposal_temperature` is \eqn{	au\in(0,1]} (default 0.5,
#' the square-root balancing function of Zanella 2020) and
#' `proposal_defensive` is the uniform weight in \eqn{[0,1)} (default 0.1).
#' Prune proposals use the negated score of each cherry's current split. The
#' acceptance ratio uses the true soft exposures and the exact forward and
#' reverse proposal probabilities, so the posterior is the same as with
#' `informed = FALSE`; only the proposal law changes. The soft exposures of
#' the two proposed children remain the cost of each move.
#'
#' PCG starts with log-gate proposal factor \eqn{L_0=\mathrm{diag}(sd\_gate)}.
#' At gate step \eqn{t=0,1,\ldots}, draw \eqn{u\sim N(0,I)} and propose
#' \eqn{\log\gamma'=\log\gamma+L_tu}. Robust adaptive Metropolis (RAM)
#' updates the proposal covariance during warm-up using the MH acceptance
#' probability \eqn{A_t}, including rejection:
#' \deqn{L_{t+1}L_{t+1}^{\mathsf T}=
#' L_t\left[I+\eta_t(A_t-a_*)\frac{uu^{\mathsf T}}{u^{\mathsf T}u}\right]
#' L_t^{\mathsf T},\qquad \eta_t=\min\{1,r\,(t+1)^{-\kappa}\}.}
#' Here \eqn{r} is the number of free gates, \eqn{a_*} is `ram_target`, and
#' \eqn{\kappa} is `ram_decay`. Require `0 < ram_target < 1`,
#' `0.5 < ram_decay <= 1`, and integer `0 <= ram_adapt <= burn`.
#' Adaptation stops after `ram_adapt` iterations; `ram_adapt = 0` fixes the
#' initial covariance. `sd_gate` contains positive initial log-scale standard
#' deviations, not the final adapted scales. Shared gates require scalar gate
#' controls. `update_gate = FALSE` keeps gates and proposal covariance fixed,
#' while the intensity and label draws still occur. Reproducible draws use
#' the existing `seed` control. Adaptation is not a convergence diagnostic
#' or an assurance of an efficiency gain.
#'
#' PCG diagnostics include `gate_joint_acceptance` (the fraction of joint gate
#' proposals accepted), its per-chain values, the per-chain mean acceptance
#' probabilities, and `ram` (final covariance and factor per chain, successful
#' update counts, and failed numerical update counts). A numerical RAM failure
#' retains the previous factor; it does not change the MH decision. The legacy
#' per-dimension gate acceptance entries repeat the joint acceptance fraction.
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
#' @references Vihola, M. (2012). Robust adaptive Metropolis algorithm with
#' coerced acceptance rate. \emph{Statistics and Computing}, 22, 997--1008.
#' \doi{10.1007/s11222-011-9269-5}.
#'
#' Zanella, G. (2020). Informed proposals for local MCMC in discrete spaces.
#' \emph{Journal of the American Statistical Association}, 115, 852--865.
#' \doi{10.1080/01621459.2019.1585255}.
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
                    sampler = c("rjmcmc", "irjmcmc", "pcg", "smc", "pgas"),
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
      "partially collapsed Gibbs (`pcg`), ",
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
