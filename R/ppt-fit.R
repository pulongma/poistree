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
#' @param x Finite numeric matrix with \eqn{n \ge 1} event locations in rows
#'   and \eqn{d \ge 1} coordinates or covariates in columns. Every row must lie
#'   within `region`. Matrix-like inputs are converted to a numeric matrix.
#' @param region Finite numeric \eqn{d\times 2} matrix, with lower bounds in
#'   column 1 and strictly larger upper bounds in column 2. Rows must follow
#'   the column order of `x`. Ordinary fits integrate intensity over this box;
#'   use [ppt_fit_quadrature()] for physical-domain integration over covariates.
#' @param gating Character scalar selecting the split mechanism:
#'   \describe{
#'     \item{`"soft"`}{S-PPT with probabilistic routing; the default.}
#'     \item{`"hard"`}{PPT with deterministic routing at each split.}
#'   }
#' @param scales Character scalar. Only `"leaf"` is supported: each terminal
#'   leaf carries an intensity parameter. Other values are rejected.
#' @param sampler Character scalar selecting posterior computation:
#'   \describe{
#'     \item{`"rjmcmc"`}{Grow/prune/change RJ-MCMC; the default. Both gates.}
#'     \item{`"irjmcmc"`}{Cached, locally informed Metropolis--Hastings tree
#'       updates with exact neighborhood scoring. Both gates.}
#'     \item{`"pcg"`}{Partially collapsed Gibbs with joint robust adaptive
#'       Metropolis gate updates. Soft gates only.}
#'     \item{`"smc"`}{Sequential Monte Carlo. Hard gates only.}
#'     \item{`"pgas"`}{Particle Gibbs: ancestor sampling for soft gates;
#'       conditional SMC without ancestor sampling for hard gates.}
#'   }
#' @param ... Named controls for the selected backend. See \link{ppt_controls}
#'   for every supported control, its default and admissible range, and an
#'   organized reference for priors, tree support, chains, gate parameters,
#'   SMC, and PGAS proposals. Controls unsupported by the selected backend
#'   are not silently ignored. All backends accept `predict_at`, `test`,
#'   `a`, `b`, `max_depth`, `cut_candidates`, and `seed`.
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
#' @section Informed Metropolis--Hastings:
#' With `sampler = "irjmcmc"`, the current RJ-MCMC target and admissible tree
#' support are retained. Each reversible tree kernel uses square-root
#' balanced neighbor weights. For target density \eqn{p}, base proposal
#' \eqn{q_0}, and neighboring states \eqn{s,t}, define
#' \deqn{r(s,t)=\frac{p(t)q_0(t,s)}{p(s)q_0(s,t)},\qquad
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
#' @section Partially collapsed Gibbs:
#' With `sampler = "pcg"`, each iteration has the following order:
#' \enumerate{
#'   \item Draw temporary leaf intensities conditional on the current tree,
#'     labels, and gates.
#'   \item Propose all free log gates jointly, with labels integrated out.
#'   \item Restore each observation label by an independent categorical draw
#'     conditional on the resulting gates and temporary intensities, including
#'     after a rejected gate proposal.
#'   \item Discard the temporary intensities and run collapsed grow/prune and
#'     change tree updates.
#'   \item Draw fresh intensities when retaining output.
#' }
#' The tree target and candidate rules are those of soft RJ-MCMC.
#'
#' With `sampler = "pcg"` and `informed = TRUE`, the grow, prune, and change
#' proposals are informed by a hard surrogate instead of uniform draws. Every
#' candidate split of a node is scored once, from the node's hard-routed
#' counts and the box exposures of its two children (the one-step lookahead
#' score of the hard SMC, with the same tree prior), and the proposal mixes a
#' tempered draw \eqn{\propto \exp(	au\,	ilde s)} with a uniform draw over
#' all candidates: `proposal_temperature` is \eqn{	au\in(0,1]} (default 0.5,
#' which takes the square root of the untempered score weight) and
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
#' @section Candidate support:
#' The default candidate count is shared by all samplers, but their existing
#' cut-generation rules and tree priors remain different. In particular,
#' soft PGAS uses a global cut grid and RJ-MCMC uses node-specific candidates.
#' The same candidate count does not make different sampler families target
#' the same tree model. Explicit `cut_candidates` values follow the
#' proposal-specific rules in \link{ppt_controls}; a supplied soft PGAS `cut_grid` overrides
#' its automatically generated quantile grid.
#'
#' @return An object of S3 class `ppt`, containing:
#'   \describe{
#'     \item{`model` and `data`}{Model/sampler identifiers, event locations,
#'       region, and optional test locations.}
#'     \item{`prediction`}{Stored prediction locations, posterior intensity
#'       draws, mean, median, and pointwise 95 percent intervals.}
#'     \item{`posterior`}{Saved tree states and posterior summaries, including
#'       `integrated_intensity_draws`, `mean_integrated_intensity`, and
#'       `mean_log_likelihood`. The integrals use the observation box for
#'       [ppt_fit()] or the weighted background for [ppt_fit_quadrature()].
#'       Unavailable evidence summaries are `NA`.}
#'     \item{`diagnostics`}{Sampler-specific traces, acceptance rates, or
#'       particle diagnostics.}
#'     \item{`prior`, `control`, `call`, and `backend`}{Prior settings,
#'       resolved fitting controls, the original call, and backend identifier.}
#'   }
#' @seealso \link{ppt_controls}, [ppt_predict()], [ppt_lambda()],
#'   \link{summary.ppt}, [ppt_diagnostics()], [ppt_fit_quadrature()]
#' @md
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
#' summary(fit)
#'
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
