#' Fitting controls for Poisson point-process trees
#'
#' Reference for named arguments passed through `...` to [ppt_fit()].
#' Supply only controls supported by the selected `gating` and `sampler`.
#' In this page, "local soft samplers" means soft `rjmcmc`, `irjmcmc`, and
#' `pcg`; "hard local samplers" means hard `rjmcmc` and `irjmcmc`.
#' All scalar numeric controls should be finite unless explicitly stated.
#'
#' @section Data, intensity prior, and reproducibility:
#' These controls are accepted by every backend.
#' \describe{
#'   \item{`predict_at`}{Numeric matrix with `ncol(x)` columns, default `x`.
#'     Rows specify locations at which posterior intensity draws are stored.
#'     Values must be finite and inside `region`; zero rows are allowed.
#'     Saved tree states also support later evaluation with [ppt_lambda()].}
#'   \item{`test`}{Optional matrix of held-out event locations, with the same
#'     columns and domain as `x`. Default `NULL`; zero rows are allowed.
#'     These points are used for posterior predictive scoring, not fitting.
#'     The score treats them as a point pattern with the same observation
#'     exposure as the fitted data; use [ppt_lppd()] to specify another scale.}
#'   \item{`a`}{Positive scalar leaf-intensity prior shape, default `0.5`.
#'     With positive `b`, leaf intensities have a Gamma(shape = `a`,
#'     rate = `b`) prior, with mean `a / b`.}
#'   \item{`b`}{Leaf-intensity prior rate. Soft backends require a positive
#'     scalar; their default `NULL` sets
#'     `a * prod(region[, 2] - region[, 1]) / nrow(x)`.
#'     Hard backends accept a nonnegative scalar, default `0`.
#'     Their `b = 0` convention uses the unnormalized leaf prior density
#'     \eqn{p(\lambda)\propto\lambda^{a-1}} and leaf score
#'     \eqn{\Gamma(a+m)/H^{a+m}}, for count \eqn{m} and exposure \eqn{H};
#'     it does not define absolute Bayesian evidence.
#'     The quadrature wrapper instead defaults to
#'     `a * sum(weights) / nrow(x)`.}
#'   \item{`seed`}{Seed passed to [set.seed()] once before fitting,
#'     default `1L`. Subsequent chains continue the resulting RNG stream.
#'     For reproducibility, use a non-missing integer seed and the same R,
#'     package, and RNG settings.}
#' }
#'
#' @section Tree prior and candidate support:
#' \describe{
#'   \item{`max_depth`}{Integer maximum leaf depth, with root depth zero.
#'     Default `8L`, except soft PGAS uses `6L`. Values `0:20` are accepted
#'     by hard backends and local soft samplers; `0` forces a root-only tree.
#'     Soft PGAS requires `1:19`. Dense hard SMC and hard Particle Gibbs
#'     additionally require
#'     \eqn{(2P+3)(2^{D+1}-1)\le 2^{24}}, where \eqn{P} is `particles`
#'     and \eqn{D} is `max_depth`. These are tree-storage guards, not total
#'     memory limits.}
#'   \item{`min_leaf_n`}{Positive integer, default `1L`, accepted by every
#'     backend except soft PGAS. A candidate split must place at least this
#'     many observations in each child under deterministic routing.
#'     For soft gates this restricts candidate support, not latent label counts.}
#'   \item{`cut_candidates`}{Integer requested candidate count, default `50L`
#'     for every backend. At least `1L` for hard SMC and hard Particle Gibbs;
#'     at least `2L` otherwise. Node-specific quantile grids request this
#'     many probability levels before duplicate and invalid cuts are removed.
#'     The local soft `"uniform"` rule uses at least 30 positions and its
#'     `"data"` rule ignores this count. Soft PGAS uses a global grid per
#'     coordinate unless `cut_grid` is supplied. Matching candidate counts
#'     across backend families does not imply matching tree priors or support.}
#'   \item{`alpha`}{Scalar in `(0, 1)`, default `0.95`, accepted by hard local
#'     samplers and every soft backend. Base split probability in
#'     \eqn{P(\mathrm{split}\mid h)=\alpha(1+h)^{-\eta}}, before depth and
#'     candidate-support restrictions, where \eqn{h} is node depth. The
#'     computed probability is clipped to `[1e-12, 1 - 1e-12]` for numerical
#'     stability. Hard SMC
#'     and hard Particle Gibbs fix the base value at `0.5`; they do not accept
#'     this control.}
#'   \item{`eta`}{Nonnegative scalar depth penalty, default `2`, accepted by
#'     the same backends as `alpha`. Larger values reduce split probabilities
#'     more rapidly with depth. Hard SMC and hard Particle Gibbs fix it at `2`.}
#'   \item{`cut_proposal`}{Character scalar for local soft samplers only:
#'     \itemize{
#'       \item `"quantile"` (default): node-specific empirical cut locations
#'         at probabilities from 0.05 to 0.95.
#'       \item `"uniform"`: equally spaced locations inside the node's box,
#'         using `max(30, cut_candidates)` positions before count filtering.
#'       \item `"data"`: midpoints between adjacent distinct observed values
#'         in the node, subject to boundary and child-count restrictions.
#'     }
#'     Node observations here are those routed by the hard split rules.
#'     Quantile and data rules exclude cuts within `1e-3` of a node boundary,
#'     in the units of the corresponding coordinate. All rules enforce
#'     `min_leaf_n`.}
#'   \item{`max_aspect_ratio`}{Numeric scalar in `[1, Inf]`, default `Inf`;
#'     hard SMC and hard Particle Gibbs only. Maximum ratio of longest to
#'     shortest side of each candidate child box, in coordinate units.
#'     `Inf` removes this ratio restriction; other shape and count checks
#'     still apply.}
#' }
#'
#' @section MCMC length and output:
#' All samplers other than hard SMC accept `chains`, `iter`, `burn`, and
#' `verbose`. Defaults differ by backend:
#' \tabular{lrrrr}{
#'   Backend \tab chains \tab iter \tab burn \tab thin\cr
#'   Hard RJ-MCMC / informed MH \tab 4 \tab 4000 \tab 1000 \tab --\cr
#'   Soft local samplers \tab 4 \tab 10000 \tab 2500 \tab 3\cr
#'   Hard Particle Gibbs \tab 1 \tab 500 \tab 100 \tab --\cr
#'   Soft PGAS \tab 1 \tab 1000 \tab 200 \tab 1
#' }
#' \describe{
#'   \item{`chains`}{Positive integer number of chains, run sequentially.}
#'   \item{`iter`}{Integer total iterations per chain, including warm-up;
#'     must be greater than `burn`.}
#'   \item{`burn`}{Nonnegative integer number of initial iterations discarded.}
#'   \item{`thin`}{Positive integer retention interval for soft backends only.
#'     Output is saved at one-based iteration numbers
#'     `burn + 1`, `burn + 1 + thin`, and so on, not exceeding `iter`.
#'     Each chain retains `ceiling((iter - burn) / thin)` states.}
#'   \item{`prediction_draws`}{Positive integer requested number of retained
#'     intensity draws and tree states per chain; default `300L`, hard local
#'     samplers only. With `K = iter - burn`, retention uses interval
#'     `max(1, floor(K / min(prediction_draws, K)))`, starting at the first
#'     post-warm-up state. Thus the actual count can exceed the requested
#'     value when the interval does not divide `K`. Structural diagnostics
#'     are recorded at every post-warm-up iteration.}
#'   \item{`verbose`}{Logical scalar, default `TRUE`. Displays chain messages
#'     and an RcppProgress bar counting completed iterations, including
#'     warm-up and states discarded by thinning. `FALSE` suppresses them.
#'     Not accepted by hard SMC.}
#' }
#'
#' @section Soft gates and their priors:
#' Every soft backend accepts the following controls except where noted.
#' The controls `gate`, `a_gate`, `b_gate`, `sd_gate`, and `gate_min` must be
#' numeric scalars or vectors of length `ncol(x)` for dimension-specific
#' gates. A scalar is recycled across coordinates. Shared gates require
#' scalar controls. All values must be finite.
#' \describe{
#'   \item{`gate_structure`}{Character scalar:
#'     \itemize{
#'       \item `"dimension"` (default): one gate parameter per coordinate,
#'         shared across all splits on that coordinate.
#'       \item `"shared"`: one common gate parameter for every split.
#'     }}
#'   \item{`gate`}{Positive initial gate sharpness, default `12`; every entry
#'     must exceed its `gate_min`. Larger values sharpen the transition at
#'     a fixed cut and width. If `update_gate = FALSE`, this is the fixed value.}
#'   \item{`a_gate`}{Positive gate-prior shape, default `36`. Each free gate
#'     \eqn{\gamma} has density proportional to
#'     \eqn{\gamma^{a_g-1}\exp(-b_g\gamma)\,1\{\gamma>g_{\min}\}}, where
#'     \eqn{a_g}, \eqn{b_g}, and \eqn{g_{\min}} are the corresponding
#'     `a_gate`, `b_gate`, and `gate_min`. Dimension-specific gates have
#'     independent priors; a shared gate contributes one prior factor.}
#'   \item{`b_gate`}{Positive gate-prior rate, default `3`. With
#'     `gate_min = 0`, the prior mean is `a_gate / b_gate`. A positive lower
#'     bound gives a truncated Gamma distribution, whose mean differs.}
#'   \item{`gate_min`}{Nonnegative lower truncation bound, default `0`.
#'     Gate proposals at or below the bound have zero target density.}
#'   \item{`sd_gate`}{Positive log-gate proposal standard deviation,
#'     default `0.07`. RJ-MCMC, informed MH, and PGAS use coordinate-wise
#'     Gaussian random walks in log gates, or one random walk for a shared
#'     gate. PCG uses these as the initial diagonal proposal standard
#'     deviations before joint RAM adaptation.}
#'   \item{`update_gate`}{Logical scalar, default `TRUE`. Enables gate
#'     Metropolis updates. `FALSE` fixes gates at `gate` and disables PCG
#'     adaptation; tree and allocation updates still run.}
#'   \item{`gate_family`}{Character scalar for local soft samplers only:
#'     \itemize{
#'       \item `"logistic"` (default): right-routing probability
#'         \eqn{\{1+\exp[-\gamma_j(x_j-c)/w]\}^{-1}}.
#'       \item `"compact"`: right-routing probability \eqn{t^2(3-2t)} for
#'         \eqn{0<t<1}, with values zero below and one above this interval;
#'         \eqn{t=(x_j-c+h)/(2h)}, \eqn{h=w/\gamma_j}, and \eqn{w} is the
#'         splitting node's width.
#'     }
#'     Soft PGAS uses logistic gates and does not accept `gate_family`.}
#'   \item{`gate_scale`}{`NULL` (default) or a character scalar:
#'     \itemize{
#'       \item `"root"`: logistic gates use `region[j, 2] - region[j, 1]`
#'         for \eqn{w}, regardless of node depth.
#'       \item `"node"`: use the splitting node's width for \eqn{w}.
#'         Supported by local soft samplers; soft PGAS rejects this value.
#'       \item `NULL`: resolves to `"root"` for logistic gates and `"node"`
#'         for compact gates. Compact gates reject `"root"`.
#'     }
#'     This controls width scaling independently of parameter sharing.
#'     Node-scaled logistic exposures use deterministic adaptive integration.}
#' }
#'
#' @section Local soft tree updates and geometry:
#' These controls are accepted by soft RJ-MCMC, informed MH, and PCG.
#' \describe{
#'   \item{`tree_moves`}{Nonnegative integer number of grow/prune tree kernels
#'     per iteration, default `3L`. Zero disables these kernels.}
#'   \item{`change_moves`}{Nonnegative integer number of terminal-split change
#'     kernels per iteration, default `8L`. Zero disables these kernels.}
#'   \item{`cache_geometry`}{Logical scalar, default `TRUE`. Reuses unchanged
#'     leaf exposures and ancestor routing factors within a fit. Supports
#'     logistic root/node scaling and compact gates. `FALSE` uses uncached
#'     calculations for numerical comparisons. No saved files are required;
#'     rounding can differ because multiplication is shared across leaves.}
#' }
#'
#' @section Additional PCG controls:
#' Accepted only with `gating = "soft", sampler = "pcg"`.
#' \describe{
#'   \item{`ram_target`}{Scalar in `(0, 1)`, default `0.234`. Target acceptance
#'     rate for the joint log-gate robust adaptive Metropolis (RAM) proposal.}
#'   \item{`ram_decay`}{Scalar in `(0.5, 1]`, default `0.7`. Exponent
#'     \eqn{\kappa} in the adaptation gain
#'     \eqn{\eta_t=\min\{1,r(t+1)^{-\kappa}\}}, where \eqn{r} is the number
#'     of free gate parameters and \eqn{t=0,1,\ldots} indexes gate steps.}
#'   \item{`ram_adapt`}{Integer in `[0, burn]`; default `NULL` resolves to
#'     `floor(0.8 * burn)`. Number of initial iterations in which the RAM
#'     factor can adapt. Zero fixes the initial proposal covariance.}
#'   \item{`informed`}{Logical scalar, default `FALSE`. With `TRUE`, propose
#'     grow/prune/change moves using a cached hard-tree surrogate and the
#'     proposal controls below. Acceptance still uses the soft target and
#'     forward/reverse proposal probabilities. This is distinct from the
#'     exact-neighborhood sampler selected by `sampler = "irjmcmc"`.}
#'   \item{`proposal_temperature`}{Scalar in `(0, 1]`, default `0.5`.
#'     For `informed = TRUE`, multiplies the hard surrogate log score
#'     (including its prior terms) before exponentiating. Smaller values
#'     flatten the informed weights; `0.5` takes their square root.}
#'   \item{`proposal_defensive`}{Scalar in `[0, 1)`, default `0.1`.
#'     For `informed = TRUE`, weight of a uniform proposal over the
#'     applicable candidate moves, mixed with the informed proposal.
#'     This is a uniform mixture, unlike the prior-action mixture in PGAS.
#'     Both proposal controls are unused when `informed = FALSE`.}
#' }
#'
#' @section Hard SMC and Particle Gibbs:
#' \describe{
#'   \item{`particles`}{Integer at least 2. Defaults: `1000L` for hard SMC,
#'     `500L` for hard Particle Gibbs, and `50L` for soft PGAS. Particle Gibbs
#'     uses this many particles per conditional-SMC sweep, including the
#'     retained reference particle.}
#'   \item{`engine`}{Character scalar for hard SMC only:
#'     \itemize{
#'       \item `"shared"` (default): particles visiting the same box share
#'         its candidate geometry and scores.
#'       \item `"dense"`: separate preallocated tree storage per particle.
#'     }
#'     Both implementations use the hard SMC target and proposal rules.
#'     The quadrature wrapper requires `"shared"`.}
#'   \item{`resample_thresh`}{Scalar in `(0, 1]`, default `0.5`, hard SMC
#'     only. Resample when particle ESS is below
#'     `resample_thresh * particles`. Hard Particle Gibbs does not expose
#'     this control. Soft PGAS instead uses `ess_threshold`.}
#' }
#'
#' @section Soft PGAS proposal and resampling:
#' These controls are accepted only with `gating = "soft", sampler = "pgas"`.
#' \describe{
#'   \item{`cut_grid`}{Default `NULL`, or a list of `ncol(x)` nonempty numeric
#'     vectors of finite cuts strictly inside the respective root interval.
#'     Supplied vectors are sorted and deduplicated and override automatic
#'     cut generation. The automatic grid uses type-1 sample quantiles at
#'     `cut_candidates` probabilities from 0.05 to 0.95; cuts within 0.001
#'     times root width of a boundary are removed. A coordinate with no
#'     remaining cut uses its region midpoint. The grid stays fixed across
#'     nodes and iterations.}
#'   \item{`label_sweeps`}{Nonnegative integer, default `1L`. Additional Gibbs
#'     sweeps over observation labels per iteration after the conditional-SMC
#'     tree/allocation update. Zero omits these additional sweeps.}
#'   \item{`ancestor_sampling`}{Logical scalar, default `TRUE`. Enables
#'     ancestor sampling for the retained reference trajectory.
#'     `FALSE` uses conditional SMC without ancestor sampling.}
#'   \item{`exact_max`}{Nonnegative integer, default `150L`. Largest allocated
#'     node count at which candidate scores and child allocations use the
#'     exact Poisson-binomial recursion. Above it, `proposal_score` chooses
#'     the action scorer and `allocation` chooses the allocation proposal.
#'     Zero still allows exact scoring of empty nodes.}
#'   \item{`proposal_score`}{Character scalar for nodes above `exact_max`:
#'     \itemize{
#'       \item `"laplace"` (default): Laplace approximation to the integrated
#'         child-allocation score.
#'       \item `"hard"`: surrogate using cumulative hard-split counts and
#'         hard-threshold exposures under the current soft parent path.
#'         True logistic gates are evaluated lazily; selected splits use
#'         true soft exposures and the full importance correction.
#'     }
#'     Exact-scored nodes are unaffected. Runtime and mixing gains are
#'     data-dependent.}
#'   \item{`proposal_temperature`}{Scalar in `(0, 1]`, default `0.5`.
#'     For `proposal_score = "hard"`, tempers only the surrogate likelihood
#'     ratio: the informed action weight is
#'     \eqn{p_0(A)\widetilde L(A)^{\tau}}, with \eqn{\tau} this control.
#'     The prior action weight \eqn{p_0} is not tempered. Smaller values
#'     reduce concentration on the highest surrogate scores.}
#'   \item{`proposal_defensive`}{Scalar in `(0, 1)`, default `0.1`.
#'     For the hard scorer, the action proposal is
#'     \eqn{(1-\epsilon)q_{\mathrm{informed}}+\epsilon p_0}, with
#'     \eqn{\epsilon} this control. Thus 0.1 assigns 10 percent of the mixture
#'     to the prior action distribution. This and `proposal_temperature`
#'     affect only hard-scored large nodes, but their ranges are validated
#'     for every soft PGAS call. They do not change the posterior target.}
#'   \item{`defensive`}{Scalar in `[0, 1)`, default `0`. Prior-action mixture
#'     weight for exact and Laplace scoring. It is separate from the hard
#'     scorer's `proposal_defensive`.}
#'   \item{`allocation`}{Character scalar for nodes above `exact_max`:
#'     \itemize{
#'       \item `"sequential"` (default): allocate observations to children
#'         one at a time using sequential imputation.
#'       \item `"rates"`: propose auxiliary child rates and allocate
#'         observations independently conditional on those rates.
#'     }
#'     Both branches include their corresponding importance correction.}
#'   \item{`resampling`}{Character scalar determining resampling opportunities:
#'     \itemize{
#'       \item `"node"` (default): after each heap node at which at least
#'         one particle advanced, except the final heap position.
#'       \item `"level"`: after each tree level except the final level.
#'     }
#'     The final tree is selected using the final particle weights.}
#'   \item{`ess_threshold`}{Scalar in `(0, 1]`, default `1`. At a resampling
#'     opportunity, resample if ESS is at most
#'     `ess_threshold * particles`. The default therefore resamples at
#'     every such opportunity.}
#' }
#'
#' @seealso [ppt_fit()], [ppt_fit_quadrature()], [ppt_diagnostics()],
#'   \link{summary.ppt}
#' @name ppt_controls
#' @md
NULL
