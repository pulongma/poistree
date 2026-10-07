# poistree

`poistree` fits Bayesian tree models for Poisson point-process intensities.
Its component API separates the gating rule from the posterior sampler:

```r
fit <- ppt_fit(
  x, region,
  gating = c("soft", "hard"),
  scales = "leaf",
  sampler = c("rjmcmc", "irjmcmc", "pcg", "smc", "pgas"),
  ...
)
```

The two models are:

| Model | Gating | Intensity | Available samplers |
|---|---|---|---|
| PPT | hard | terminal leaf | SMC, RJ-MCMC, informed MH (`irjmcmc`), Particle Gibbs (`pgas` token) |
| S-PPT | soft | terminal leaf | RJ-MCMC, informed MH (`irjmcmc`), partially collapsed Gibbs (`pcg`), Particle Gibbs with ancestor sampling (`pgas`) |

`scales` is reserved and accepts only `"leaf"`. The default call fits S-PPT
by RJ-MCMC:

```r
sppt <- ppt_fit(
  x, region,
  predict_at = prediction_grid,
  chains = 4, iter = 10000, burn = 2500
)
```

MCMC fits show an RcppProgress bar for each chain by default. Each bar counts
completed iterations, including burn-in and iterations discarded by thinning.
Pass `verbose = FALSE` to `ppt_fit()` to suppress bars and chain messages.

Explicit `gating = "hard"` selects PPT:

```r
ppt <- ppt_fit(
  x, region,
  gating = "hard", sampler = "smc",
  predict_at = prediction_grid, particles = 1000
)

ppt_mcmc <- ppt_fit(
  x, region,
  gating = "hard", sampler = "rjmcmc",
  predict_at = prediction_grid,
  chains = 4, iter = 10000, burn = 2500
)
```

Run `demo("soft-tree-surface", package = "poistree")` for a complete S-PPT
example that simulates a two-bump Poisson process, fits without `predict_at`,
and evaluates a posterior surface and transect afterward with `ppt_lambda()`.

Hard-leaf SMC and Particle Gibbs impose no aspect-ratio restriction by default
(`max_aspect_ratio = Inf`). A finite value can be supplied when compact cells
are scientifically required; the minimum child width and observation-count
constraints remain active in either case. All hard samplers and soft RJ-MCMC
informed MH, and PCG accept the positive integer `min_leaf_n` (default 1), which
screens candidate splits by their hard-routed child counts. Soft PGAS does
not accept this argument.

Hard RJ-MCMC and informed MH accept `prediction_draws` to control retained
predictions and tree states. The soft versions accept `thin`, `tree_moves`,
and `change_moves` instead, and reject `prediction_draws`. PCG accepts the
soft RJ-MCMC controls. Soft RJ-MCMC, informed MH, and PCG also accept `gate_family = "logistic"` (the default) or
`"compact"`. Soft PGAS uses logistic gates and does not accept `gate_family`.
All soft backends accept `gate_structure`, `gate`, and gate-prior controls.

All soft backends accept `gate_scale = NULL`, independently of
`gate_structure`. For logistic gates, `NULL` retains the historical
`"root"` scaling: the gate on input `j` is
`plogis(gate[j] * (x[j] - cut) / root_width[j])`.
Use `gate_scale = "node"` to divide by the splitting node's local
width instead, so the same gate parameter gives a sharper transition
in a narrower node. `gate_structure = "dimension"` continues to share
one parameter per input across its nodes; `"shared"` uses one parameter
across all inputs. Compact gates keep their historical node scaling
when `gate_scale` is `NULL` or `"node"`; `"root"` is not supported.
Node-relative logistic path integrals use deterministic adaptive
quadrature. Stored posterior states retain the scale choice, and older
fitted objects continue to use their stored gate mode.
Node scaling for logistic gates is available with `sampler = "pcg"`,
`"rjmcmc"`, and `"irjmcmc"`. Soft PGAS currently supports only root
scaling; it rejects `gate_scale = "node"` before fitting because its
fixed-grid ancestor sampler requires a different derivation for
node-dependent widths.

```r
fit_node <- ppt_fit(
  x, region, gating = "soft", sampler = "pcg",
  gate_structure = "dimension", gate_scale = "node"
)
```

Omitting `gate_scale` keeps existing calls and logistic results unchanged.
The effective setting is recorded in `fit_node$model$gate_scale`,
`fit_node$prior$gate$scale`, and `fit_node$control$gate_scale`.


`max_depth` must be one finite integer. All hard samplers and soft RJ-MCMC
informed MH, and PCG accept 0 through 20; zero gives a root-only tree. Soft PGAS
requires at least 1 and its current reference-tree storage limit permits at
most 19. Dense hard SMC and hard Particle Gibbs additionally limit the
combination of depth and particle count. If a tree-storage limit is exceeded,
reduce `max_depth` or, for those hard backends, `particles`. These guards
protect tree indexing and preallocated tree storage; cache and output memory
remain additional costs. Existing default depths are unchanged.

Every sampler defaults to `cut_candidates = 50`, including hard Particle
Gibbs. For quantile proposals, this counts requested quantile probabilities;
duplicate and inadmissible cuts are removed, so the number of valid cuts can
be smaller. The count must be an integer of at least 1 for hard SMC and hard
Particle Gibbs, and at least 2 for the other backends. Soft PGAS constructs
its grid globally; the RJ-MCMC samplers and
hard particle samplers construct cuts within each node. For soft RJ-MCMC,
informed MH, and PCG, `cut_proposal = "uniform"` retains the existing grid-size floor
of 30 before filtering, so a smaller explicit `cut_candidates` still uses
30 grid locations. `cut_proposal = "data"` uses all admissible data midpoints
and ignores the count. The default quantile and uniform grids request 50
locations. A supplied soft PGAS `cut_grid` overrides its generated quantile grid.

Use `sampler = "irjmcmc"` for cached locally informed Metropolis--Hastings
with either hard or soft gating:

```r
informed <- ppt_fit(
  x, region, gating = "hard", sampler = "irjmcmc",
  chains = 4, iter = 10000, burn = 2500, cut_candidates = 50
)
```

This sampler preserves the corresponding RJ-MCMC target and returns ordinary
posterior draws. For each reversible kernel with base proposal `q0`, neighbor
scores are `q0(s,t) * sqrt(pi(t) * q0(t,s) / (pi(s) * q0(s,t)))`. Their sum
is `Z(s)`. A neighbor is drawn proportionally to these scores and accepted
with probability `min(1, Z(s) / Z(t))`; rejected transitions retain the current
state. This is an MH sampler, and downstream means, intervals, and predictive
densities use ordinary posterior weights.

The hard sampler uses exact collapsed Gamma-Poisson likelihoods over grow,
prune, and terminal-split change neighbors. It caches candidate cuts and local
scores and keeps the current neighborhood after rejection. The soft sampler
uses the augmented tree and allocation state at fixed gates, with exact
Poisson-binomial count recursion to aggregate allocation scores and sample
the selected allocation. This recursion costs O(m²) per candidate for m
allocated observations at the node. Its `tree_moves` and `change_moves` still control
separate grow/prune and change kernels. Allocation Gibbs updates and gate
Metropolis updates remain in the iteration. Exact informed scoring costs more
per transition, so a gain in effective samples per second depends on the data.

Use `sampler = "pcg"` for an allocation-collapsed joint gate update with
robust adaptive Metropolis (RAM) warm-up:

```r
sppt_pcg <- ppt_fit(
  x, region, gating = "soft", sampler = "pcg",
  chains = 4, iter = 10000, burn = 2500,
  sd_gate = 0.07, ram_target = 0.234, ram_decay = 0.7,
  ram_adapt = 2000
)
```

Each iteration draws temporary leaf intensities, updates all free log gates
jointly after integrating out the labels, then redraws every label from its
categorical full conditional, even if the gate proposal was rejected. It next
discards those intensities and performs standard collapsed grow/prune and
change updates. Retained output uses fresh intensities from the updated tree.
Both `logistic` and `compact` gates and `dimension` and `shared` gate structures
are supported. The tree target and candidate rules match soft RJ-MCMC.

RAM adapts the joint proposal covariance toward `ram_target` using the MH
acceptance probability. `ram_decay` must lie in `(0.5, 1]`.
`ram_adapt` is an integer from zero through `burn`, defaulting to
`floor(0.8 * burn)`; the covariance is fixed afterward. `sd_gate` gives its
initial log-scale standard deviations. Use `ram_adapt = 0` for a fixed
proposal or `update_gate = FALSE` for fixed gates. Adaptation does not certify
convergence or improved sampling efficiency. `ppt_diagnostics(sppt_pcg)` exposes
joint gate acceptance and per-chain final proposal covariance, factor, and
adaptation counts in `ram`. PCG uses ordinary posterior weights.

Soft RJ-MCMC, informed MH, and PCG use `cache_geometry = TRUE` by default.
They reuse shared ancestor log-memberships and cache leaf exposures for both
`gate_scale = "root"` and `"node"`. The optimization also supports compact
gates. Set `cache_geometry = FALSE` to use direct calculations for validation.
The cache is confined to a fit/chain; proposed gates are evaluated separately,
and changes to ancestor rules or local widths invalidate affected entries.
This changes computation only: priors, tree proposals, update counts, and
random-number generation are unchanged. In the private quadrature adapter,
background memberships serve both weighted exposure and grid prediction.

The hard PPT SMC sampler (`sampler = "smc"`) runs by default on a shared-path
store (`engine = "shared"`): particles that reach the same box share one stored
node with its candidate cuts and marginal-likelihood scores, scores are computed
from sorted coordinates with binary-search counts, and resampling copies maps of
node references rather than trees. The sampler is identical in law to the
original per-particle implementation (`engine = "dense"`) and two orders of
magnitude faster when `d` is large; see `CLAUDE/shared_path_smc_hard_ppt.pdf`
in the SoftPPT project.

For the hard PPT, `sampler = "pgas"` selects the exact conditional-SMC
Particle-Gibbs backend without ancestor sampling. For S-PPT it selects
Particle Gibbs with ancestor sampling (PGAS): the conditional SMC runs over
the tree and the latent allocation variables of the observations jointly, one heap node per step,
with the exact one-step-ahead proposal (a Poisson-binomial expectation over
the allocations of the node's observations to its children) for nodes with at most `exact_max`
points and, above that, a Laplace approximation of the same expectation written
as a one-dimensional Beta integral (O(m) per candidate).
For exact-scored nodes, logistic cut locations share one Poisson-binomial
recursion per node and coordinate. Each candidate count distribution is
obtained by an exact exponential tilt of that reference distribution,
reducing allocation-scoring work from O(d M m^2) to O(d m^2 + d M m)
for m allocated observations and M cuts per coordinate. Exposures remain
cut-specific. Conditional allocation tables are built lazily once per selected
coordinate and shared across its cuts. The exact importance ratio cancels the
allocation proposal, so a forced reference allocation needs no table or replay.
Reference count coefficients and count normalizers are reused for sampling;
zero/all-left counts bypass allocation tables. Root scores and lazy coordinate
tables persist across sweeps, invalidating only coordinates whose gate changes.
The default large-node proposal remains `proposal_score = "laplace"`,
and `exact_max` remains 150. An optional `proposal_score = "hard"` uses
cumulative hard-split counts and interval exposures under the current soft
parent path to rank all cuts. It borrows GS-BART's reuse of additive summaries;
it is a Gamma-Poisson proposal surrogate, not a change to the soft model or
an implementation of GS-BART's importance-tempering sampler. Its action scores
use `proposal_temperature = 0.5` and a prior mixture `proposal_defensive = 0.1`.
These controls apply only above `exact_max`; the original `defensive` remains
the mixture control for exact/Laplace scoring. True logistic gates and child
exposures are computed lazily for selected or reference splits, with the full
importance correction. The full candidate grid is retained, with 50 requested
cuts per coordinate by default (duplicate or invalid cuts are removed).
Count preparation costs O(d (m + M)), after fixed-bin preprocessing; path
exposure integration is additional work. Hard surrogates can be inaccurate
for broad soft gates, so compare effective samples per second before choosing
this option. See `inst/benchmarks/pgas-hard-proposal-benchmark.R`.
Above `exact_max` the allocation of the node's points is proposed either by
sequential imputation or through auxiliary child rates (`allocation`), with
exact importance weights in both cases.
Conditional multinomial resampling with ancestor sampling takes place after
every heap node at which a particle advanced (`resampling = "node"`) or after
every tree level (`"level"`), optionally only when the ESS falls below
`ess_threshold * particles`. Ancestor sampling is made valid by
working on a complete decision tree with a fixed per-input cut grid (global
quantiles by default, or a user-supplied `cut_grid`), so that the full
suffix target ratio is defined for every particle. The gating parameters
(one per input by default, `gate_structure = "dimension"`) are updated by a
systematic Metropolis--Hastings scan and the allocation variables by Gibbs
sweeps inside the cycle.
The particle system uses a shared-path store in the spirit of the `ResTree`
package: particles hold handles to path nodes (geometric path plus
the observations allocated to the node), the candidate expansion of each distinct
frontier node is computed once per level and shared by every particle holding
it, and resampling moves or copies handle maps rather than trees
(`fit$diagnostics$expanded_nodes` reports the distinct expansions per sweep).

```r
sppt_pg <- ppt_fit(
  x, region,
  gating = "soft", scales = "leaf", sampler = "pgas",
  particles = 30, iter = 1000, burn = 200, max_depth = 6
)
```

Hard splits consistently route `x < cut` to the left child and ties to the
right child on every dimension. Points on the fitted region's outer boundary
remain included. This convention applies to fitting, saved predictions, and
evaluation of retained states.

Every backend returns class `ppt` with the same major fields. Use
`ppt_lambda()`, `ppt_predict()`, `ppt_marginal()`, `plot()`, `ppt_summary()`,
`ppt_logLik()`, `ppt_lppd()`, `ppt_integral()`, `ppt_diagnostics()`, and
`ppt_sim()` for downstream work. For SMC, `posterior$log_relative_normalizer`
stores the log estimate of the target normalizer divided by the unsplit root
score. With a proper intensity prior (`b > 0`),
`posterior$log_target_normalizer` adds the normalized Gamma-Poisson root
factor; `ppt_summary()` reports it explicitly as a target normalizer.
With `b = 0`, only the relative quantity is available, under the package's
formal improper-prior leaf-score convention.

`posterior$log_evidence` is `NA`: the current tree-prior weights assign each
valid split on dimension `j` weight `rho / (d * K_j)` and do not renormalize
over dimensions without valid cuts. The sum of these tree weights is not
computed, so even with `b > 0` the target normalizer is not absolute Bayesian
evidence. The MCMC backends do not estimate either normalizer. For backward
compatibility, diagnostics named `log_evidence_increment` and
`log_evidence_running` still contain increments and cumulative sums of the
**relative** log normalizer.

The hard backends retain their existing model specifications: RJ-MCMC
defaults to split parameter `alpha = 0.95`, whereas SMC and Particle Gibbs
use `0.5`. Their candidate-cut and cell-geometry rules also differ. Selecting
the same gating mechanism does not by itself make these samplers target the
same posterior. The shared and dense **SMC** engines do share the same target
for identical controls, including `cut_candidates` and `max_depth = 0`.

The posterior intensity of any fit can be evaluated after fitting at
arbitrary locations: `ppt_lambda()` returns pointwise summaries over an
automatic grid (ready for `ggplot2` surface or curve plots) or the raw
locations-by-draws matrix, `ppt_predict()` accepts any `newdata` inside the
fitted region, and `ppt_lppd(fit, test = pattern)` scores held-out point
patterns. Supplying `predict_at` at fit time remains available as a
precomputed fast path.

`ppt_marginal()` integrates each retained posterior state by default: exactly
over axis-aligned boxes for hard gates and using the fitted compact/logistic
path integrals for soft gates. A user-supplied uniform reference design remains
available as an explicit approximation; see `?ppt_marginal`.

### Physical-domain quadrature and informed PCG

`ppt_fit_quadrature()` integrates intensity over supplied background covariates
and physical exposure weights. In version 0.5.7 this API is part of the main
package as well as the real-data adapter. For example:

```r
fit <- ppt_fit_quadrature(
  x = X_train, region = region,
  background = Q, weights = q_weights,
  gating = "soft", sampler = "pcg",
  informed = TRUE,
  proposal_temperature = 0.5,
  proposal_defensive = 0.1,
  gate_scale = "node"
)
fit$control$informed
```

Other controls, including chains, iterations, seed, and prediction locations,
are passed through to `ppt_fit()`. `informed = FALSE` remains the default.
Both root- and node-scaled gates are supported. The wrapper supports soft PCG
and hard shared SMC; informed RJ-MCMC and PGAS are separate algorithms and are
not enabled by this quadrature wrapper. `ppt_integral()` and `ppt_lppd()` use
the saved weighted integrals. `ppt_marginal()` still integrates over the
covariate box and is not a physical-domain marginal for these fits.

## Prediction performance (0.5.2)

Post-hoc soft prediction shares ancestor calculations across leaves, for root
and node logistic scales and compact gates. It needs only one scalar per
active node as temporary working memory. Existing saved fits are supported.
`ppt_lambda()` and `ppt_lppd()` reuse retained prediction draws for matching
locations, including reordered or repeated rows, and evaluate only missing
locations. `ppt_predict(type = "mean")` skips unrequested quantiles; intensity
summary quantiles share one sort while keeping the same weighted definition.
These changes require no refitting or changes to chain settings.

### Fitting performance (0.5.3)

Soft RJ-MCMC and PCG cache candidate cuts and split support for unchanged
node paths. The cache is local to each chain, fills axes only as needed, and
discards obsolete paths after every tree proposal. PCG also reuses the event
log-normalizers from the selected gate evaluation when restoring allocations.
These optimizations apply automatically to root- and node-scaled logistic
gates and compact gates, including the private quadrature adapter. Candidate
ordering, acceptance probabilities and random draws are unchanged; no new
settings or changes to chain length are required. The informed sampler keeps
its existing candidate cache.

### Membership reuse during fitting (0.5.4)

Soft RJ-MCMC and PCG grow/change proposals reuse cached parent log-memberships
and evaluate only the proposed split for affected observations. With the
private quadrature adapter, exact matching training rows also reuse background
memberships; repeated events retain their individual labels and counts.
Unmatched rows use direct gate calculations, without approximate matching.
The immutable row mapping is shared by gate proposals within each chain.
Both optimizations support root/node logistic and compact gates, preserve
random-number order and the existing probability calculations, and are enabled
by the existing default `cache_geometry = TRUE`.
