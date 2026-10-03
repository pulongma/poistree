# poistree

`poistree` fits Bayesian tree models for Poisson point-process intensities.
Its component API separates the gating rule from the posterior sampler:

```r
fit <- ppt_fit(
  x, region,
  gating = c("soft", "hard"),
  scales = "leaf",
  sampler = c("rjmcmc", "smc", "pgas"),
  ...
)
```

The two models are:

| Model | Gating | Intensity | Available samplers |
|---|---|---|---|
| PPT | hard | terminal leaf | SMC, RJ-MCMC, Particle Gibbs (`pgas` token) |
| S-PPT | soft | terminal leaf | RJ-MCMC, Particle Gibbs with ancestor sampling (`pgas`) |

`scales` is reserved and accepts only `"leaf"`. The default call fits S-PPT
by RJ-MCMC:

```r
sppt <- ppt_fit(
  x, region,
  predict_at = prediction_grid,
  chains = 4, iter = 10000, burn = 2500
)
```

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
constraints remain active in either case. The default minimum child occupancy
is `min_leaf_n = 1` for every model and sampler; larger values can be supplied
as an explicit regularization or computational control.

For the hard PPT, `sampler = "pgas"` selects the exact conditional-SMC
Particle-Gibbs backend without ancestor sampling. For S-PPT it selects
Particle Gibbs with ancestor sampling (PGAS): the conditional SMC runs over
the tree and the latent allocation variables of the observations jointly, one heap node per step,
with the exact one-step-ahead proposal (a Poisson-binomial expectation over
the allocations of the node's observations to its children) for nodes with at most `exact_max`
points and, above that, a Laplace approximation of the same expectation written
as a one-dimensional Beta integral (O(m) per candidate, correct tails).
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

Every backend returns class `ppt` with the same major fields. Use
`ppt_lambda()`, `ppt_predict()`, `ppt_marginal()`, `plot()`, `ppt_summary()`,
`ppt_logLik()`, `ppt_lppd()`, `ppt_integral()`, `ppt_diagnostics()`, and
`ppt_sim()` for downstream work. SMC fits additionally store their log
marginal-likelihood (evidence) estimate in `fit$posterior$log_evidence`,
reported by `ppt_summary()`; the MCMC backends do not estimate it.

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
