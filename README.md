# poistree

`poistree` fits Bayesian tree models for Poisson point-process intensities.
Its component API separates the gating rule, intensity representation, scale
prior, and posterior sampler:

```r
fit <- ppt_fit(
  x, region,
  gating = c("soft", "hard"),
  scales = c("multiscale", "leaf"),
  sampler = c("rjmcmc", "irjmcmc", "smc", "pgas"),
  scale_prior = c("independent", "markov"),
  ...
)
```

The four primary models are:

| Model | Gating | Intensity scales | Available samplers |
|---|---|---|---|
| PPT | hard | terminal leaf | SMC, RJ-MCMC, Particle Gibbs (`pgas` token) |
| S-PPT | soft | terminal leaf | RJ-MCMC |
| MPPT | hard | multiscale | RJ-MCMC; iRJ-MCMC for independent scales |
| S-MPPT | soft | multiscale | RJ-MCMC; iRJ-MCMC for independent scales |

For multiscale models, `scale_prior` selects independent or Markov scales.
The default call fits S-MPPT with an independent scale prior. Setting only
`scales = "leaf"` fits S-PPT:

```r
smppt <- ppt_fit(
  x, region,
  predict_at = prediction_grid,
  chains = 4, iter = 10000, burn = 2500
)

sppt <- ppt_fit(
  x, region,
  scales = "leaf",
  predict_at = prediction_grid,
  chains = 4, iter = 10000, burn = 2500
)
```

Explicit `gating = "hard"` choices retain the hard-gated models:

```r
ppt <- ppt_fit(
  x, region,
  gating = "hard", scales = "leaf", sampler = "smc",
  predict_at = prediction_grid, particles = 1000
)

mppt <- ppt_fit(
  x, region,
  gating = "hard", scales = "multiscale", sampler = "rjmcmc",
  predict_at = prediction_grid,
  chains = 4, iter = 10000, burn = 2500
)
```

For independent-scale MPPT and S-MPPT, `sampler = "irjmcmc"` selects the
informed reversible-jump kernel. It always combines likelihood-informed cut
probabilities with sequential conditional allocation proposals. Their
normalized forward and reverse probabilities enter the Metropolis--Hastings
ratio, so the kernel retains the ordinary posterior target.

The current implementation informs the cut only after the grow/change action,
eligible node, and split dimension have been selected by their ordinary
proposal laws. It is therefore an exact locally informed RJ-MCMC proposal,
not the full-neighborhood informed-importance-tempering (IIT) transition.

```r
fit_irj <- ppt_fit(
  x, region,
  gating = "soft", scales = "multiscale",
  scale_prior = "independent", sampler = "irjmcmc",
  chains = 4, iter = 10000, burn = 2500
)
```

Run `demo("soft-tree-surface", package = "poistree")` for a complete S-MPPT
example that simulates a two-bump Poisson process, fits without `predict_at`,
and evaluates a posterior surface and transect afterward with `ppt_lambda()`.

Hard-leaf SMC and Particle Gibbs impose no aspect-ratio restriction by default
(`max_aspect_ratio = Inf`). A finite value can be supplied when compact cells
are scientifically required; the minimum child width and observation-count
constraints remain active in either case. The default minimum child occupancy
is `min_leaf_n = 1` for every model and sampler; larger values can be supplied
as an explicit regularization or computational control.

For backward compatibility, `sampler = "pgas"` selects the exact
conditional-SMC Particle-Gibbs backend. Ancestor sampling is currently
disabled: its valid tree-specific implementation requires a full suffix
target ratio, rather than only the probability of the next reference action.

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
