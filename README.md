# poistree

`poistree` fits Bayesian tree models for Poisson point-process intensities.
The public interface separates three model components:

```r
fit <- ppt_fit(
  x,
  region,
  gating = c("hard", "soft"),
  scales = c("leaf", "multiscale"),
  sampler = c("rjmcmc", "smc", "pgas"),
  scale_prior = c("markov", "independent"),
  ...
)
```

The only implemented model is the hard-gated terminal-leaf PPT. It can be
fitted with sequential Monte Carlo, reversible-jump MCMC, or Particle Gibbs
with ancestor sampling:

```r
fit <- ppt_fit(
  x,
  region,
  gating = "hard",
  scales = "leaf",
  sampler = "rjmcmc",
  predict_at = prediction_grid,
  chains = 2, iter = 4000, burn = 1000
)
```

The result has S3 class `ppt`. Use `ppt_predict()`, `plot()`,
`ppt_summary()`, `ppt_logLik()`, `ppt_evidence()`, `ppt_lppd()`,
`ppt_diagnostics()`, and `ppt_sim()`.

The general component selectors remain part of `ppt_fit()` so future models,
including soft-gated or multiscale trees, can be added without changing user
code. Unsupported component combinations currently give an informative error.

```r
# Hard terminal-leaf PPT fitted by PGAS
ppt_pgas <- ppt_fit(
  x, region,
  gating = "hard", scales = "leaf", sampler = "pgas",
  predict_at = prediction_grid,
  particles = 500, iter = 1000, burn = 200
)
```

```r
# Hard terminal-leaf PPT
ppt <- ppt_fit(
  x, region,
  gating = "hard", scales = "leaf", sampler = "smc",
  predict_at = prediction_grid
)
```
