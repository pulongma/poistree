# poistree

`poistree` fits Bayesian tree models for Poisson point-process intensities.
Use it to estimate intensity surfaces, summarize uncertainty, explore marginal
intensity curves, and evaluate predictions on held-out point patterns.

## Implemented models

| Model | Gating | Available samplers |
| --- | --- | --- |
| PPT | Hard | SMC (`smc`), RJ-MCMC (`rjmcmc`) |
| S-PPT | Soft | partially collapsed Gibbs (`pcg`) |

Both models use terminal-leaf intensities. `ppt_fit()` defaults to S-PPT with
RJ-MCMC. For hard PPT, the `pgas` option uses Particle Gibbs without ancestor
sampling.

## Main functions

| Function | Purpose |
| --- | --- |
| `ppt_fit()` | Fit a model over a rectangular observation region. |
| `ppt_fit_quadrature()` | Fit using background covariates and exposure weights. |
| `ppt_sim()` | Simulate point patterns from a supplied intensity function. |
| `ppt_predict()` | Predict intensity summaries at specified locations. |
| `ppt_lambda()` | Obtain intensity surfaces, credible intervals, or posterior draws. |
| `ppt_marginal()` | Summarize intensity along one covariate, averaging over the others. |
| `ppt_lppd()` | Compute a test point pattern's log predictive density. |
| `ppt_diagnostics()` | Inspect sampler diagnostics and produce diagnostic plots. |

## Simulate and fit

The following example generates a point pattern on the unit square and fits
an S-PPT model. The runs below are short demonstrations.

```r
library(poistree)

region <- rbind(c(0, 1), c(0, 1))
intensity <- function(z) {
  40 + 160 * exp(-20 * ((z[, 1] - 0.3)^2 + (z[, 2] - 0.7)^2))
}
x <- ppt_sim(intensity, region, lambda_max = 200, seed = 7)
colnames(x) <- c("x1", "x2")

fit <- ppt_fit(
  x, region, gating = "soft", sampler = "rjmcmc",
  chains = 2, iter = 1000, burn = 300, max_depth = 3,
  seed = 7, verbose = FALSE
)
summary(fit)
```

For soft PCG or hard PPT with SMC:

```r
fit_pcg <- ppt_fit(
  x, region, gating = "soft", sampler = "pcg",
  chains = 2, iter = 1000, burn = 300, max_depth = 3,
  seed = 7, verbose = FALSE
)

fit_ppt <- ppt_fit(
  x, region, gating = "hard", sampler = "smc",
  particles = 200, max_depth = 3, seed = 7
)
```

## Predict and plot intensity

Evaluate the fitted intensity at new locations, or obtain posterior draws for
further analysis:

```r
at <- cbind(x1 = seq(0, 1, length.out = 100), x2 = 0.5)
prediction <- ppt_predict(fit, newdata = at, type = "interval")
head(prediction)

posterior <- ppt_lambda(fit, at = at, type = "draws")
dim(posterior$draws)
```

Create an intensity surface with `ggplot2`:

```r
library(ggplot2)

surface <- ppt_lambda(fit, n = 40)
ggplot(surface, aes(x1, x2, fill = mean)) +
  geom_raster() +
  scale_fill_viridis_c(option = "plasma") +
  coord_equal() +
  labs(fill = "Intensity", title = "Posterior mean intensity") +
  theme_minimal()
```

## Marginal intensity curves

`ppt_marginal()` averages over the fitted ranges of the other covariates.
It returns a curve with posterior mean, median, and pointwise credible limits.

```r
marginal <- ppt_marginal(fit, variable = "x1", level = 0.95)
ggplot(marginal, aes(value, mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() +
  labs(x = "x1", y = "Average intensity") +
  theme_minimal()
```

## Diagnostics and held-out prediction

Inspect the available acceptance rates, tree traces, or particle diagnostics:

```r
diagnostics <- ppt_diagnostics(fit)
print(diagnostics)
ppt_diagnostics(fit, plot = TRUE)
```

Score an independent point pattern observed over the same region and exposure:

```r
test <- ppt_sim(intensity, region, lambda_max = 200, seed = 8)
ppt_lppd(fit, test = test)
```

## Fit with background covariates and exposure weights

`ppt_fit_quadrature()` supports S-PPT with PCG and hard PPT with SMC. Here,
background locations are cell centers on the unit square, and their weights
are cell areas.

```r
centers <- seq(0.025, 0.975, length.out = 20)
background <- as.matrix(expand.grid(x1 = centers, x2 = centers))
weights <- rep(1 / nrow(background), nrow(background))

fit_q <- ppt_fit_quadrature(
  x, region, background = background, weights = weights,
  gating = "soft", sampler = "pcg",
  chains = 2, iter = 1000, burn = 300, max_depth = 3,
  seed = 7, verbose = FALSE
)
summary(fit_q)
ppt_predict(fit_q, newdata = at, type = "mean")
```

For these fits, `ppt_marginal()` describes the covariate-box average, rather
than a marginal weighted by physical exposure.

## Complete demo and help

The included demo simulates a two-bump point process and plots the fitted
surface, uncertainty, and a transect. It requires `ggplot2` and `patchwork`.

```r
demo("soft-tree-surface", package = "poistree")
```

See `?ppt_fit`, `?ppt_fit_quadrature`, `?ppt_controls`, and the individual
function help pages for full documentation.
