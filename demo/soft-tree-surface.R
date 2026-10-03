## S-MPPT surface estimation with post-hoc intensity evaluation.
## This demonstration deliberately does not use `predict_at` in `ppt_fit()`.

needed <- c("ggplot2", "patchwork")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(
    "This demo requires the suggested package(s): ",
    paste(missing, collapse = ", "),
    call. = FALSE
  )
}

library(poistree)
library(ggplot2)
library(patchwork)

set.seed(2026)

## True intensity on [0,1]^2: two smooth bumps with exactly 400 expected
## points. The normalizer accounts for Gaussian mass truncated by the region.
bump1_mass <- diff(pnorm(c(0, 1), mean = 0.3, sd = 0.12))^2
bump2_mass <- diff(pnorm(c(0, 1), mean = 0.75, sd = 0.10)) *
  diff(pnorm(c(0, 1), mean = 0.70, sd = 0.15))
lambda_normalizer <- bump1_mass + 0.6 * bump2_mass

lambda_true <- function(x) {
  400 * (
    dnorm(x[, 1], 0.3, 0.12) * dnorm(x[, 2], 0.3, 0.12) +
      0.6 * dnorm(x[, 1], 0.75, 0.10) * dnorm(x[, 2], 0.70, 0.15)
  ) / lambda_normalizer
}

## Simulate the inhomogeneous Poisson process by thinning.
region <- matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE)
lam_max <- 400 * (
  dnorm(0, 0, 0.12)^2 +
    0.6 * dnorm(0, 0, 0.10) * dnorm(0, 0, 0.15)
) / lambda_normalizer
n_prop <- rpois(1, lam_max)
prop <- cbind(runif(n_prop), runif(n_prop))
keep <- runif(n_prop) < lambda_true(prop) / lam_max
x <- prop[keep, , drop = FALSE]
colnames(x) <- c("x1", "x2")
cat("simulated points:", nrow(x), "\n")

## Fit S-MPPT (soft multiscale PPT with independent scales).
## No `predict_at` is needed.
fit <- ppt_fit(
  x, region,
  gating = "soft", scales = "leaf",
  scale_prior = "independent", sampler = "rjmcmc",
  max_depth = 6, min_leaf_n = 1,
  chains = 4, iter = 3000, burn = 1000, thin = 4,
  cut_candidates = 15, gate = 15, update_gate = TRUE,
  verbose = FALSE
)
print(fit)

## Post-hoc grid evaluation: 90 x 90 lattice and a ggplot2-ready data frame.
surface <- ppt_lambda(fit, n = 90)
grid <- as.matrix(surface[, c("x1", "x2")])
surface$true <- lambda_true(grid)
surface$ci_width <- surface$upper - surface$lower

theme_sci <- theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10),
        legend.position = "right",
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey35"))
fill_scale <- scale_fill_viridis_c(
  name = expression(lambda(x)),
  limits = c(0, max(surface$true, surface$mean))
)

p_true <- ggplot(surface, aes(x1, x2, fill = true)) +
  geom_raster() + fill_scale + coord_equal(expand = FALSE) + theme_sci +
  labs(title = "True intensity", x = expression(x[1]), y = expression(x[2]))

p_mean <- ggplot(surface, aes(x1, x2, fill = mean)) +
  geom_raster() + fill_scale + coord_equal(expand = FALSE) + theme_sci +
  labs(title = "Posterior mean",
       #subtitle = "ppt_lambda(fit, n = 90); no predict_at at fit time",
       x = expression(x[1]), y = expression(x[2])) +
  geom_point(data = as.data.frame(x), aes(x1, x2), inherit.aes = FALSE,
             size = 0.15, alpha = 0.35, colour = "white")

p_ci <- ggplot(surface, aes(x1, x2, fill = ci_width)) +
  geom_raster() +
  scale_fill_viridis_c(name = "95% CI\nwidth", option = "magma") +
  coord_equal(expand = FALSE) + theme_sci +
  labs(title = "Pointwise uncertainty", x = expression(x[1]),
       y = expression(x[2]))

# 1-D transect at x2 = 0.3 straight through the main bump -- again post hoc.
transect <- cbind(seq(0, 1, length.out = 300), 0.3)
colnames(transect) <- c("x1", "x2")
line <- ppt_lambda(fit, at = transect)
line$true <- lambda_true(transect)

p_line <- ggplot(line, aes(x1)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue",
              alpha = 0.25) +
  geom_line(aes(y = mean, colour = "Posterior mean"), linewidth = 0.7) +
  geom_line(aes(y = true, colour = "Truth"), linewidth = 0.6,
            linetype = "22") +
  scale_colour_manual(NULL, values = c("Posterior mean" = "steelblue4",
                                       "Truth" = "black")) +
  theme_sci +
  labs(title = expression("Transect at" ~ x[2] == 0.3),
       subtitle = "shaded: pointwise 95% credible band",
       x = expression(x[1]), y = expression(lambda(x[1], 0.3)))

figure <- (p_true | p_mean) / (p_ci | p_line) +
  plot_annotation(
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  )

print(figure)
