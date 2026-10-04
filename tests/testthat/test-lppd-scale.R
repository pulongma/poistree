test_that("scaled and plug-in lppd match direct computations", {
  set.seed(3)
  region <- matrix(c(0, 1), 2, 2, byrow = TRUE)
  x <- matrix(runif(400), ncol = 2)
  test <- matrix(runif(100), ncol = 2)
  r <- 0.25
  fits <- list(
    ppt_fit(x, region, gating = "hard", sampler = "smc", test = test,
            a = 0.5, b = 0.5 / nrow(x), particles = 20, max_depth = 3, seed = 1),
    ppt_fit(x, region, gating = "soft", sampler = "rjmcmc", test = test,
            a = 0.5, b = 0.5 / nrow(x), chains = 1, iter = 60, burn = 10,
            thin = 1, max_depth = 3, seed = 1, verbose = FALSE)
  )
  for (fit in fits) {
    d <- ppt_lambda(fit, at = test, type = "draws")
    lp <- colSums(log(r * d$draws)) - r * ppt_integral(fit, "draws")
    posterior <- max(lp) + log(sum(d$weights * exp(lp - max(lp))))
    plugin <- sum(log(r * ppt_lambda(fit, at = test)$mean)) - r * ppt_integral(fit)
    expect_equal(as.numeric(ppt_lppd(fit, test)), fit$posterior$lppd, tolerance = 1e-8)
    expect_equal(as.numeric(ppt_lppd(fit, scale = r)), posterior, tolerance = 1e-8)
    expect_equal(as.numeric(ppt_lppd(fit, scale = r, type = "plugin")), plugin, tolerance = 1e-8)
  }
})
