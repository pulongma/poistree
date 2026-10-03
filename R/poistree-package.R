#' poistree: Bayesian Poisson point-process tree models
#'
#' The `poistree` package estimates Poisson point-process intensity functions
#' with adaptive axis-aligned trees. The unified [ppt_fit()] interface
#' separates modeling and computational choices:
#'
#' \itemize{
#'   \item `gating`: hard partitions or soft logistic/compact gates;
#'   \item `scales`: terminal-leaf or additive multiscale intensities;
#'   \item `sampler`: sequential Monte Carlo, reversible-jump MCMC, informed
#'     reversible-jump MCMC, or Particle Gibbs with exact conditional SMC. The
#'     legacy token `pgas` selects the Particle-Gibbs backend; ancestor
#'     sampling is currently disabled.
#' }
#'
#' The fitted object has class `ppt`. Use [ppt_predict()] for posterior
#' intensity summaries, [ppt_marginal()] for one-input marginal intensity
#' curves, [plot.ppt()] for fitted intensities, [ppt_summary()]
#' for model summaries, [ppt_logLik()] and [ppt_lppd()] for likelihood and
#' predictive evaluation, and [ppt_diagnostics()] for sampler diagnostics.
#'
#' @section Model configurations:
#' \tabular{llll}{
#' Model \tab `gating` \tab `scales` / prior \tab `sampler` \cr
#' PPT \tab `hard` \tab `leaf` \tab `smc` \cr
#' PPT \tab `hard` \tab `leaf` \tab `rjmcmc` \cr
#' PPT \tab `hard` \tab `leaf` \tab `pgas` \cr
#' S-PPT \tab `soft` \tab `leaf` \tab `rjmcmc` \cr
#' MPPT \tab `hard` \tab `multiscale`, Markov or independent \tab `rjmcmc`; `irjmcmc` for independent scales \cr
#' S-MPPT \tab `soft` \tab `multiscale`, Markov or independent \tab `rjmcmc`; `irjmcmc` for independent scales \cr
#' }
#'
#' The default is S-MPPT: soft gating, independent multiscale scales, and
#' RJ-MCMC. Setting only `scales = "leaf"` selects S-PPT. Explicit
#' `gating = "hard"` calls retain MPPT and PPT.
#'
#' @section Basic workflow:
#' Supply the observed event locations or covariate vectors as an `n` by `d`
#' numeric matrix and the bounded observation window as a `d` by 2 matrix.
#' The optional `predict_at` argument precomputes posterior intensities at
#' selected locations. It is not required: arbitrary locations can be
#' evaluated after fitting with [ppt_lambda()] or [ppt_predict()].
#'
#' @examples
#' \dontrun{
#' set.seed(2026)
#' x <- matrix(runif(400), ncol = 2)
#' region <- matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE)
#' grid <- as.matrix(expand.grid(
#'   x1 = seq(0, 1, length.out = 30),
#'   x2 = seq(0, 1, length.out = 30)
#' ))
#'
#' # Hard terminal-leaf PPT fitted by SMC
#' fit_ppt <- ppt_fit(
#'   x, region, gating = "hard", scales = "leaf", sampler = "smc",
#'   predict_at = grid, particles = 500
#' )
#'
#' # Hard terminal-leaf PPT fitted by RJ-MCMC
#' fit_ppt_mcmc <- ppt_fit(
#'   x, region, gating = "hard", scales = "leaf", sampler = "rjmcmc",
#'   predict_at = grid, chains = 2, iter = 2000, burn = 500
#' )
#'
#' ppt_summary(fit_ppt_mcmc)
#' intensity <- ppt_predict(fit_ppt_mcmc, type = "interval")
#' plot(fit_ppt_mcmc)
#' ppt_diagnostics(fit_ppt_mcmc)
#'
#' # Hard terminal-leaf PPT fitted by conditional-SMC Particle Gibbs
#' fit_ppt_pgas <- ppt_fit(
#'   x, region, gating = "hard", scales = "leaf", sampler = "pgas",
#'   predict_at = grid, particles = 200, iter = 500, burn = 100
#' )
#'
#' # Soft multiscale PPT; this is also the default model. The remaining
#' # defaults select independent scales and RJ-MCMC.
#' fit_smppt <- ppt_fit(
#'   x, region, predict_at = grid,
#'   chains = 2, iter = 2000, burn = 500
#' )
#'
#' # Soft terminal-leaf PPT; soft gating and RJ-MCMC remain at their defaults.
#' fit_sppt <- ppt_fit(
#'   x, region, scales = "leaf", predict_at = grid,
#'   chains = 2, iter = 2000, burn = 500
#' )
#' }
#'
#' @docType package
#' @name poistree
#' @aliases poistree-package
#' @keywords internal
"_PACKAGE"
