#' poistree: Bayesian Poisson point-process tree models
#'
#' The `poistree` package estimates Poisson point-process intensity functions
#' with adaptive axis-aligned trees. The unified [ppt_fit()] interface
#' separates modeling and computational choices:
#'
#' \itemize{
#'   \item `gating`: the partition mechanism; hard partitions are currently
#'     implemented and other gates are reserved for future backends;
#'   \item `scales`: the intensity representation; terminal-leaf intensities
#'     are currently implemented and multiscale intensities are reserved for
#'     future backends;
#'   \item `sampler`: sequential Monte Carlo, reversible-jump MCMC, or Particle
#'     Gibbs with ancestor sampling.
#' }
#'
#' The fitted object has class `ppt`. Use [ppt_predict()] for posterior
#' intensity summaries, [plot.ppt()] for fitted intensities, [ppt_summary()]
#' for model summaries, [ppt_logLik()] and [ppt_lppd()] for likelihood and
#' predictive evaluation, and [ppt_diagnostics()] for sampler diagnostics.
#'
#' @section Model configurations:
#' \tabular{llll}{
#' Model \tab `gating` \tab `scales` / prior \tab `sampler` \cr
#' PPT \tab `hard` \tab `leaf` \tab `smc` \cr
#' PPT \tab `hard` \tab `leaf` \tab `rjmcmc` \cr
#' PPT \tab `hard` \tab `leaf` \tab `pgas` \cr
#' }
#'
#' The component-based wrapper is intentionally retained. Future backends can
#' be registered internally without changing calls to [ppt_fit()].
#'
#' @section Basic workflow:
#' Supply the observed event locations or covariate vectors as an `n` by `d`
#' numeric matrix and the bounded observation window as a `d` by 2 matrix.
#' Prediction locations are passed through `predict_at`.
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
#' # Hard terminal-leaf PPT fitted by PGAS
#' fit_ppt_pgas <- ppt_fit(
#'   x, region, gating = "hard", scales = "leaf", sampler = "pgas",
#'   predict_at = grid, particles = 200, iter = 500, burn = 100
#' )
#' }
#'
#' @docType package
#' @name poistree
#' @aliases poistree-package
#' @keywords internal
"_PACKAGE"
