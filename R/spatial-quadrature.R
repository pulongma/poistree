#' Fit a point-process tree with a weighted integration rule
#'
#' Model intensity as a function of covariates while integrating over a physical
#' domain using background covariates and their area or space-time exposures.
#'
#' @param x Matrix of event covariates, one event per row.
#' @param region Matrix with one row per covariate and lower/upper bounds in
#'   its two columns. These bounds define tree splits, not physical exposure.
#' @param background Matrix of covariates at integration locations.
#' @param weights One finite positive physical-domain exposure per background
#'   row. The integrated intensity is the weighted sum over these rows.
#' @param gating Either `"soft"` (default) or `"hard"`.
#' @param sampler `"pcg"` for soft gates or `"smc"` for hard gates.
#' @param a,b Gamma leaf-intensity prior shape and rate. The default rate is
#'   `a * sum(weights) / nrow(x)`.
#' @param ... Other controls passed to [ppt_fit()], including `predict_at`,
#'   `test`, `gate_scale`, chain settings, `seed`, and `verbose`.
#' @param informed Logical; use hard-surrogate informed tree proposals for
#'   soft PCG. The default `FALSE` preserves standard PCG proposals.
#' @param proposal_temperature Informed PCG score temperature in `(0, 1]`,
#'   default `0.5`. Used only when `informed = TRUE`.
#' @param proposal_defensive Uniform proposal mixture weight in `[0, 1)`,
#'   default `0.1`. Used only when `informed = TRUE`.
#'
#' @details Soft PCG supports both root- and node-scaled gates. Informed
#'   proposals use hard-routed counts and weighted background exposures for
#'   proposal scores; acceptance uses the soft quadrature target and the
#'   forward/reverse proposal probabilities. Hard SMC uses the shared engine.
#'   Other samplers are not supported by this wrapper.
#'
#'   Fits run sequentially within one R process. The temporary integration
#'   rule is cleared after fitting, including on errors. [ppt_integral()] and
#'   [ppt_lppd()] use the stored weighted integrals; [ppt_lambda()] evaluates
#'   saved tree states at new covariates. `ppt_marginal()` still integrates
#'   over covariate-box coordinates and is not a physical-domain marginal
#'   for these fits.
#' @return A `ppt` fit with the integration rule in `data$quadrature` and
#'   the selected proposal controls in `control`.
#' @md
#' @export
#' @examples
#' set.seed(7)
#' x <- matrix(runif(40), ncol = 2)
#' background <- as.matrix(expand.grid(seq(0, 1, length.out = 5),
#'                                    seq(0, 1, length.out = 5)))
#' region <- rbind(c(0, 1), c(0, 1))
#' fit <- ppt_fit_quadrature(
#'   x, region, background, rep(1 / nrow(background), nrow(background)),
#'   sampler = "pcg", informed = TRUE, gate_scale = "node",
#'   chains = 1, iter = 30, burn = 10, thin = 2,
#'   max_depth = 2, cut_candidates = 3, verbose = FALSE
#' )
#' fit$control$informed
ppt_fit_quadrature <- function(x, region, background, weights,
                               gating = c("soft", "hard"),
                               sampler = if (gating == "soft") "pcg" else "smc",
                               a = 0.5, b = NULL, ...,
                               informed = FALSE, proposal_temperature = 0.5,
                               proposal_defensive = 0.1) {
  gating <- match.arg(gating)
  if ((gating == "soft" && sampler != "pcg") ||
      (gating == "hard" && sampler != "smc")) {
    stop("Quadrature fitting supports soft PCG or hard SMC only.", call. = FALSE)
  }
  if (!is.logical(informed) || length(informed) != 1L || is.na(informed)) {
    stop("`informed` must be a logical scalar.", call. = FALSE)
  }
  if (gating == "hard" && informed) {
    stop("`informed = TRUE` requires soft PCG quadrature fitting.", call. = FALSE)
  }
  x <- .ppt_validate_points(x, name = "x")
  region <- .ppt_validate_region(region, ncol(x))
  background <- .ppt_validate_points(background, ncol(x), region, "background")
  if (!is.numeric(weights) || length(weights) != nrow(background) ||
      any(!is.finite(weights)) || any(weights <= 0)) {
    stop("Supply one finite, positive physical-domain weight per background row.",
         call. = FALSE)
  }
  if (is.null(b)) b <- a * sum(weights) / nrow(x)
  dots <- list(...)
  if (gating == "hard") {
    if (!is.null(dots$engine) && dots$engine != "shared") {
      stop("Quadrature fitting requires the hard shared SMC engine.", call. = FALSE)
    }
    dots$engine <- "shared"
  } else {
    dots$informed <- informed
    dots$proposal_temperature <- proposal_temperature
    dots$proposal_defensive <- proposal_defensive
  }
  qpp_set_quadrature(background, as.numeric(weights), region)
  on.exit(qpp_clear_quadrature(), add = TRUE)
  fit <- do.call(ppt_fit, c(list(x = x, region = region, gating = gating,
                                sampler = sampler, a = a, b = b), dots))
  fit$data$quadrature <- list(background = background, weights = weights,
                             total_exposure = sum(weights))
  fit$model$integration <- "weighted physical-domain quadrature"
  fit$model$adapter_version <- "1"
  fit$call <- match.call()
  fit
}

#' Inspect a one-split quadrature geometry
#'
#' @param background,weights,region Integration rule and covariate bounds, as
#'   in [ppt_fit_quadrature()].
#' @param axis One-based split coordinate.
#' @param cut Split location.
#' @param side `-1L` for the left child, `1L` for the right child.
#' @param gate Logistic gate sharpness.
#' @return Hard weighted exposure and soft leaf-geometry diagnostics.
#' @md
#' @export
ppt_quadrature_geometry <- function(background, weights, region, axis, cut,
                                    side = -1L, gate = 10) {
  qpp_set_quadrature(as.matrix(background), as.numeric(weights), as.matrix(region))
  on.exit(qpp_clear_quadrature(), add = TRUE)
  box <- region
  box[axis, if (side == -1L) 2L else 1L] <- cut
  soft <- ppstree_geometry(as.integer(axis), cut, diff(region[axis, ]),
                           as.integer(side), as.matrix(background), region, gate)
  list(hard_exposure = qpp_box_exposure_r(box), soft = soft)
}
