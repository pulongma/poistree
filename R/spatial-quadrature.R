#' Fit a point-process tree with a weighted integration rule
#'
#' Model intensity as a function of event covariates and approximate its
#' physical-domain integral using background covariates and exposure weights.
#' For intensity \eqn{\lambda}, background rows \eqn{u_k}, and weights
#' \eqn{w_k}, the fitted likelihood uses
#' \eqn{\int\lambda\approx\sum_k w_k\lambda(u_k)}.
#'
#' @param x Finite numeric matrix of event covariates with at least one row
#'   and one column. Rows are events and columns are covariates in the same
#'   order as `region` and `background`. All rows must lie inside `region`.
#' @param region Finite numeric matrix with one row per covariate and two
#'   columns giving lower and strictly larger upper bounds. These bounds
#'   define the tree's covariate domain and gate scaling; their product is
#'   not used as the physical-domain exposure in this wrapper.
#' @param background Finite numeric matrix of covariates at integration
#'   locations, with at least one row and `ncol(x)` columns in the same
#'   order as `x`. All rows must lie within `region`. Integration locations
#'   need not coincide with events, and repeated covariate rows are allowed.
#' @param weights Numeric vector of length `nrow(background)`, containing one
#'   finite, strictly positive physical-domain exposure per integration row.
#'   For example, weights may be cell areas or space-time volumes. They are
#'   used as supplied, without normalization; `sum(weights)` is the total
#'   exposure. Intensity is expressed per unit of this exposure.
#' @param gating Character scalar:
#'   \describe{
#'     \item{`"soft"`}{Probabilistic split gates; the default.}
#'     \item{`"hard"`}{Deterministic split gates.}
#'   }
#' @param sampler Character scalar determined by `gating` by default:
#'   \describe{
#'     \item{`"pcg"`}{Required for soft gates; partially collapsed Gibbs.}
#'     \item{`"smc"`}{Required for hard gates; shared-engine SMC.}
#'   }
#'   Other gating/sampler combinations are rejected by this wrapper.
#' @param a Positive finite scalar Gamma leaf-intensity prior shape,
#'   default `0.5`.
#' @param b Gamma leaf-intensity prior rate. Default `NULL` resolves to
#'   `a * sum(weights) / nrow(x)`. Require a positive finite scalar for
#'   soft PCG; hard SMC also accepts zero under its unnormalized prior
#'   convention. See \link{ppt_controls} for the distinction.
#' @param ... Named controls passed to [ppt_fit()] for the selected backend;
#'   see \link{ppt_controls} for defaults, ranges, and usage. For soft PCG these
#'   include `gate_family`, `gate_scale`, `gate_structure`, chain controls,
#'   `predict_at`, and `test`. Hard SMC accepts particle controls, but its
#'   `engine` must be `"shared"`. Prediction and test matrices contain
#'   covariates with the same column order and domain as `x`.
#' @param informed Logical scalar, default `FALSE`. For soft PCG, `TRUE`
#'   enables cached hard-surrogate tree proposals with the same weighted
#'   exposure rule. Hard SMC requires `FALSE`.
#' @param proposal_temperature Finite numeric scalar in `(0, 1]`, default
#'   `0.5`. With informed soft PCG, multiplies the surrogate log score before
#'   exponentiation; smaller values flatten the informed proposal weights.
#'   Unused unless `informed = TRUE`.
#' @param proposal_defensive Finite numeric scalar in `[0, 1)`, default
#'   `0.1`. With informed soft PCG, the fraction of the proposal mixture
#'   assigned to uniform candidate selection, with the remainder assigned
#'   to informed weights. Unused unless `informed = TRUE`.
#'
#' @details
#' Soft PCG supports root- and node-scaled logistic gates and node-scaled
#' compact gates. Informed proposals use hard-routed event counts and
#' weighted background exposures for their surrogate scores; the acceptance
#' ratio uses the soft weighted-integration target and forward/reverse
#' proposal probabilities. This wrapper changes the integration measure;
#' it does not change the definition of covariate split candidates.
#'
#' Fits run sequentially within one R process. Nested weighted-integration
#' fits are unsupported. The temporary integration rule is cleared after
#' fitting, including on errors. Post-fit functions use the following rules:
#' \itemize{
#'   \item `fit$posterior$integrated_intensity_draws` and
#'     `fit$posterior$mean_integrated_intensity` store the weighted integrals;
#'     [ppt_lppd()] uses these integrals when scoring a test pattern.
#'   \item [ppt_lambda()] evaluates saved tree states at new covariates.
#'   \item [ppt_marginal()] integrates over covariate-box coordinates;
#'     it is not a physical-domain marginal for these fits.
#' }
#' @return A `ppt` object with the usual fit components described in
#'   [ppt_fit()]. In addition, `data$quadrature` contains `background`,
#'   `weights`, and `total_exposure = sum(weights)`; `model$integration`
#'   identifies weighted physical-domain integration. Resolved controls
#'   are stored in `control`.
#' @seealso [ppt_fit()], \link{ppt_controls}
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
#' Compare the exposure of a hard child box with the weighted exposure of
#' the corresponding single logistic gate over the full integration rule.
#'
#' @param background Finite numeric matrix of covariates at integration
#'   locations, with at least one row. Columns correspond to rows of `region`.
#'   Every row must lie inside `region`.
#' @param weights Numeric vector of length `nrow(background)` containing
#'   strictly positive, finite exposure weights. Values are not normalized.
#' @param region Finite numeric matrix with `ncol(background)` rows and two
#'   columns containing lower and strictly larger upper bounds.
#' @param axis Integer scalar from `1` to `ncol(background)` identifying the
#'   split coordinate, using R's one-based indexing.
#' @param cut Finite numeric scalar split location, strictly inside the
#'   interval `region[axis, ]`.
#' @param side Integer scalar selecting the child:
#'   \describe{
#'     \item{`-1L`}{Left child; the default. Hard routing uses `x[axis] < cut`.}
#'     \item{`1L`}{Right child. Hard routing uses `x[axis] >= cut`.}
#'   }
#' @param gate Positive finite logistic sharpness, default `10`. Supply a
#'   scalar or a vector of length `ncol(background)`; the selected coordinate
#'   uses `gate[axis]` after scalar recycling. The right-gate probability is
#'   `plogis(gate[axis] * (x[axis] - cut) / diff(region[axis, ]))`.
#'   This helper uses root-relative logistic scaling.
#' @return A list containing:
#'   \describe{
#'     \item{`hard_exposure`}{Sum of weights for background rows routed to
#'       the selected hard child. The outer upper region boundary is included.}
#'     \item{`soft`}{A list with `phi` (selected-child membership probability
#'       at each background row), `log_phi` (its log probability), and
#'       `H = sum(weights * phi)` (the soft child exposure).}
#'   }
#' @seealso [ppt_fit_quadrature()]
#' @md
#' @keywords internal
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
