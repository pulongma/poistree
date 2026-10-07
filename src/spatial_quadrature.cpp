#include "spatial_quadrature.h"
// [[Rcpp::depends(RcppArmadillo)]]

bool qpp_active = false;
arma::mat qpp_background;
arma::mat qpp_region;
arma::vec qpp_weights;

// [[Rcpp::export]]
void qpp_set_quadrature(const arma::mat& background, const arma::vec& weights,
                        const arma::mat& region) {
  if (qpp_active) Rcpp::stop("Nested spatial-quadrature fits are unsupported.");
  if (background.n_rows == 0 || background.n_rows != weights.n_elem ||
      region.n_rows != background.n_cols || region.n_cols != 2 ||
      !background.is_finite() || !weights.is_finite() || !region.is_finite() ||
      arma::any(weights <= 0.0))
    Rcpp::stop("Supply finite background rows and one positive weight per row.");
  for (arma::uword j = 0; j < region.n_rows; ++j) {
    if (region(j, 1) <= region(j, 0) ||
        arma::any(background.col(j) < region(j, 0)) ||
        arma::any(background.col(j) > region(j, 1)))
      Rcpp::stop("Background must lie inside the covariate region.");
  }
  qpp_background = background;
  qpp_weights = weights;
  qpp_region = region;
  qpp_active = true;
}

// [[Rcpp::export]]
void qpp_clear_quadrature() {
  qpp_active = false;
  qpp_background.reset(); qpp_region.reset(); qpp_weights.reset();
}

double qpp_box_exposure(const arma::mat& box) {
  if (!qpp_active) return arma::prod(box.col(1) - box.col(0));
  double area = 0.0;
  for (arma::uword i = 0; i < qpp_weights.n_elem; ++i) {
    bool inside = true;
    for (arma::uword j = 0; j < box.n_rows; ++j) {
      // Use the same left-if-x<cut rule as hard-tree routing. Include the
      // outer upper boundary, which belongs to the final interval only.
      const bool last = box(j, 1) == qpp_region(j, 1);
      const double x = qpp_background(i, j);
      if (x < box(j, 0) || (last ? x > box(j, 1) : x >= box(j, 1))) {
        inside = false;
        break;
      }
    }
    if (inside) area += qpp_weights[i];
  }
  return area;
}

// [[Rcpp::export]]
double qpp_box_exposure_r(const arma::mat& box) {
  return qpp_box_exposure(box);
}
