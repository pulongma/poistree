#ifndef POISTREE_SPATIAL_QUADRATURE_H
#define POISTREE_SPATIAL_QUADRATURE_H
#include <RcppArmadillo.h>

// Quadrature fitting runs synchronously within one R process. The R
// wrapper installs one integration rule for the duration of a fit and clears
// it with on.exit(), including on errors. Separate R processes have independent
// rules. No original sampler proposal or acceptance step is changed.
extern bool qpp_active;
extern arma::mat qpp_background;
extern arma::mat qpp_region;
extern arma::vec qpp_weights;
double qpp_box_exposure(const arma::mat& box);
#endif
