#ifndef POISTREE_SPATIAL_QUADRATURE_H
#define POISTREE_SPATIAL_QUADRATURE_H
#include <RcppArmadillo.h>

extern bool qpp_active;
extern arma::mat qpp_background;
extern arma::mat qpp_region;
extern arma::vec qpp_weights;
double qpp_box_exposure(const arma::mat& box);
#endif
