// Rcpp entry point for the soft terminal-leaf PPT Particle Gibbs sampler
// (conditional SMC with ancestor sampling); see SMCtree.cpp and soft_smc.h.
#ifndef _USE_Armadillo
#define _USE_Armadillo
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#endif

#ifndef _USE_ProgressBar
#define _USE_ProgressBar
#include <R_ext/Utils.h>
#include <iostream>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#endif

using namespace Rcpp;

#include "PPT.h"
#include "SMCtree.h"

// [[Rcpp::export]]
Rcpp::List SPPT_fit_PGAS(const arma::mat& X, const arma::mat& grid, const arma::mat& Xtest,
                         const arma::mat& region, Rcpp::List cut_grid,
                         double a, double b, arma::vec gate,
                         arma::vec a_gate, arma::vec b_gate, arma::vec sd_gate, arma::vec gate_min,
                         bool gate_shared, double rho, double eta, int max_depth,
                         int P, int niter, int burn, int thin, int label_sweeps,
                         bool update_gate, bool ancestor_sampling, int exact_max,
                         double defensive, bool resample_node, double ess_threshold,
                         bool allocation_rates, bool verbose) {
  const int d = X.n_cols;
  if (X.n_rows == 0 || d == 0) stop("X must be a non-empty matrix");
  if ((int)region.n_rows != d || region.n_cols != 2) stop("region must be a d by 2 matrix");
  if (cut_grid.size() != d) stop("cut_grid must have one vector per axis");
  if (a <= 0.0 || b <= 0.0) stop("a and b must be positive for the soft PPT");
  if (rho <= 0.0 || rho >= 1.0 || eta < 0.0) stop("require 0 < rho < 1 and eta >= 0");
  if (max_depth < 1 || P < 2 || niter <= burn || burn < 0 || thin < 1 || label_sweeps < 0)
    stop("invalid Particle-Gibbs controls");
  if (defensive < 0.0 || defensive >= 1.0) stop("defensive must lie in [0, 1)");
  if (!(ess_threshold > 0.0 && ess_threshold <= 1.0)) stop("ess_threshold must lie in (0, 1]");
  if ((int)gate.n_elem != d || (int)a_gate.n_elem != d || (int)b_gate.n_elem != d ||
      (int)sd_gate.n_elem != d || (int)gate_min.n_elem != d)
    stop("gate parameters must have length d");

  SoftModel M;
  M.X = X; M.region = region; M.a = a; M.b = b; M.rho = rho; M.eta = eta;
  M.Dmax = max_depth; M.n_nodes = 1 << (max_depth + 1);
  M.exact_max = exact_max; M.defensive = defensive;
  M.init_tables();
  M.grid.resize(d);
  for (int j = 0; j < d; ++j) {
    NumericVector g = cut_grid[j];
    if (g.size() == 0) stop("every axis needs at least one grid cut");
    for (double c : g) {
      if (!(c > region(j, 0) && c < region(j, 1))) stop("grid cuts must lie strictly inside the region");
      M.grid[j].push_back(c);
    }
  }

  SoftSMCtree smc(M, P, ancestor_sampling);
  smc.resample_node = resample_node;
  smc.ess_threshold = ess_threshold;
  smc.alloc_rates = allocation_rates;
  return smc.PGAS(grid, Xtest, gate, a_gate, b_gate, sd_gate, gate_min, gate_shared,
                  niter, burn, thin, label_sweeps, update_gate, verbose);
}
