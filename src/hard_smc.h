// Shared-path SMC for the hard terminal-leaf Poisson process tree (PPT).
//
// A particle is a sparse map from heap position to (store node, decision).
// Store nodes are boxes with their observation index sets; because the
// box determines the observation set, two particles that reach the same box
// share one node, its candidate cuts and their scores (shared_path_smc_hard_ppt.tex).
#ifndef POISTREE_HARD_SMC_H
#define POISTREE_HARD_SMC_H

#include <RcppArmadillo.h>
#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <limits>

struct HardNode {
  int heap = 1, depth = 0, parent = -1, cand = -1;   // cand: candidate index at the parent
  arma::mat box;                                     // d x 2
  std::vector<int> pts;                              // observation indices (released after expansion)
  int m = 0;
  double area = 1.0;
  bool expanded = false;
  std::vector<int> cand_axis, cand_nL;               // valid candidate splits
  std::vector<double> cand_cut, score;               // score[k]: log prior x marginal likelihood of candidate k
  std::vector<double> logq;                          // index 0 = stop, then 1 + candidate
  double logQ0 = 0.0, logPhi = 0.0, log_inc = 0.0;
  std::unordered_map<int, std::pair<int, int> > children;   // candidate -> (left, right)
};

struct HardRecord { int node; signed char S; int act; };   // S: -1 undecided, 1 split, 0 leaf

struct HardParticle {
  std::map<int, HardRecord> rec;
  double logw = 0.0;
};

class HardSMCtree {
public:
  const arma::mat& X;
  arma::mat region;
  int d, n, Dmax, min_leaf, cut_grid_n, P;
  double a, b, rho, eta, max_aspect, resample_thresh;
  std::vector<HardNode> store;
  std::vector<HardParticle> particles;
  arma::vec weights, ESS_hist, logZ_inc, logZ_run;
  double logZ = 0.0;
  int n_expanded = 0, n_resampled = 0;

  HardSMCtree(const arma::mat& X_, const arma::mat& region_, int Dmax_, int min_leaf_, int cut_grid_n_,
              int P_, double a_, double b_, double rho_, double eta_, double max_aspect_, double resample_thresh_)
    : X(X_), region(region_), d(X_.n_cols), n(X_.n_rows), Dmax(Dmax_), min_leaf(min_leaf_),
      cut_grid_n(cut_grid_n_), P(P_), a(a_), b(b_), rho(rho_), eta(eta_), max_aspect(max_aspect_),
      resample_thresh(resample_thresh_) {}

  // Gamma-Poisson marginal of a leaf with m observations and volume `area`;
  // b = 0 keeps the package's improper-prior limit (PPT::PPT_base_mloglik).
  double logQ0(int m, double area) const {
    if (b > 0.0) return std::lgamma(m + a) - std::lgamma(a) + a * std::log(b) - (m + a) * std::log(b + area);
    return std::lgamma(m + a) - (m + a) * std::log(area);
  }
  double rho_depth(int depth) const {
    return std::min(1.0 - 1e-12, std::max(1e-12, rho * std::pow(1.0 + depth, -eta)));
  }

  int make_root();
  int make_child(int parent, int cand, int side, std::vector<int>&& pts);
  void expand_node(int v);
  bool sample_position(int t);
  void resample();
  void sweep();
  Rcpp::List export_particles(const arma::mat& grid, arma::mat& lam_draws, arma::vec& loglik,
                              arma::vec& lppd, arma::vec& integral);
};

#endif
