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
                         bool gate_shared, double rho, double eta, double max_depth,
                         int P, int niter, int burn, int thin, int label_sweeps,
                         bool update_gate, bool ancestor_sampling, int exact_max,
                         double defensive, bool resample_node, double ess_threshold,
                         bool allocation_rates, bool verbose) {
  const int depth = ppt_checked_depth(max_depth, 1);
  ppt_check_soft_storage(depth);
  const int d = X.n_cols;
  if (X.n_rows == 0 || d == 0) stop("X must be a non-empty matrix");
  if ((int)region.n_rows != d || region.n_cols != 2) stop("region must be a d by 2 matrix");
  if (cut_grid.size() != d) stop("cut_grid must have one vector per axis");
  if (a <= 0.0 || b <= 0.0) stop("a and b must be positive for the soft PPT");
  if (rho <= 0.0 || rho >= 1.0 || eta < 0.0) stop("require 0 < rho < 1 and eta >= 0");
  if (depth < 1 || P < 2 || niter <= burn || burn < 0 || thin < 1 || label_sweeps < 0)
    stop("invalid Particle-Gibbs controls");
  if (defensive < 0.0 || defensive >= 1.0) stop("defensive must lie in [0, 1)");
  if (!(ess_threshold > 0.0 && ess_threshold <= 1.0)) stop("ess_threshold must lie in (0, 1]");
  if ((int)gate.n_elem != d || (int)a_gate.n_elem != d || (int)b_gate.n_elem != d ||
      (int)sd_gate.n_elem != d || (int)gate_min.n_elem != d)
    stop("gate parameters must have length d");

  SoftModel M;
  M.X = X; M.region = region; M.a = a; M.b = b; M.rho = rho; M.eta = eta;
  M.Dmax = depth; M.n_nodes = ppt_tree_slots(depth) + 1;
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

// Internal diagnostic for the exact cut-tilt calculation.  The direct branch
// repeats the Poisson-binomial recursion independently at every cut so tests
// can compare the optimized scores and all count probabilities.
// [[Rcpp::export]]
Rcpp::List SPPT_exact_cut_scores(Rcpp::NumericVector x,
                                Rcpp::NumericVector cuts,
                                double gate, double width,
                                Rcpp::NumericVector H_left,
                                Rcpp::NumericVector H_right,
                                double a, double b) {
  const int m = x.size(), ncut = cuts.size();
  if (ncut == 0 || H_left.size() != ncut || H_right.size() != ncut)
    stop("cuts and child exposures must have the same positive length");
  if (!std::isfinite(gate) || gate < 0.0 || !std::isfinite(width) || width <= 0.0 ||
      !std::isfinite(a) || a <= 0.0 || !std::isfinite(b) || b <= 0.0)
    stop("require finite gate >= 0 and width, a, b > 0");
  for (double xi : x) if (!std::isfinite(xi)) stop("x must be finite");
  for (int c = 0; c < ncut; ++c) {
    if (!std::isfinite(cuts[c]) || !std::isfinite(H_left[c]) ||
        !std::isfinite(H_right[c]) || H_left[c] < 0.0 || H_right[c] < 0.0)
      stop("cuts must be finite and child exposures finite and nonnegative");
  }

  SoftModel M;
  M.X.set_size(m, 1);
  M.a = a; M.b = b;
  M.init_tables();
  const double lower = *std::min_element(cuts.begin(), cuts.end());
  const double upper = *std::max_element(cuts.begin(), cuts.end());
  const double reference = lower + 0.5 * (upper - lower);
  std::vector<double> logits(m);
  for (int i = 0; i < m; ++i) logits[i] = gate * (reference - x[i]) / width;
  const std::shared_ptr<SoftExactAxis> axis = soft_exact_axis(logits, reference, gate, width, M);
  const std::vector<double>& reference_pmf = axis->logpmf;
  NumericVector log_scores(ncut), direct_log_scores(ncut), original_log_scores(ncut);
  NumericMatrix log_pmf(ncut, m + 1), direct_log_pmf(ncut, m + 1);
  NumericMatrix log_allocation_count(ncut, m + 1);

  for (int c = 0; c < ncut; ++c) {
    const double delta = gate * (cuts[c] - reference) / width;
    const double lbL = std::log(b + H_left[c]), lbR = std::log(b + H_right[c]);
    double count_norm;
    log_scores[c] = soft_exact_log_psi(*axis, delta, M, lbL, lbR, &count_norm);
    original_log_scores[c] = soft_tilt_log_psi(reference_pmf, delta, M, lbL, lbR);
    std::vector<double> tilted(m + 1);
    for (int k = 0; k <= m; ++k) tilted[k] = reference_pmf[k] + k * delta;
    const double log_norm = soft_lse(tilted);
    for (int k = 0; k <= m; ++k) {
      log_pmf(c, k) = tilted[k] - log_norm;
      log_allocation_count(c, k) = axis->coef[k] + k * (delta + lbR - lbL) - count_norm;
    }

    for (int i = 0; i < m; ++i) logits[i] = gate * (cuts[c] - x[i]) / width;
    const std::vector<double> direct = soft_pb_logpmf_logodds(logits);
    std::vector<double> terms(m + 1);
    for (int k = 0; k <= m; ++k) {
      direct_log_pmf(c, k) = direct[k];
      terms[k] = direct[k] + M.logQ_lb(k, lbL) + M.logQ_lb(m - k, lbR);
    }
    direct_log_scores[c] = soft_lse(terms);
  }
  return List::create(_["log_scores"] = log_scores,
                      _["original_log_scores"] = original_log_scores,
                      _["direct_log_scores"] = direct_log_scores,
                      _["log_pmf"] = log_pmf,
                      _["direct_log_pmf"] = direct_log_pmf,
                      _["log_allocation_count"] = log_allocation_count);
}

// Internal diagnostic: exercise the production coefficient count sampler and
// lazy coordinate-level conditional-allocation table across multiple cuts.
// [[Rcpp::export]]
Rcpp::List SPPT_exact_reuse_draws(Rcpp::NumericVector logits,
                                 Rcpp::NumericVector shifts,
                                 Rcpp::NumericVector H_left,
                                 Rcpp::NumericVector H_right,
                                 double a, double b, int draws,
                                 int fixed_count = -1) {
  const int m = logits.size(), ncut = shifts.size();
  if (ncut == 0 || H_left.size() != ncut || H_right.size() != ncut || draws < 1 ||
      fixed_count < -1 || fixed_count > m)
    stop("invalid dimensions, draw count, or fixed allocation count");
  if (!(a > 0.0) || !(b > 0.0)) stop("a and b must be positive");
  std::vector<double> odds(logits.begin(), logits.end());
  for (double x : odds) if (!std::isfinite(x)) stop("logits must be finite");
  SoftModel M;
  M.X.set_size(m, 1); M.a = a; M.b = b; M.init_tables();
  std::shared_ptr<SoftExactAxis> axis = soft_exact_axis(odds, 0.0, 1.0, 1.0, M);
  List allocations(ncut);
  IntegerVector prefix_builds(ncut);
  NumericVector log_scores(ncut);
  NumericMatrix log_count(ncut, m + 1);
  for (int c = 0; c < ncut; ++c) {
    if (!std::isfinite(shifts[c]) || !std::isfinite(H_left[c]) || H_left[c] < 0.0 ||
        !std::isfinite(H_right[c]) || H_right[c] < 0.0)
      stop("shifts must be finite and exposures finite and nonnegative");
    const double lbL = std::log(b + H_left[c]), lbR = std::log(b + H_right[c]);
    double norm;
    log_scores[c] = soft_exact_log_psi(*axis, shifts[c], M, lbL, lbR, &norm);
    for (int k = 0; k <= m; ++k)
      log_count(c, k) = axis->coef[k] + k * (shifts[c] + lbR - lbL) - norm;
    IntegerMatrix out(draws, m);
    std::vector<char> bits(m);
    for (int r = 0; r < draws; ++r) {
      const int k = fixed_count < 0 ?
        soft_exact_sample_count(*axis, shifts[c], lbL, lbR, norm) : fixed_count;
      if (soft_exact_draw_bits(*axis, k, bits)) ++prefix_builds[c];
      for (int i = 0; i < m; ++i) out(r, i) = bits[i];
    }
    allocations[c] = out;
  }
  int prefix_cells = 0;
  for (const auto& row : axis->prefix) prefix_cells += row.size();
  return List::create(_["allocations"] = allocations, _["log_count"] = log_count,
                      _["log_scores"] = log_scores, _["prefix_builds"] = prefix_builds,
                      _["prefix_cells"] = prefix_cells);
}

// Internal diagnostic for replay of a fixed allocation, including allocations
// whose nonpreferred routing probability underflows on the probability scale.
// [[Rcpp::export]]
Rcpp::List SPPT_exact_allocation_logprob(Rcpp::NumericVector logits,
                                        Rcpp::IntegerVector bits) {
  const int m = logits.size();
  if (bits.size() != m) stop("logits and bits must have the same length");
  std::vector<double> log_left(m), log_right(m), odds(m);
  std::vector<char> allocation(m);
  arma::mat gate_points(m, 1);
  SoftGateTable G;
  G.n_cand = 1;
  G.cand_axis = {0}; G.cand_cut = {0.0};
  G.points = &gate_points;
  G.gate = arma::vec(1, arma::fill::ones);
  G.width = arma::vec(1, arma::fill::ones);
  G.loglft.resize(m);
  NumericVector gate_log_right(m);
  int k = 0;
  double log_joint = 0.0;
  for (int i = 0; i < m; ++i) {
    if (!std::isfinite(logits[i]) || (bits[i] != 0 && bits[i] != 1))
      stop("logits must be finite and bits must be zero or one");
    odds[i] = logits[i];
    log_left[i] = pst_logistic_log_right(logits[i]);
    log_right[i] = pst_logistic_log_right(-logits[i]);
    gate_points(i, 0) = -logits[i];
    G.loglft[i] = log_left[i];
    gate_log_right[i] = G.right(i, 0);
    allocation[i] = bits[i];
    k += bits[i];
    log_joint += bits[i] ? log_left[i] : log_right[i];
  }
  const std::vector<std::vector<double> > table = soft_pb_table_logprobs(log_left, log_right);
  const std::vector<double> pmf = soft_pb_logpmf_logodds(odds);
  const double replay = soft_cond_bernoulli_logprobs(table, log_left, log_right, k, allocation, true);
  return List::create(_["log_conditional"] = replay,
                      _["direct_log_conditional"] = log_joint - pmf[k],
                      _["log_count"] = pmf[k],
                      _["log_joint"] = log_joint,
                      _["gate_log_right"] = gate_log_right);
}
