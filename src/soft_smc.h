// Soft terminal-leaf PPT (S-PPT): particle structures and helpers shared by the
// SMC, conditional SMC and Particle Gibbs code in SMCtree.cpp.
//
// Model (recursive logistic gates, rates collapsed, colouring labels kept):
//   lambda(x) = sum_{leaves v} lambda_v phi_v(x),  lambda_v ~ Ga(a, b) (shape-rate)
//   phi_v = product of the gates on the root-to-v path,
//   right gate sigma(gate_j (x_j - c) / (b_j - a_j)), left gate = 1 - right.
// Each point i carries a label z_i = heap id of the leaf that generated it.
//   p(x, z | T, gate) = prod_v Q(m_v, H_v) prod_i phi_{z_i}(x_i),
//   Q(m, H) = b^a Gamma(a + m) / {Gamma(a) (b + H)^(a + m)}.
// Tree prior: S_h ~ Bern(rho_d), rho_d = rho (1 + depth)^-eta is the SPLIT
// probability (S = 1 split, as in PPT.cpp); axis uniform; cut uniform on a
// fixed per-axis grid (box free, required by ancestor sampling).
#ifndef POISTREE_SOFT_SMC_H
#define POISTREE_SOFT_SMC_H

#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <map>
#include <memory>
#include <utility>
#include <cmath>
#include <limits>
#include "soft_logistic.h"

struct SoftGateStep {
  int axis;
  double cut;
  int side;                       // -1 left, +1 right
};

struct SoftNodeP {
  bool active = false;
  int S = -1;                     // -1 undecided, 1 split, 0 leaf
  int J = -1;
  double L = NA_REAL;
  int depth = 0;
  int m = 0;                      // points coloured to this node
  double logH = 0.0;              // log exposure under the particle's gate
  std::vector<SoftGateStep> path;
};

// ---- shared-path store -------------------------------------------------------
// A coloured path node: heap id, geometric path (through `parent` and `gate`)
// and the set of points coloured to it.  Everything the forward step needs is
// a deterministic function of the node, so particles hold handles into the
// store and the expansion is computed once per distinct node.
struct SoftPathNode {
  int heap = 1, depth = 0, parent = -1;
  SoftGateStep gate{-1, 0.0, 0};  // edge from the parent (axis -1 at the root)
  std::vector<int> pts;
  std::vector<double> axisH;      // per-axis log integrals of the path
  double logH = 0.0;
  bool expanded = false;
  std::vector<double> axisL, axisR;      // per candidate: one-axis integral with the extra gate
  std::vector<double> logHL, logHR;      // per candidate: child log exposures
  std::vector<double> score, logprior, logq;   // index 0 = stop, then 1 + candidate
  double logQA = 0.0, logPhi = 0.0, log_inc = 0.0;
  std::unordered_map<int, std::shared_ptr<std::vector<std::vector<double> > > > pb;
  std::unordered_multimap<unsigned long long, std::pair<int, int> > children;
  double glued = 0.0;
  int glued_stamp = -1;
};

struct SoftRecord {               // a particle's view of one heap id
  int node;                       // store index
  signed char S;                  // -1 undecided, 1 split, 0 leaf
  int act;                        // candidate index drawn (split only)
};

struct SoftParticleS {
  std::map<int, SoftRecord> rec;  // heap id -> record
  double logw = 0.0;
};

// log left-gate values for every point and every grid cut, for one gate vector
struct SoftGateTable {
  int n_cand = 0;
  std::vector<int> offset, cand_axis;
  std::vector<double> cand_cut, loglft, lft;   // log and plain left-gate values
  double left(int i, int c) const { return loglft[(size_t)i * n_cand + c]; }
  double r(int i, int c) const { return lft[(size_t)i * n_cand + c]; }
  double right(int i, int c) const {
    double x = left(i, c);
    return x > -0.6931471805599453 ? std::log(-std::expm1(x)) : std::log1p(-std::exp(x));
  }
  int index(int axis, double cut) const {
    for (int c = offset[axis]; c < offset[axis + 1]; ++c) if (std::abs(cand_cut[c] - cut) < 1e-12) return c;
    return -1;
  }
};

// Reference trajectory for conditional SMC.  Inactive decisions and bits are
// prior draws given the active part; they are materialized on demand and the
// cache is cleared at every Gibbs iteration.
struct SoftRef {
  std::vector<SoftNodeP> node;
  std::vector<int> z;
  std::unordered_map<long long, char> bits;
  bool valid = false;
};

struct SoftModel {
  arma::mat X;
  arma::mat region;
  double a, b;                    // leaf Gamma(a, b)
  double rho, eta;                // split probability rho (1 + d)^-eta
  int Dmax;
  int n_nodes;                    // 2^(Dmax + 1)
  std::vector<std::vector<double> > grid;   // fixed cut grid per axis
  int exact_max;                  // exact Poisson-binomial proposal when m_A <= exact_max
  double defensive;               // mixture weight on the prior action proposal
  std::vector<double> lgam;       // lgamma(a + k), k = 0..n
  double a_log_b;                 // a log b

  void init_tables() {
    lgam.resize(n() + 1);
    for (int k = 0; k <= n(); ++k) lgam[k] = R::lgammafn(a + k);
    a_log_b = a * std::log(b);
  }
  // log Q(m, H) = log int lambda^m e^{-lambda H} Ga(lambda; a, b) d lambda
  double logQ(int m, double H) const { return a_log_b - lgam[0] + lgam[m] - (a + m) * std::log(b + H); }
  double logQ_lb(int m, double log_b_plus_H) const { return a_log_b - lgam[0] + lgam[m] - (a + m) * log_b_plus_H; }

  int d() const { return X.n_cols; }
  int n() const { return X.n_rows; }
  double width(int j) const { return region(j, 1) - region(j, 0); }
  double rho_depth(int depth) const {
    double r = rho * std::pow(1.0 + depth, -eta);
    return std::min(1.0 - 1e-12, std::max(1e-12, r));
  }
  int depth_of(int h) const { int dpt = 0; while (h > 1) { h /= 2; ++dpt; } return dpt; }
};

static inline double soft_log_Q(int m, double H, double a, double b) {
  return a * std::log(b) - R::lgammafn(a) + R::lgammafn(a + m) - (a + m) * std::log(b + H);
}

static inline double soft_logaddexp(double a, double b) {
  if (a == -std::numeric_limits<double>::infinity()) return b;
  if (b == -std::numeric_limits<double>::infinity()) return a;
  return std::max(a, b) + std::log1p(std::exp(-std::abs(a - b)));
}

// log of the left gate value at x for a split (axis, cut)
static inline double soft_log_left(const SoftModel& M, const arma::vec& gate,
                                   int axis, double cut, double xj) {
  double zz = gate[axis] * (xj - cut) / M.width(axis);
  return pst_logistic_log_right(-zz);
}

static inline double soft_log_gate(const SoftModel& M, const arma::vec& gate,
                                   const SoftGateStep& g, double xj) {
  double zz = gate[g.axis] * (xj - g.cut) / M.width(g.axis);
  return g.side < 0 ? pst_logistic_log_right(-zz) : pst_logistic_log_right(zz);
}

static inline double soft_log_phi(const SoftModel& M, const arma::vec& gate,
                                  const std::vector<SoftGateStep>& path, int i) {
  double out = 0.0;
  for (const SoftGateStep& g : path) out += soft_log_gate(M, gate, g, M.X(i, g.axis));
  return out;
}

// log of the one-dimensional path integral on axis j, optionally with one
// extra gate (extra_side = -1 left, +1 right, 0 none) appended to the path.
static inline double soft_axis_log_integral(const SoftModel& M, const arma::vec& gate,
                                            const std::vector<SoftGateStep>& path, int j,
                                            double extra_cut = 0.0, int extra_side = 0) {
  std::vector<double> cuts;
  std::vector<int> sides;
  for (const SoftGateStep& g : path) if (g.axis == j) { cuts.push_back(g.cut); sides.push_back(g.side); }
  if (extra_side != 0) { cuts.push_back(extra_cut); sides.push_back(extra_side); }
  double Hj = pst_logistic_path_axis_integral(cuts, sides, M.region(j, 0), M.region(j, 1), gate[j]);
  return std::log(std::max(Hj, 1e-300));
}

static inline double soft_log_exposure(const SoftModel& M, const arma::vec& gate,
                                       const std::vector<SoftGateStep>& path) {
  double logH = 0.0;
  for (int j = 0; j < M.d(); ++j) logH += soft_axis_log_integral(M, gate, path, j);
  return logH;
}

// log Poisson-binomial pmf of the sum of independent Bern(r) bits (log r given),
// rolling O(m) memory; log space keeps the far tails, which the collapsed
// leaf marginals can amplify by hundreds of nats.
static inline std::vector<double> soft_pb_logpmf(const std::vector<double>& logr) {
  const double ninf = -std::numeric_limits<double>::infinity();
  std::vector<double> p(logr.size() + 1, ninf);
  p[0] = 0.0;
  for (size_t s = 0; s < logr.size(); ++s) {
    const double l1 = logr[s], l0 = l1 > -0.6931471805599453 ? std::log(-std::expm1(l1)) : std::log1p(-std::exp(l1));
    for (size_t k = s + 1; k-- > 0;) {
      p[k + 1] = soft_logaddexp(p[k + 1], l1 + p[k]);
      p[k] += l0;
    }
  }
  return p;
}

// Psi(j, c) = E_B[Q(N, H_l) Q(m - N, H_r)] for m > exact_max, through the exact
// Beta-integral representation
//   Psi = b^{2a} Gamma(2a + m) / Gamma(a)^2 * c_l^a c_r^a c_r^m
//         * int_0^1 u^{a-1} (1-u)^{a-1} prod_i (alpha_i + beta_i u) du,
//   c_l = 1/(b + H_l), c_r = 1/(b + H_r), alpha_i = 1 - r_i, beta_i = r_i c_l/c_r - (1 - r_i),
// whose integrand is log-concave in u; the product is integrated by a Laplace
// approximation at its mode (quadratic expansion, exact Gaussian integral over
// [0, 1]) and the Beta factor is evaluated at the mode.  O(m) per candidate and
// with the correct tail behaviour, unlike a normal approximation of N.
static inline double soft_log_psi_laplace(const SoftModel& M, const SoftGateTable& G,
                                          const std::vector<int>& pts, int c,
                                          double lbL, double lbR) {
  const int m = pts.size();
  if (m == 0) return M.logQ_lb(0, lbL) + M.logQ_lb(0, lbR);
  const double rho = std::exp(lbR - lbL);               // c_l / c_r
  std::vector<double> alpha(m), beta(m);
  double mu = 0.0;
  for (int s = 0; s < m; ++s) {
    const double r = G.r(pts[s], c);
    alpha[s] = 1.0 - r; beta[s] = r * rho - (1.0 - r); mu += r;
  }
  auto slope = [&](double u, double& g1, double& g2) {   // h0'(u), h0''(u)
    g1 = 0.0; g2 = 0.0;
    for (int s = 0; s < m; ++s) { const double q = beta[s] / (alpha[s] + beta[s] * u); g1 += q; g2 -= q * q; }
  };
  double g1, g2, lo = 0.0, hi = 1.0, u = std::min(1.0 - 1e-9, std::max(1e-9, mu / m));
  slope(0.0, g1, g2);
  if (g1 <= 0.0) u = 0.0;
  else {
    slope(1.0, g1, g2);
    if (g1 >= 0.0) u = 1.0;
    else for (int it = 0; it < 30; ++it) {                // safeguarded Newton on the concave h0
      slope(u, g1, g2);
      if (g1 > 0.0) lo = u; else hi = u;
      double step = -g1 / g2, un = u + step;
      if (!(un > lo && un < hi)) un = 0.5 * (lo + hi);
      if (std::abs(un - u) < 1e-7) { u = un; break; }
      u = un;
    }
  }
  slope(u, g1, g2);
  double h0 = m * (-lbR);
  for (int s = 0; s < m; ++s) h0 += std::log(alpha[s] + beta[s] * u);
  // int_0^1 exp(g1 (v - u) + g2 (v - u)^2 / 2) dv, g2 < 0
  const double kappa = std::max(-g2, 1e-12), sd = 1.0 / std::sqrt(kappa), shift = g1 / kappa;
  const double zlo = (0.0 - u - shift) / sd, zhi = (1.0 - u - shift) / sd;
  // log{Phi(zhi) - Phi(zlo)} through the tail that keeps precision
  double logP;
  if (zlo > 0.0) { const double a1 = R::pnorm(zlo, 0.0, 1.0, 0, 1), a2 = R::pnorm(zhi, 0.0, 1.0, 0, 1); logP = a1 + std::log1p(-std::exp(a2 - a1)); }
  else if (zhi < 0.0) { const double a1 = R::pnorm(zhi, 0.0, 1.0, 1, 1), a2 = R::pnorm(zlo, 0.0, 1.0, 1, 1); logP = a1 + std::log1p(-std::exp(a2 - a1)); }
  else logP = std::log(R::pnorm(zhi, 0.0, 1.0, 1, 0) - R::pnorm(zlo, 0.0, 1.0, 1, 0));
  const double logZ = 0.5 * std::log(2.0 * M_PI / kappa) + 0.5 * g1 * g1 / kappa + logP;
  const double width = std::min(0.25, std::min(sd, 1.0 / std::max(std::abs(g1), 1e-300)));
  const double ut = std::min(1.0 - width, std::max(width, u));   // Beta factor evaluated here
  return 2.0 * M.a_log_b - 2.0 * M.lgam[0] + R::lgammafn(2.0 * M.a + m) - M.a * (lbL + lbR) +
    h0 + (M.a - 1.0) * (std::log(ut) + std::log1p(-ut)) + logZ;
}

static inline double soft_lse(const std::vector<double>& v) {
  double mx = -std::numeric_limits<double>::infinity();
  for (double x : v) if (x > mx) mx = x;
  if (!std::isfinite(mx)) return mx;
  double s = 0.0;
  for (double x : v) s += std::exp(x - mx);
  return mx + std::log(s);
}

static inline int soft_sample_log(const std::vector<double>& logp) {
  double norm = soft_lse(logp), u = R::unif_rand(), cum = 0.0;
  for (size_t k = 0; k < logp.size(); ++k) {
    cum += std::exp(logp[k] - norm);
    if (u <= cum) return (int)k;
  }
  return (int)logp.size() - 1;
}

// Poisson-binomial prefix table in log space: row s holds log P(sum of the
// first s bits = k).  Log space avoids the underflow of long near-deterministic
// colourings, whose replay for the reference particle needs finite values.
static inline std::vector<std::vector<double> > soft_pb_table(const std::vector<double>& logr) {
  const double ninf = -std::numeric_limits<double>::infinity();
  size_t m = logr.size();
  std::vector<std::vector<double> > tab(m + 1, std::vector<double>(m + 1, ninf));
  tab[0][0] = 0.0;
  for (size_t s = 1; s <= m; ++s) {
    const double l1 = logr[s - 1], l0 = l1 > -0.6931471805599453 ? std::log(-std::expm1(l1)) : std::log1p(-std::exp(l1));
    for (size_t k = 0; k <= s; ++k)
      tab[s][k] = soft_logaddexp(l0 + tab[s - 1][k], k > 0 ? l1 + tab[s - 1][k - 1] : ninf);
  }
  return tab;
}

// log P(B | sum B = k) for independent Bern(r) bits, drawing or replaying B from
// the last bit backwards with the log prefix table.
static inline double soft_cond_bernoulli(const std::vector<std::vector<double> >& tab,
                                         const std::vector<double>& logr, int k,
                                         std::vector<char>& B, bool replay) {
  int m = logr.size(); double logq = 0.0;
  for (int s = m; s >= 1; --s) {
    double logp1 = k > 0 ? std::min(0.0, logr[s - 1] + tab[s - 1][k - 1] - tab[s][k]) : -std::numeric_limits<double>::infinity();
    if (k == s) logp1 = 0.0;                       // all remaining bits must be 1
    if (!replay) B[s - 1] = std::log(R::unif_rand()) < logp1 ? 1 : 0;
    logq += B[s - 1] ? logp1 : (logp1 > -0.6931471805599453 ? std::log(-std::expm1(logp1)) : std::log1p(-std::exp(logp1)));
    k -= B[s - 1];
  }
  return logq;
}

#endif
