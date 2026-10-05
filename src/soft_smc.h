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
#include "tree_limits.h"
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

// Exact count and allocation data for one node/coordinate.  The reference
// distribution changes with the node's observations and the coordinate gate,
// but is shared by every cut on that coordinate.  The prefix table remains
// empty until an interior allocation count is actually drawn.
struct SoftExactAxis {
  double refcut = 0.0, gate = 0.0, width = 1.0;
  bool finite = true;
  std::vector<double> logits, log_left, log_right, logpmf, coef;
  std::vector<std::vector<double> > prefix;
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
  bool hard_scored = false;
  std::vector<double> axisL, axisR;      // per candidate: one-axis integral with the extra gate
  std::vector<double> logHL, logHR;      // per candidate: child log exposures
  std::vector<char> true_exposure_ready;
  std::vector<int> proxy_count;
  std::vector<double> proxy_logHL, proxy_logHR; // hard surrogate only; never used by target
  std::vector<double> score, logprior, logq;   // index 0 = stop, then 1 + candidate
  double logQA = 0.0, logPhi = 0.0, log_inc = 0.0;
  std::unordered_map<int, std::shared_ptr<SoftExactAxis> > exact_axis;
  std::vector<double> count_norm;       // candidate-specific coefficient log normalizers
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
  // Only used when a log-left probability becomes subnormal or zero. Keep the data
  // and O(d) parameters, rather than a third O(n * candidates) probability array.
  const arma::mat* points = nullptr;
  arma::vec gate, width;
  struct Value { double logleft, left; };
  bool lazy = false;               // dense fallback supports existing internal callers
  size_t max_axis_entries = 65536;
  mutable std::vector<std::unordered_map<size_t, Value> > sparse;
  mutable int last_i = -1, last_c = -1;
  mutable Value last_value{0.0, 0.0};
  mutable size_t evaluations = 0, cache_hits = 0;
  void invalidate_axis(int j) {
    if ((int)sparse.size() > j) sparse[j].clear();
    if (last_c >= 0 && cand_axis[last_c] == j) last_i = last_c = -1;
  }
  Value value(int i, int c) const {
    if (i == last_i && c == last_c) { ++cache_hits; return last_value; }
    const int j = cand_axis[c];
    auto& entries = sparse[j];
    const size_t key = (size_t)i * (offset[j + 1] - offset[j]) + c - offset[j];
    const auto found = entries.find(key);
    Value out;
    if (found != entries.end()) {
      out = found->second; ++cache_hits;
    } else {
      const double z = gate[j] * ((*points)(i, j) - cand_cut[c]) / width[j];
      out.logleft = pst_logistic_log_right(-z); out.left = std::exp(out.logleft);
      ++evaluations;
      if (max_axis_entries > 0) {
        // Eviction affects work only. Never retain references into this map.
        if (entries.size() >= max_axis_entries) entries.clear();
        entries.emplace(key, out);
      }
    }
    last_i = i; last_c = c; last_value = out;
    return out;
  }
  double left(int i, int c) const {
    return lazy ? value(i, c).logleft : loglft[(size_t)i * n_cand + c];
  }
  double r(int i, int c) const {
    return lazy ? value(i, c).left : lft[(size_t)i * n_cand + c];
  }
  double right(int i, int c) const {
    double x = left(i, c);
    if (x >= -std::numeric_limits<double>::min() && points != nullptr) {
      const int j = cand_axis[c];
      return pst_logistic_log_right(gate[j] * ((*points)(i, j) - cand_cut[c]) / width[j]);
    }
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
  bool hard_proposal = false;     // hard surrogate only above exact_max
  double proposal_temperature = 0.5, proposal_defensive = 0.1;
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

// A reference count distribution for all cuts on one axis. Compute both log
// probabilities from the log odds, so rare tails survive even when expit rounds
// to 0 or 1. The array is updated backwards and uses O(m) working memory.
static inline std::vector<double> soft_pb_logpmf_logodds(const std::vector<double>& logodds) {
  const double ninf = -std::numeric_limits<double>::infinity();
  std::vector<double> p(logodds.size() + 1, ninf);
  p[0] = 0.0;
  for (size_t s = 0; s < logodds.size(); ++s) {
    const double l1 = pst_logistic_log_right(logodds[s]);
    const double l0 = pst_logistic_log_right(-logodds[s]);
    for (size_t k = s + 1; k > 0; --k)
      p[k] = soft_logaddexp(l0 + p[k], l1 + p[k - 1]);
    p[0] += l0;
  }
  return p;
}

// Psi(j, c) = E_B[Q(N, H_l) Q(m - N, H_r)] for m > exact_max, through the exact
// Beta-integral representation
//   Psi = b^{2a} Gamma(2a + m) / Gamma(a)^2 * c_l^a c_r^a c_r^m
//         * int_0^1 u^{a-1} (1-u)^{a-1} prod_i (alpha_i + beta_i u) du,
//   c_l = 1/(b + H_l), c_r = 1/(b + H_r), alpha_i = 1 - r_i, beta_i = r_i c_l/c_r - (1 - r_i),
// whose product of affine factors is log-concave in u. The product is integrated
// by a Laplace
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

// Exact exponential tilting across logistic cut locations. For a reference
// cut s0, delta = gate_j * (s - s0) / width_j and p_s(k) is proportional to
// p_s0(k) exp(k * delta). Exposures remain specific to the candidate cut.
static inline double soft_tilt_log_psi(const std::vector<double>& logpmf,
    double delta, const SoftModel& M, double lbL, double lbR) {
  const int m = (int)logpmf.size() - 1;
  std::vector<double> tilted(m + 1), joint(m + 1);
  for (int k = 0; k <= m; ++k) {
    tilted[k] = logpmf[k] + (k == 0 ? 0.0 : k * delta);
    joint[k] = tilted[k] + M.logQ_lb(k, lbL) + M.logQ_lb(m - k, lbR);
  }
  return soft_lse(joint) - soft_lse(tilted);
}

// LSE_k(v[k] + k * slope), without allocating a candidate-specific vector.
static inline double soft_lse_affine(const std::vector<double>& v, double slope) {
  double mx = -std::numeric_limits<double>::infinity();
  for (size_t k = 0; k < v.size(); ++k) {
    const double x = v[k] + (k == 0 ? 0.0 : k * slope);
    if (x > mx) mx = x;
  }
  if (!std::isfinite(mx)) return mx;
  double sum = 0.0;
  for (size_t k = 0; k < v.size(); ++k)
    sum += std::exp(v[k] + (k == 0 ? 0.0 : k * slope) - mx);
  return mx + std::log(sum);
}

static inline std::shared_ptr<SoftExactAxis> soft_exact_axis(
    const std::vector<double>& logits, double refcut, double gate, double width,
    const SoftModel& M) {
  std::shared_ptr<SoftExactAxis> out = std::make_shared<SoftExactAxis>();
  out->refcut = refcut; out->gate = gate; out->width = width;
  out->logits = logits;
  const int m = logits.size();
  out->log_left.resize(m); out->log_right.resize(m);
  out->logpmf.assign(m + 1, -std::numeric_limits<double>::infinity());
  out->logpmf[0] = 0.0;
  for (int s = 0; s < m; ++s) {
    if (!std::isfinite(logits[s])) out->finite = false;
    const double l1 = out->log_left[s] = pst_logistic_log_right(logits[s]);
    const double l0 = out->log_right[s] = pst_logistic_log_right(-logits[s]);
    for (int k = s + 1; k > 0; --k)
      out->logpmf[k] = soft_logaddexp(l0 + out->logpmf[k], l1 + out->logpmf[k - 1]);
    out->logpmf[0] += l0;
  }
  out->coef.resize(m + 1);
  for (int k = 0; k <= m; ++k)
    out->coef[k] = out->logpmf[k] + M.lgam[k] + M.lgam[m - k];
  return out;
}

// The coefficient array A(k) = log p_0(k) + log Gamma(a+k) +
// log Gamma(a+m-k) is independent of the candidate.  The returned normalizer
// is reused when sampling K for a selected action.
static inline double soft_exact_log_psi(const SoftExactAxis& axis, double delta,
    const SoftModel& M, double lbL, double lbR, double* log_count_norm = nullptr) {
  const int m = (int)axis.logpmf.size() - 1;
  const double norm = soft_lse_affine(axis.coef, delta + lbR - lbL);
  if (log_count_norm != nullptr) *log_count_norm = norm;
  const double constant = 2.0 * (M.a_log_b - M.lgam[0]) -
    M.a * lbL - (M.a + m) * lbR;
  return constant + norm - soft_lse_affine(axis.logpmf, delta);
}

static inline int soft_exact_sample_count(const SoftExactAxis& axis, double delta,
    double lbL, double lbR,
    double log_count_norm = std::numeric_limits<double>::quiet_NaN()) {
  const double slope = delta + lbR - lbL;
  if (!std::isfinite(log_count_norm))
    log_count_norm = soft_lse_affine(axis.coef, slope);
  const double u = R::unif_rand();
  double cumulative = 0.0;
  for (size_t k = 0; k < axis.coef.size(); ++k) {
    cumulative += std::exp(axis.coef[k] + (k == 0 ? 0.0 : k * slope) - log_count_norm);
    if (u <= cumulative) return (int)k;
  }
  return (int)axis.coef.size() - 1;
}

static inline bool soft_exact_prepare_prefix(SoftExactAxis& axis) {
  if (!axis.prefix.empty()) return false;
  const size_t m = axis.logits.size();
  axis.prefix.resize(m + 1);
  axis.prefix[0] = {0.0};
  const double ninf = -std::numeric_limits<double>::infinity();
  for (size_t s = 1; s <= m; ++s) {
    axis.prefix[s].resize(s + 1);
    for (size_t k = 0; k <= s; ++k)
      axis.prefix[s][k] = soft_logaddexp(
        k < s ? axis.log_right[s - 1] + axis.prefix[s - 1][k] : ninf,
        k > 0 ? axis.log_left[s - 1] + axis.prefix[s - 1][k - 1] : ninf);
  }
  return true;
}

// Draw B | sum(B)=k using the coordinate's reference probabilities.  The
// common cut shift cancels in this conditional law.  No proposal probability
// is evaluated because the exact allocation factor cancels from the weight.
// Return whether this draw built the lazy prefix table.
static inline bool soft_exact_draw_bits(SoftExactAxis& axis, int k,
    std::vector<char>& B) {
  const int m = axis.logits.size();
  B.resize(m);
  if (k == 0 || k == m) {
    std::fill(B.begin(), B.end(), k == m ? 1 : 0);
    return false;
  }
  const bool built = soft_exact_prepare_prefix(axis);
  for (int s = m; s > 0; --s) {
    if (k == 0 || k == s) {
      std::fill(B.begin(), B.begin() + s, k == s ? 1 : 0);
      break;
    }
    const double l1 = axis.log_left[s - 1] + axis.prefix[s - 1][k - 1];
    const double l0 = axis.log_right[s - 1] + axis.prefix[s - 1][k];
    B[s - 1] = std::log(R::unif_rand()) < l1 - soft_logaddexp(l0, l1) ? 1 : 0;
    k -= B[s - 1];
  }
  return built;
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

// Two-log-probability version used for exact allocation proposals/replay.
static inline std::vector<std::vector<double> > soft_pb_table_logprobs(
    const std::vector<double>& log_left, const std::vector<double>& log_right) {
  const double ninf = -std::numeric_limits<double>::infinity();
  const size_t m = log_left.size();
  std::vector<std::vector<double> > tab(m + 1, std::vector<double>(m + 1, ninf));
  tab[0][0] = 0.0;
  for (size_t s = 1; s <= m; ++s)
    for (size_t k = 0; k <= s; ++k)
      tab[s][k] = soft_logaddexp(log_right[s - 1] + tab[s - 1][k],
          k > 0 ? log_left[s - 1] + tab[s - 1][k - 1] : ninf);
  return tab;
}

static inline double soft_cond_bernoulli_logprobs(
    const std::vector<std::vector<double> >& tab,
    const std::vector<double>& log_left, const std::vector<double>& log_right,
    int k, std::vector<char>& B, bool replay) {
  const double ninf = -std::numeric_limits<double>::infinity();
  double logq = 0.0;
  for (int s = (int)log_left.size(); s >= 1; --s) {
    // Normalize the two branches directly: subtracting a probability rounded
    // to one would erase rare right allocations in a forced reference path.
    const double l1 = k > 0 ? log_left[s - 1] + tab[s - 1][k - 1] : ninf;
    const double l0 = k < s ? log_right[s - 1] + tab[s - 1][k] : ninf;
    const double norm = soft_logaddexp(l0, l1);
    if (!replay) B[s - 1] = std::log(R::unif_rand()) < l1 - norm ? 1 : 0;
    logq += (B[s - 1] ? l1 : l0) - norm;
    k -= B[s - 1];
  }
  return logq;
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
