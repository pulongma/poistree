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
  int side;
};

struct SoftNodeP {
  bool active = false;
  int S = -1;
  int J = -1;
  double L = NA_REAL;
  int depth = 0;
  int m = 0;
  double logH = 0.0;
  std::vector<SoftGateStep> path;
};

struct SoftExactAxis {
  double refcut = 0.0, gate = 0.0, width = 1.0;
  bool finite = true;
  std::vector<double> logits, log_left, log_right, logpmf, coef;
  std::vector<std::vector<double> > prefix;
};

struct SoftPathNode {
  int heap = 1, depth = 0, parent = -1;
  SoftGateStep gate{-1, 0.0, 0};
  std::vector<int> pts;
  std::vector<double> axisH;
  double logH = 0.0;
  bool expanded = false;
  bool hard_scored = false;
  std::vector<double> axisL, axisR;
  std::vector<double> logHL, logHR;
  std::vector<char> true_exposure_ready;
  std::vector<int> proxy_count;
  std::vector<double> proxy_logHL, proxy_logHR;
  std::vector<double> score, logprior, logq;
  double logQA = 0.0, logPhi = 0.0, log_inc = 0.0;
  std::unordered_map<int, std::shared_ptr<SoftExactAxis> > exact_axis;
  std::vector<double> count_norm;
  std::unordered_multimap<unsigned long long, std::pair<int, int> > children;
  double glued = 0.0;
  int glued_stamp = -1;
};

struct SoftRecord {
  int node;
  signed char S;
  int act;
};

struct SoftParticleS {
  std::map<int, SoftRecord> rec;
  double logw = 0.0;
};

struct SoftGateTable {
  int n_cand = 0;
  std::vector<int> offset, cand_axis;
  std::vector<double> cand_cut, loglft, lft;

  const arma::mat* points = nullptr;
  arma::vec gate, width;
  struct Value { double logleft, left; };
  bool lazy = false;
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

struct SoftRef {
  std::vector<SoftNodeP> node;
  std::vector<int> z;
  std::unordered_map<long long, char> bits;
  bool valid = false;
};

struct SoftModel {
  arma::mat X;
  arma::mat region;
  double a, b;
  double rho, eta;
  int Dmax;
  int n_nodes;
  std::vector<std::vector<double> > grid;
  int exact_max;
  double defensive;
  bool hard_proposal = false;
  double proposal_temperature = 0.5, proposal_defensive = 0.1;
  std::vector<double> lgam;
  double a_log_b;

  void init_tables() {
    lgam.resize(n() + 1);
    for (int k = 0; k <= n(); ++k) lgam[k] = R::lgammafn(a + k);
    a_log_b = a * std::log(b);
  }
  // Return the collapsed Gamma-Poisson leaf log marginal.
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

// Return the log left-gate probability for a split.
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

// Return a coordinate path log integral, optionally appending a gate.

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

// Return the Poisson-binomial log PMF from log Bernoulli probabilities.

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

// Return the Poisson-binomial log PMF from Bernoulli log odds.

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

// Approximate a large-node split score by the Laplace method.

static inline double soft_log_psi_laplace(const SoftModel& M, const SoftGateTable& G,
                                          const std::vector<int>& pts, int c,
                                          double lbL, double lbR) {
  const int m = pts.size();
  if (m == 0) return M.logQ_lb(0, lbL) + M.logQ_lb(0, lbR);
  const double rho = std::exp(lbR - lbL);
  std::vector<double> alpha(m), beta(m);
  double mu = 0.0;
  for (int s = 0; s < m; ++s) {
    const double r = G.r(pts[s], c);
    alpha[s] = 1.0 - r; beta[s] = r * rho - (1.0 - r); mu += r;
  }
  auto slope = [&](double u, double& g1, double& g2) {
    g1 = 0.0; g2 = 0.0;
    for (int s = 0; s < m; ++s) { const double q = beta[s] / (alpha[s] + beta[s] * u); g1 += q; g2 -= q * q; }
  };
  double g1, g2, lo = 0.0, hi = 1.0, u = std::min(1.0 - 1e-9, std::max(1e-9, mu / m));
  slope(0.0, g1, g2);
  if (g1 <= 0.0) u = 0.0;
  else {
    slope(1.0, g1, g2);
    if (g1 >= 0.0) u = 1.0;
    else for (int it = 0; it < 30; ++it) {
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

  const double kappa = std::max(-g2, 1e-12), sd = 1.0 / std::sqrt(kappa), shift = g1 / kappa;
  const double zlo = (0.0 - u - shift) / sd, zhi = (1.0 - u - shift) / sd;

  double logP;
  if (zlo > 0.0) { const double a1 = R::pnorm(zlo, 0.0, 1.0, 0, 1), a2 = R::pnorm(zhi, 0.0, 1.0, 0, 1); logP = a1 + std::log1p(-std::exp(a2 - a1)); }
  else if (zhi < 0.0) { const double a1 = R::pnorm(zhi, 0.0, 1.0, 1, 1), a2 = R::pnorm(zlo, 0.0, 1.0, 1, 1); logP = a1 + std::log1p(-std::exp(a2 - a1)); }
  else logP = std::log(R::pnorm(zhi, 0.0, 1.0, 1, 0) - R::pnorm(zlo, 0.0, 1.0, 1, 0));
  const double logZ = 0.5 * std::log(2.0 * M_PI / kappa) + 0.5 * g1 * g1 / kappa + logP;
  const double width = std::min(0.25, std::min(sd, 1.0 / std::max(std::abs(g1), 1e-300)));
  const double ut = std::min(1.0 - width, std::max(width, u));
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

// Compute an exact split score by tilting a reference count distribution.

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

// Return log sum exp(v[k] + k * slope).
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

// Compute an exact split score and its count-sampling normalizer.

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

// Draw allocations conditional on their total; report whether a table was built.

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

// Build the Poisson-binomial prefix table in log space.

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

// Build a prefix table from left and right log probabilities.
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

    const double l1 = k > 0 ? log_left[s - 1] + tab[s - 1][k - 1] : ninf;
    const double l0 = k < s ? log_right[s - 1] + tab[s - 1][k] : ninf;
    const double norm = soft_logaddexp(l0, l1);
    if (!replay) B[s - 1] = std::log(R::unif_rand()) < l1 - norm ? 1 : 0;
    logq += (B[s - 1] ? l1 : l0) - norm;
    k -= B[s - 1];
  }
  return logq;
}

// Draw or replay conditional Bernoulli allocations and return their log probability.

static inline double soft_cond_bernoulli(const std::vector<std::vector<double> >& tab,
                                         const std::vector<double>& logr, int k,
                                         std::vector<char>& B, bool replay) {
  int m = logr.size(); double logq = 0.0;
  for (int s = m; s >= 1; --s) {
    double logp1 = k > 0 ? std::min(0.0, logr[s - 1] + tab[s - 1][k - 1] - tab[s][k]) : -std::numeric_limits<double>::infinity();
    if (k == s) logp1 = 0.0;
    if (!replay) B[s - 1] = std::log(R::unif_rand()) < logp1 ? 1 : 0;
    logq += B[s - 1] ? logp1 : (logp1 > -0.6931471805599453 ? std::log(-std::expm1(logp1)) : std::log1p(-std::exp(logp1)));
    k -= B[s - 1];
  }
  return logq;
}

#endif
