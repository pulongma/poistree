// ============================================================================
// ppt_state.cpp -- post-hoc intensity evaluation from serialized posterior
// state draws.
//
// Every RJ-MCMC backend stores, for each retained draw, one matrix with a row
// per node and columns (heap id, axis, cut, lambda, reserved, m), plus the gate
// vector of the draw for soft fits.  Given that state, the intensity
//   lambda(s) = sum_v lambda_v phi_v(s)
// is reproducible EXACTLY at arbitrary locations: the node boxes follow from
// the heap topology and the cuts, hard routing uses "left iff x_j < cut", and
// the soft bases are products of the ancestor gates (recursive logistic or
// compact cubic), so nothing beyond this state is needed.  This is what frees
// `ppt_predict()`, `ppt_marginal()`, and `ppt_lppd()` from the requirement
// that every location of interest be listed in `predict_at` before fitting.
//
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <limits>
#include "soft_logistic.h"
using namespace Rcpp;

namespace pstate {

struct SNode {
  long hid;
  int axis;          // -1 leaf
  double cut;
  double lambda;
  arma::mat box;     // d x 2, derived from the heap topology
};

struct SGate {
  int axis;
  double cut;
  int side;          // -1 left, +1 right
  double eff_width;  // compact family: parent width / (1+depth)^gate_depth
};

static inline int heap_depth(long hid) {
  int depth = 0;
  while (hid > 1) { hid >>= 1; depth++; }
  return depth;
}

// compact (cubic) gate value; h = eff_width / gate
static inline double compact_value(const SGate&g, double x, double gate) {
  double h = g.eff_width / gate;
  double t = (x - (g.cut - h)) / (2.0 * h);
  double right = t <= 0.0 ? 0.0 : (t >= 1.0 ? 1.0 : t * t * (3.0 - 2.0 * t));
  return g.side < 0 ? 1.0 - right : right;
}

static std::vector<double> compact_poly_multiply(
    const std::vector<double>&a, const std::vector<double>&b) {
  std::vector<double> out(a.size() + b.size() - 1, 0.0);
  for (size_t i = 0; i < a.size(); i++)
    for (size_t j = 0; j < b.size(); j++) out[i + j] += a[i] * b[j];
  return out;
}

// Exact integral of the product of the compact-cubic gates on one axis.
// Splitting at every transition endpoint makes each factor polynomial on
// every segment, so their product can be integrated coefficient by
// coefficient.
static double compact_axis_integral(const std::vector<SGate>&path, int axis,
                                    double dom_lo, double dom_hi,
                                    double gate) {
  std::vector<const SGate*> gates;
  std::vector<double> breaks = {dom_lo, dom_hi};
  for (const SGate&g : path) if (g.axis == axis) {
    gates.push_back(&g);
    double h = g.eff_width / gate;
    if (g.cut - h > dom_lo && g.cut - h < dom_hi)
      breaks.push_back(g.cut - h);
    if (g.cut + h > dom_lo && g.cut + h < dom_hi)
      breaks.push_back(g.cut + h);
  }
  if (gates.empty()) return dom_hi - dom_lo;

  std::sort(breaks.begin(), breaks.end());
  std::vector<double> uniq;
  uniq.reserve(breaks.size());
  for (double z : breaks)
    if (uniq.empty() ||
        std::abs(z - uniq.back()) > 1e-13 * (1.0 + std::abs(z)))
      uniq.push_back(z);

  double total = 0.0;
  for (size_t s = 0; s + 1 < uniq.size(); s++) {
    double lo = uniq[s], hi = uniq[s + 1], dx = hi - lo;
    if (!(dx > 0.0)) continue;
    std::vector<double> poly(1, 1.0);
    bool zero = false;
    for (const SGate*gp : gates) {
      const SGate&g = *gp;
      double h = g.eff_width / gate;
      double tr_lo = g.cut - h, tr_hi = g.cut + h;
      if (hi <= tr_lo) {
        if (g.side > 0) { zero = true; break; }
        continue;
      }
      if (lo >= tr_hi) {
        if (g.side < 0) { zero = true; break; }
        continue;
      }
      double t0 = (lo - tr_lo) / (2.0 * h);
      double q = dx / (2.0 * h);
      double t02 = t0 * t0, t03 = t02 * t0;
      double q2 = q * q, q3 = q2 * q;
      std::vector<double> factor(4);
      factor[0] = 3.0 * t02 - 2.0 * t03;
      factor[1] = 6.0 * t0 * q - 6.0 * t02 * q;
      factor[2] = 3.0 * q2 - 6.0 * t0 * q2;
      factor[3] = -2.0 * q3;
      if (g.side < 0) {
        factor[0] = 1.0 - factor[0];
        for (size_t k = 1; k < factor.size(); k++) factor[k] = -factor[k];
      }
      poly = compact_poly_multiply(poly, factor);
    }
    if (zero) continue;
    double integral = 0.0;
    for (size_t k = 0; k < poly.size(); k++)
      integral += poly[k] / (double)(k + 1);
    total += dx * integral;
  }
  if (total < 0.0 && total > -1e-12 * (dom_hi - dom_lo)) total = 0.0;
  return std::max(total, std::numeric_limits<double>::min());
}

static double soft_axis_value(const std::vector<SGate>&path, int axis,
                              double x, const arma::mat&region,
                              const arma::vec&gate, int gate_mode) {
  double log_value = 0.0;
  for (const SGate&g : path) if (g.axis == axis) {
    if (gate_mode == 2) {
      double width = region(axis, 1) - region(axis, 0);
      double z = gate[axis] * (x - g.cut) / width;
      log_value += g.side < 0 ? pst_logistic_log_right(-z)
                              : pst_logistic_log_right(z);
    } else {
      double q = compact_value(g, x, gate[axis]);
      if (q <= 0.0) return 0.0;
      log_value += std::log(q);
    }
  }
  return log_value < -745.0 ? 0.0 : std::exp(log_value);
}

static double soft_axis_integral(const std::vector<SGate>&path, int axis,
                                 const arma::mat&region,
                                 const arma::vec&gate, int gate_mode) {
  if (gate_mode == 1)
    return compact_axis_integral(
      path, axis, region(axis, 0), region(axis, 1), gate[axis]
    );
  std::vector<double> cuts;
  std::vector<int> sides;
  for (const SGate&g : path) if (g.axis == axis) {
    cuts.push_back(g.cut);
    sides.push_back(g.side);
  }
  return pst_logistic_path_axis_integral(
    cuts, sides, region(axis, 0), region(axis, 1), gate[axis]
  );
}

}  // namespace pstate

// [[Rcpp::export]]
arma::mat ppt_eval_state(List state_nodes, arma::mat state_gate,
                         arma::mat region, arma::mat newdata,
                         int gate_mode, double gate_depth) {
  using namespace pstate;
  const int S = state_nodes.size();
  const arma::uword np = newdata.n_rows, d = newdata.n_cols;
  if (region.n_rows != d || region.n_cols != 2)
    stop("`region` must be a d by 2 matrix matching `newdata`.");
  arma::mat out(np, S, arma::fill::zeros);

  for (int s = 0; s < S; s++) {
    arma::mat M = as<arma::mat>(state_nodes[s]);
    if (M.n_cols < 4) stop("malformed state draw: need >= 4 columns");
    const int nn = M.n_rows;

    // sort rows by heap id so parents precede children
    arma::uvec ord = arma::sort_index(M.col(0));
    std::vector<SNode> nd(nn);
    std::unordered_map<long, int> at;
    at.reserve(2 * nn);
    for (int r = 0; r < nn; r++) {
      const arma::uword j = ord[r];
      SNode z;
      z.hid = (long)M(j, 0);
      z.axis = (int)M(j, 1);
      z.cut = M(j, 2);
      z.lambda = std::isfinite(M(j, 3)) ? M(j, 3) : 0.0;
      nd[r] = z;
      at[z.hid] = r;
    }
    if (!at.count(1L)) stop("state draw has no root node");
    // derive boxes from the topology
    for (int r = 0; r < nn; r++) {
      SNode&z = nd[r];
      if (z.hid == 1L) { z.box = region; continue; }
      auto ip = at.find(z.hid / 2);
      if (ip == at.end()) stop("state draw has an orphan node");
      const SNode&pa = nd[ip->second];
      if (pa.axis < 0) stop("state draw parent is a leaf");
      z.box = pa.box;
      if (z.hid % 2 == 0) z.box(pa.axis, 1) = pa.cut;
      else                z.box(pa.axis, 0) = pa.cut;
    }
    // gate vector of the draw (soft only)
    arma::vec gate;
    if (gate_mode != 0) {
      if ((int)state_gate.n_rows <= s || state_gate.n_cols != d)
        stop("`state_gate` must be draws by d for soft fits");
      gate = state_gate.row(s).t();
      if (!gate.is_finite() || arma::any(gate <= 0.0))
        stop("soft state draw has an invalid gate vector");
    }
    // contributing nodes and, for soft fits, their ancestor gate paths
    std::vector<int> contrib;
    contrib.reserve(nn);
    for (int r = 0; r < nn; r++)
      if (nd[r].axis < 0) contrib.push_back(r);
    std::vector<std::vector<SGate> > paths;
    if (gate_mode != 0) {
      paths.resize(contrib.size());
      for (size_t k = 0; k < contrib.size(); k++) {
        long h = nd[contrib[k]].hid;
        while (h > 1) {
          const SNode&pa = nd[at[h / 2]];
          SGate g;
          g.axis = pa.axis;
          g.cut = pa.cut;
          g.side = (h % 2 == 0) ? -1 : 1;
          g.eff_width = (pa.box(pa.axis, 1) - pa.box(pa.axis, 0)) /
            std::pow(1.0 + (double)heap_depth(pa.hid), gate_depth);
          paths[k].push_back(g);
          h /= 2;
        }
      }
    }

    for (arma::uword i = 0; i < np; i++) {
      const arma::rowvec x = newdata.row(i);
      double val = 0.0;
      if (gate_mode == 0) {
        // Hard routing takes the rate of the terminal leaf.
        long h = 1;
        for (;;) {
          const SNode&z = nd[at[h]];
          if (z.axis < 0) { val = z.lambda; break; }
          h = (x[z.axis] < z.cut) ? 2 * h : 2 * h + 1;
          if (!at.count(h)) stop("state draw routing reached a missing node");
        }
      } else {
        for (size_t k = 0; k < contrib.size(); k++) {
          const SNode&z = nd[contrib[k]];
          if (z.lambda <= 0.0) continue;
          double lphi = 0.0;
          bool zero = false;
          for (const SGate&g : paths[k]) {
            if (gate_mode == 2) {
              double width = region(g.axis, 1) - region(g.axis, 0);
              double zz = gate[g.axis] * (x[g.axis] - g.cut) / width;
              lphi += g.side < 0 ? pst_logistic_log_right(-zz)
                                 : pst_logistic_log_right(zz);
            } else {
              double q = compact_value(g, x[g.axis], gate[g.axis]);
              if (q <= 0.0) { zero = true; break; }
              lphi += std::log(q);
            }
          }
          if (!zero && lphi > -745.0) val += z.lambda * std::exp(lphi);
        }
      }
      out(i, s) = val;
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::mat ppt_marginal_state(List state_nodes, arma::mat state_gate,
                             arma::mat region, arma::vec grid,
                             int variable, int gate_mode,
                             double gate_depth, bool average) {
  using namespace pstate;
  const int S = state_nodes.size();
  const arma::uword d = region.n_rows;
  if (S <= 0) stop("`state_nodes` must contain at least one draw.");
  if (region.n_cols != 2 || d == 0 || !region.is_finite())
    stop("`region` must be a finite d by 2 matrix.");
  for (arma::uword j = 0; j < d; j++)
    if (!(region(j, 1) > region(j, 0)))
      stop("Every fitted-region upper bound must exceed its lower bound.");
  if (variable < 0 || variable >= (int)d)
    stop("`variable` is outside the fitted dimension.");
  if (!grid.is_finite()) stop("`grid` must be finite.");
  if (gate_mode < 0 || gate_mode > 2)
    stop("`gate_mode` must be 0 (hard), 1 (compact), or 2 (logistic).");
  if (!std::isfinite(gate_depth) || gate_depth < 0.0)
    stop("`gate_depth` must be finite and nonnegative.");

  double normalizer = 1.0;
  if (average) for (arma::uword j = 0; j < d; j++)
    if ((int)j != variable)
      normalizer *= region(j, 1) - region(j, 0);
  arma::mat out(grid.n_elem, S, arma::fill::zeros);

  for (int s = 0; s < S; s++) {
    arma::mat M = as<arma::mat>(state_nodes[s]);
    if (M.n_cols < 4 || M.n_rows == 0)
      stop("malformed state draw: need nonempty rows and >= 4 columns");
    const int nn = M.n_rows;

    arma::uvec ord = arma::sort_index(M.col(0));
    std::vector<SNode> nd(nn);
    std::unordered_map<long, int> at;
    at.reserve(2 * nn);
    for (int r = 0; r < nn; r++) {
      const arma::uword j = ord[r];
      SNode z;
      z.hid = (long)M(j, 0);
      z.axis = (int)M(j, 1);
      z.cut = M(j, 2);
      z.lambda = M(j, 3);
      if (z.hid < 1 || z.axis < -1 || z.axis >= (int)d)
        stop("state draw contains an invalid node identifier or axis");
      if (!std::isfinite(z.lambda) || z.lambda < 0.0)
        stop("state draw contains an invalid node rate");
      nd[r] = z;
      at[z.hid] = r;
    }
    if (!at.count(1L)) stop("state draw has no root node");
    for (int r = 0; r < nn; r++) {
      SNode&z = nd[r];
      if (z.hid == 1L) { z.box = region; continue; }
      auto ip = at.find(z.hid / 2);
      if (ip == at.end()) stop("state draw has an orphan node");
      const SNode&pa = nd[ip->second];
      if (pa.axis < 0 || !std::isfinite(pa.cut))
        stop("state draw parent is a leaf or has an invalid cut");
      z.box = pa.box;
      if (z.hid % 2 == 0) z.box(pa.axis, 1) = pa.cut;
      else                z.box(pa.axis, 0) = pa.cut;
    }

    arma::vec gate;
    if (gate_mode != 0) {
      if ((int)state_gate.n_rows <= s || state_gate.n_cols != d)
        stop("`state_gate` must be draws by d for soft fits");
      gate = state_gate.row(s).t();
      if (!gate.is_finite() || arma::any(gate <= 0.0))
        stop("soft state draw has an invalid gate vector");
    }

    std::vector<int> contrib;
    contrib.reserve(nn);
    for (int r = 0; r < nn; r++)
      if (nd[r].axis < 0) contrib.push_back(r);
    std::vector<std::vector<SGate> > paths(contrib.size());
    if (gate_mode != 0) {
      for (size_t k = 0; k < contrib.size(); k++) {
        long h = nd[contrib[k]].hid;
        while (h > 1) {
          auto ip = at.find(h / 2);
          if (ip == at.end()) stop("state draw has an orphan path");
          const SNode&pa = nd[ip->second];
          SGate g;
          g.axis = pa.axis;
          g.cut = pa.cut;
          g.side = (h % 2 == 0) ? -1 : 1;
          g.eff_width = (pa.box(pa.axis, 1) - pa.box(pa.axis, 0)) /
            std::pow(1.0 + (double)heap_depth(pa.hid), gate_depth);
          paths[k].push_back(g);
          h /= 2;
        }
      }
    }

    for (size_t k = 0; k < contrib.size(); k++) {
      const SNode&z = nd[contrib[k]];
      if (z.lambda == 0.0) continue;
      double integrated_other_axes = 1.0;
      for (arma::uword j = 0; j < d; j++) if ((int)j != variable) {
        integrated_other_axes *= gate_mode == 0
          ? z.box(j, 1) - z.box(j, 0)
          : soft_axis_integral(paths[k], (int)j, region, gate, gate_mode);
      }
      integrated_other_axes /= normalizer;

      if (gate_mode == 0) {
        double upper = region(variable, 1);
        double tolerance = std::sqrt(std::numeric_limits<double>::epsilon()) *
          std::max(1.0, std::abs(upper));
        bool is_last = std::abs(z.box(variable, 1) - upper) <= tolerance;
        for (arma::uword q = 0; q < grid.n_elem; q++) {
          bool active = grid[q] >= z.box(variable, 0) &&
            (grid[q] < z.box(variable, 1) ||
             (is_last && grid[q] <= z.box(variable, 1) + tolerance));
          if (active)
            out(q, s) += z.lambda * integrated_other_axes;
        }
      } else {
        for (arma::uword q = 0; q < grid.n_elem; q++)
          out(q, s) += z.lambda * integrated_other_axes *
            soft_axis_value(
              paths[k], variable, grid[q], region, gate, gate_mode
            );
      }
    }
  }
  return out;
}
