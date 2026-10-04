// Shared-path SMC for the hard terminal-leaf PPT; see hard_smc.h and
// CLAUDE/shared_path_smc_hard_ppt.tex.  The sampler follows Ma (2026),
// Sect. 3: breadth-first growth, one-step lookahead proposal, adaptive
// multinomial resampling, unbiased evidence increments.
#ifndef _USE_Armadillo
#define _USE_Armadillo
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#endif
#include "hard_smc.h"

static inline double hard_lse(const std::vector<double>& v) {
  double mx = -std::numeric_limits<double>::infinity();
  for (double x : v) if (x > mx) mx = x;
  if (!std::isfinite(mx)) return mx;
  double s = 0.0;
  for (double x : v) s += std::exp(x - mx);
  return mx + std::log(s);
}

static inline int hard_sample_log(const std::vector<double>& logp) {
  const double norm = hard_lse(logp), u = R::unif_rand();
  double cum = 0.0;
  for (size_t k = 0; k < logp.size(); ++k) {
    cum += std::exp(logp[k] - norm);
    if (u <= cum) return (int)k;
  }
  return (int)logp.size() - 1;
}

// box constraints of PPT::find_valid_cuts: minimal side 1e-2, aspect ratio cap
static inline bool hard_good_shape(const arma::mat& box, double rmax) {
  arma::vec side = box.col(1) - box.col(0);
  if (side.min() < 1e-2) return false;
  if (std::isinf(rmax) && rmax > 0) return true;
  if (!std::isfinite(rmax) || rmax < 1.0) return false;
  return side.max() / side.min() <= rmax;
}

int HardSMCtree::make_root() {
  HardNode root;
  root.box = region;
  root.pts.resize(n);
  for (int i = 0; i < n; ++i) root.pts[i] = i;
  root.m = n;
  root.area = arma::prod(region.col(1) - region.col(0));
  store.push_back(std::move(root));
  return 0;
}

int HardSMCtree::make_child(int parent, int cand, int side, std::vector<int>&& pts) {
  const HardNode& par = store[parent];
  HardNode c;
  c.heap = 2 * par.heap + (side < 0 ? 0 : 1);
  c.depth = par.depth + 1;
  c.parent = parent;
  c.cand = cand;
  c.box = par.box;
  c.box(par.cand_axis[cand], side < 0 ? 1 : 0) = par.cand_cut[cand];
  c.m = pts.size();
  c.pts = std::move(pts);
  c.area = arma::prod(c.box.col(1) - c.box.col(0));
  store.push_back(std::move(c));
  return (int)store.size() - 1;
}

// Candidate cuts and scores of one box (the one-step lookahead ingredients of
// Ma 2026, eqs. (4)-(6)), computed once per distinct node.  Per axis: sort the
// coordinate, take type-1 quantiles of the distinct interior values at
// cut_grid_n probabilities in [0.05, 0.95], count the left observations by
// binary search; a node with no valid cut stops with probability one.
void HardSMCtree::expand_node(int v) {
  HardNode& A = store[v];
  if (A.expanded) return;
  const double buffer = 1e-3, rho_d = rho_depth(A.depth);
  A.logQ0 = logQ0(A.m, A.area);
  std::vector<double> xs(A.m), probs(cut_grid_n);
  for (int q = 0; q < cut_grid_n; ++q) probs[q] = cut_grid_n == 1 ? 0.05 : 0.05 + 0.9 * q / (cut_grid_n - 1.0);
  std::vector<int> axis_count(d, 0);
  if (A.m >= 2 * min_leaf) {
    for (int j = 0; j < d; ++j) {
      for (int s = 0; s < A.m; ++s) xs[s] = X(A.pts[s], j);
      std::sort(xs.begin(), xs.end());
      const double lo = A.box(j, 0) + buffer, hi = A.box(j, 1) - buffer;
      std::vector<double> uniq;
      for (int s = 0; s < A.m; ++s)
        if (xs[s] > lo && xs[s] < hi && (uniq.empty() || xs[s] != uniq.back())) uniq.push_back(xs[s]);
      std::vector<double> cuts;
      if (!uniq.empty()) {
        for (double p : probs) {
          const int i = std::max(1, (int)std::floor((uniq.size() - 1) * p + 1.0));
          cuts.push_back(uniq[i - 1]);
        }
        std::sort(cuts.begin(), cuts.end());
        cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
      }
      auto try_cut = [&](double cut) {
        arma::mat boxL = A.box, boxR = A.box;
        boxL(j, 1) = cut; boxR(j, 0) = cut;
        if (!hard_good_shape(boxL, max_aspect) || !hard_good_shape(boxR, max_aspect)) return;
        // left child is [lo, cut) except on the last axis, where it is [lo, cut]
        // (the box convention of in_region_nd in the dense implementation)
        const int nL = (int)((j == d - 1 ? std::upper_bound(xs.begin(), xs.end(), cut)
                                         : std::lower_bound(xs.begin(), xs.end(), cut)) - xs.begin()), nR = A.m - nL;
        const double areaL = arma::prod(boxL.col(1) - boxL.col(0)), areaR = arma::prod(boxR.col(1) - boxR.col(0));
        if (nL < min_leaf || nR < min_leaf || areaL <= 0.0 || areaR <= 0.0) return;
        A.cand_axis.push_back(j); A.cand_cut.push_back(cut); A.cand_nL.push_back(nL);
        A.score.push_back(logQ0(nL, areaL) + logQ0(nR, areaR));
        ++axis_count[j];
      };
      for (double cut : cuts) try_cut(cut);
      if (axis_count[j] == 0) {                       // median fallback of find_valid_cuts
        const double med = A.m % 2 ? xs[A.m / 2] : 0.5 * (xs[A.m / 2 - 1] + xs[A.m / 2]);
        if (med > lo && med < hi) try_cut(med);
      }
    }
  }
  // prior of a split: rho_d x (1/d) x (1/K_j); the increment is action independent
  const int C = A.cand_axis.size();
  A.logq.assign(C + 1, 0.0);
  std::vector<double> full(C + 1);
  full[0] = std::log(1.0 - rho_d) + A.logQ0;
  for (int k = 0; k < C; ++k) {
    A.score[k] += std::log(rho_d) - std::log((double)d) - std::log((double)axis_count[A.cand_axis[k]]);
    full[k + 1] = A.score[k];
  }
  A.logPhi = hard_lse(full);
  for (int k = 0; k <= C; ++k) A.logq[k] = full[k] - A.logPhi;
  A.log_inc = C == 0 ? 0.0 : A.logPhi - A.logQ0;
  A.expanded = true;
  ++n_expanded;
}

// Every particle with an undecided record at heap position t draws its
// decision; children are created once per (node, candidate) and shared.
bool HardSMCtree::sample_position(int t) {
  bool advanced = false;
  std::vector<int> touched;
  for (int p = 0; p < P; ++p) {
    HardParticle& par = particles[p];
    auto it = par.rec.find(t);
    if (it == par.rec.end() || it->second.S != -1) continue;
    const int v = it->second.node;
    expand_node(v);
    touched.push_back(v);
    HardNode& A = store[v];
    advanced = true;
    par.logw += A.log_inc;
    if (A.cand_axis.empty()) { it->second.S = 0; continue; }
    const int act = hard_sample_log(A.logq);
    if (act == 0) { it->second.S = 0; continue; }
    const int c = act - 1;
    int lid, rid;
    auto ch = A.children.find(c);
    if (ch != A.children.end()) { lid = ch->second.first; rid = ch->second.second; }
    else {
      std::vector<int> left, right;
      left.reserve(A.cand_nL[c]); right.reserve(A.m - A.cand_nL[c]);
      const int j = A.cand_axis[c];
      const double cut = A.cand_cut[c];
      for (int i : A.pts) ((j == d - 1 ? X(i, j) <= cut : X(i, j) < cut) ? left : right).push_back(i);
      lid = make_child(v, c, -1, std::move(left));
      rid = make_child(v, c, +1, std::move(right));
      store[v].children.emplace(c, std::make_pair(lid, rid));
    }
    HardRecord& R = particles[p].rec[t];                 // store may have grown
    R.S = 1; R.act = c;
    particles[p].rec[2 * t] = HardRecord{lid, -1, -1};
    particles[p].rec[2 * t + 1] = HardRecord{rid, -1, -1};
  }
  // observation sets of decided nodes are no longer needed
  for (int v : touched) { std::vector<int>().swap(store[v].pts); }
  return advanced;
}

// Multinomial resampling when ESS < resample_thresh * P; offspring take the
// ancestor's record map (first by move, others by copy).
void HardSMCtree::resample() {
  std::vector<double> lw(P);
  for (int p = 0; p < P; ++p) lw[p] = particles[p].logw;
  const double norm = hard_lse(lw);
  Rcpp::NumericVector W(P);
  double sumsq = 0.0;
  for (int p = 0; p < P; ++p) { W[p] = std::exp(lw[p] - norm); sumsq += W[p] * W[p]; }
  if (1.0 / sumsq >= resample_thresh * P) return;
  ++n_resampled;
  Rcpp::IntegerVector anc = Rcpp::sample(P, P, true, W, false);
  std::vector<HardParticle> next(P);
  std::vector<int> first_use(P, -1);
  for (int p = 0; p < P; ++p) {
    const int a_ = anc[p];
    if (first_use[a_] < 0) { first_use[a_] = p; next[p].rec = std::move(particles[a_].rec); }
    else next[p].rec = next[first_use[a_]].rec;
  }
  particles.swap(next);
}

void HardSMCtree::sweep() {
  const int T = (1 << Dmax) - 1;
  store.clear();
  n_expanded = 0; n_resampled = 0;
  const int root = make_root();
  particles.assign(P, HardParticle());
  for (int p = 0; p < P; ++p) particles[p].rec[1] = HardRecord{root, -1, -1};
  ESS_hist.zeros(T); logZ_inc.zeros(T); logZ_run.zeros(T);
  logZ = 0.0;
  std::vector<double> lw(P);
  for (int t = 1; t <= T; ++t) {
    for (int p = 0; p < P; ++p) lw[p] = particles[p].logw;
    const double before = hard_lse(lw);
    sample_position(t);
    for (int p = 0; p < P; ++p) lw[p] = particles[p].logw;
    const double after = hard_lse(lw);
    logZ += after - before;
    logZ_inc[t - 1] = after - before;
    logZ_run[t - 1] = logZ;
    double sumsq = 0.0;
    for (int p = 0; p < P; ++p) sumsq += std::exp(2.0 * (lw[p] - after));
    ESS_hist[t - 1] = 1.0 / sumsq;
    resample();
    if ((t & 255) == 0) Rcpp::checkUserInterrupt();
  }
  for (int p = 0; p < P; ++p) lw[p] = particles[p].logw;
  const double norm = hard_lse(lw);
  weights.set_size(P);
  for (int p = 0; p < P; ++p) weights[p] = std::exp(lw[p] - norm);
}

// Draw leaf rates, evaluate the intensity at `grid`, and serialize the
// particles as lists of their nodes (region, depth, is_leaf, S, J, L, lambda, m).
Rcpp::List HardSMCtree::export_particles(const arma::mat& grid, arma::mat& lam_draws, arma::vec& loglik,
                                         arma::vec& lppd, arma::vec& integral) {
  const int np = grid.n_rows;
  lam_draws.set_size(np, P); loglik.set_size(P); lppd.set_size(P); integral.set_size(P);
  Rcpp::List out(P);
  for (int p = 0; p < P; ++p) {
    const HardParticle& par = particles[p];
    std::unordered_map<int, double> lambda;       // heap -> rate of a terminal node
    double ll = 0.0, tot = 0.0;
    Rcpp::List nodes(par.rec.size());
    int k = 0;
    for (const auto& kv : par.rec) {
      const HardNode& v = store[kv.second.node];
      const bool leaf = kv.second.S != 1;
      double lam = NA_REAL;
      int J = -1; double L = NA_REAL;
      if (leaf) {
        lam = R::rgamma(a + v.m, 1.0 / (b + v.area));
        lambda[kv.first] = lam;
        ll += logQ0(v.m, v.area);
        tot += lam * v.area;
      } else { J = v.cand_axis[kv.second.act]; L = v.cand_cut[kv.second.act]; }
      nodes[k++] = Rcpp::List::create(
        Rcpp::Named("region") = v.box, Rcpp::Named("depth") = v.depth,
        Rcpp::Named("is_leaf") = leaf, Rcpp::Named("S") = leaf ? 0 : 1,
        Rcpp::Named("J") = J, Rcpp::Named("L") = L,
        Rcpp::Named("lambda") = lam, Rcpp::Named("m") = v.m);
    }
    out[p] = nodes;
    loglik[p] = ll; integral[p] = tot;
    double lp = -tot;
    for (int g = 0; g < np; ++g) {
      int h = 1;
      for (;;) {
        const HardRecord& R = par.rec.at(h);
        if (R.S != 1) break;
        const HardNode& v = store[R.node];
        const int j = v.cand_axis[R.act];
        const bool go_left = j == d - 1 ? grid(g, j) <= v.cand_cut[R.act] : grid(g, j) < v.cand_cut[R.act];
        h = go_left ? 2 * h : 2 * h + 1;
      }
      lam_draws(g, p) = lambda[h];
      lp += std::log(lambda[h]);
    }
    lppd[p] = lp;
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List PPT_fit_SMC_shared(const arma::mat& pts, const arma::mat& grid, const arma::mat& region,
                              int max_depth, int P, int min_leaf_n, double resample_thresh,
                              double a, double b, double max_aspect_ratio, int cut_grid_n) {
  const int d = pts.n_cols;
  if (pts.n_rows == 0 || d == 0) Rcpp::stop("pts must be a non-empty matrix");
  if ((int)region.n_rows != d || region.n_cols != 2) Rcpp::stop("region must be a d by 2 matrix");
  if (max_depth < 1 || P < 1 || min_leaf_n < 1 || cut_grid_n < 1) Rcpp::stop("invalid SMC controls");
  HardSMCtree smc(pts, region, max_depth, min_leaf_n, cut_grid_n, P, a, b, 0.5, 2.0, max_aspect_ratio, resample_thresh);
  smc.sweep();

  arma::mat lam_draws; arma::vec loglik, lppd, integral;
  Rcpp::List particle = smc.export_particles(grid, lam_draws, loglik, lppd, integral);

  const int np = grid.n_rows;
  arma::vec lam_mean = lam_draws * smc.weights;
  arma::mat qmat(3, np);
  const arma::vec alphas = {0.025, 0.5, 0.975};
  std::vector<std::pair<double, double> > vw(P);
  for (int g = 0; g < np; ++g) {
    for (int p = 0; p < P; ++p) vw[p] = {lam_draws(g, p), smc.weights[p]};
    std::sort(vw.begin(), vw.end());
    double csum = 0.0; int k = 0;
    for (const auto& pr : vw) {
      csum += pr.second;
      while (k < 3 && csum >= alphas[k]) qmat(k, g) = pr.first, ++k;
      if (k == 3) break;
    }
    for (; k < 3; ++k) qmat(k, g) = vw.back().first;
  }
  return Rcpp::List::create(
    Rcpp::_["loglik"] = arma::dot(smc.weights, loglik), Rcpp::_["lppd"] = arma::dot(smc.weights, lppd),
    Rcpp::_["lambda"] = Rcpp::List::create(Rcpp::_["mean"] = lam_mean, Rcpp::_["median"] = qmat.row(1).t(),
                                           Rcpp::_["lower95"] = qmat.row(0).t(), Rcpp::_["upper95"] = qmat.row(2).t(),
                                           Rcpp::_["draws"] = lam_draws),
    Rcpp::_["particle"] = particle, Rcpp::_["weights"] = smc.weights,
    Rcpp::_["integrated_intensity"] = integral,
    Rcpp::_["resample_thresh"] = resample_thresh, Rcpp::_["ESS"] = smc.ESS_hist,
    Rcpp::_["logZ"] = smc.logZ, Rcpp::_["logZ_inc"] = smc.logZ_inc, Rcpp::_["logZ_run"] = smc.logZ_run,
    Rcpp::_["expanded"] = smc.n_expanded, Rcpp::_["resampled"] = smc.n_resampled);
}
