// ===========================================================================
// MCMCtree.cpp -- RcppArmadillo reversible-jump MCMC for the Poisson Process
// Tree model (comparison with the SMC). Leaf intensities lambda_l ~ Ga(a,b) are
// integrated out, so the chain moves over TREE space with target
//     pi(T|x) propto prior(T) * m(T),  m(T) = prod_leaf Q0(A_leaf),
//     Q0(A) = Gamma(n+a)/area^{n+a}    (Jeffreys a=0.5,b=0; general b>0 supported).
// Moves: GROW / PRUNE / CHANGE (Green 1995; Chipman-George-McCulloch 1998).
// Acceptance ratios validated against exact enumeration (see verification/).
//
//   Rcpp::sourceCpp("MCMCtree.cpp")        # standalone, or build as part of PPTree
// ===========================================================================
#ifndef _USE_Armadillo
#define _USE_Armadillo
#include <RcppArmadillo.h>
#include "tree_limits.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#endif
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <utility>
using namespace Rcpp;

// ---- leaf marginal Q0 (log) and depth-dependent split prior ----------------
static inline double mlogQ0(int n, double area, double a, double b) {
  if (b <= 0.0) return R::lgammafn(n + a) - (n + a) * std::log(area);      // Jeffreys/improper limit
  return R::lgammafn(a + n) - R::lgammafn(a) + a*std::log(b) - (a + n)*std::log(b + area);
}
static inline double mrho(int depth, double alpha, double eta) {
  double r = alpha * std::pow(1.0 + depth, -eta);
  if (r < 1e-12) r = 1e-12; if (r > 1.0 - 1e-12) r = 1.0 - 1e-12; return r;
}
static inline int runif_int(int k) { return (int)(R::unif_rand() * k); }     // 0..k-1

// Geometry and the hyperparameters are fixed for the lifetime of a fit. These
// caches are shared by copies of a node in the current and proposed trees only;
// rejected tree states are never retained in a global history cache.
struct MLocalSplit {
  int axis;
  double cut;
  double log_weight;                 // split prior times the two terminal leaves
  bool left_growable, right_growable;
};
struct MNodeCache {
  bool cuts_ready = false;
  std::vector< std::vector<double> > cuts;
  bool scores_ready = false;
  double leaf_log_weight = 0.0;
  std::vector<MLocalSplit> splits;
};

// ---- pointer-based tree node (so GROW/PRUNE can add/remove nodes) -----------
struct MNode {
  arma::mat box;       // d x 2 : [lower, upper] per axis
  arma::uvec idx;      // indices of points in this node
  int depth;
  MNode* L; MNode* R; MNode* parent; int axis; double cut;
  double lambda;
  std::shared_ptr<MNodeCache> cache;
  MNode(const arma::mat& b, const arma::uvec& i, int d)
    : box(b), idx(i), depth(d), L(nullptr), R(nullptr), parent(nullptr), axis(-1),
      cut(NA_REAL), lambda(NA_REAL), cache(std::make_shared<MNodeCache>()) {}
  ~MNode() { if (L) delete L; if (R) delete R; }
  bool   leaf() const { return L == nullptr; }
  double area() const { return arma::prod(box.col(1) - box.col(0)); }
  int    n()    const { return (int)idx.n_elem; }
};
static void collect_leaves(MNode* nd, std::vector<MNode*>& out) {
  if (nd->leaf()) out.push_back(nd); else { collect_leaves(nd->L, out); collect_leaves(nd->R, out); }
}
static void collect_prunable(MNode* nd, std::vector<MNode*>& out) {           // internal, both children leaves
  if (!nd->leaf()) {
    if (nd->L->leaf() && nd->R->leaf()) out.push_back(nd);
    collect_prunable(nd->L, out); collect_prunable(nd->R, out);
  }
}
static int maximum_leaf_depth(MNode* root) {
  std::vector<MNode*> leaves;
  collect_leaves(root, leaves);
  int depth = 0;
  for (size_t k = 0; k < leaves.size(); ++k)
    depth = std::max(depth, leaves[k]->depth);
  return depth;
}
static void collect_nodes(MNode* nd, std::vector<MNode*>& out) {
  out.push_back(nd);
  if (!nd->leaf()) {
    collect_nodes(nd->L, out);
    collect_nodes(nd->R, out);
  }
}
static Rcpp::List mtree_to_R(MNode* root) {
  std::vector<MNode*> nodes;
  collect_nodes(root, nodes);
  Rcpp::List out(nodes.size());
  for (size_t k = 0; k < nodes.size(); ++k) {
    MNode* node = nodes[k];
    out[k] = Rcpp::List::create(
      Rcpp::_["region"] = node->box,
      Rcpp::_["idx"] = node->idx + 1,
      Rcpp::_["depth"] = node->depth,
      Rcpp::_["is_empty"] = false,
      Rcpp::_["is_leaf"] = node->leaf(),
      Rcpp::_["S"] = node->leaf() ? 0 : 1,
      Rcpp::_["J"] = node->axis,
      Rcpp::_["L"] = node->cut,
      Rcpp::_["lambda"] = node->lambda
    );
  }
  return out;
}

// ---- valid quantile cuts at a node (filtered by min_leaf & a small buffer) --
static inline double qtype1(const arma::vec& s, double p) {                   // type-1 quantile (matches PPTree)
  int n = s.n_elem; double h = (n - 1) * p + 1.0; int i = std::max(1, (int)std::floor(h)); return s(i - 1);
}
static const std::vector< std::vector<double> >&
mvalid_cuts(MNode* nd, const arma::mat& pts, int min_leaf, int cut_grid_n) {
  if (nd->cache->cuts_ready) return nd->cache->cuts;
  int d = nd->box.n_rows, n = nd->n();
  std::vector< std::vector<double> > out(d);
  if (n < 2 * min_leaf) {
    nd->cache->cuts.swap(out);
    nd->cache->cuts_ready = true;
    return nd->cache->cuts;
  }
  for (int j = 0; j < d; ++j) {
    arma::vec colj = pts.col(j); arma::vec s = arma::sort(colj.elem(nd->idx));
    double lo = nd->box(j,0), hi = nd->box(j,1); double last = arma::datum::nan;
    for (int q = 0; q < cut_grid_n; ++q) {
      double p = cut_grid_n == 1 ? 0.05 : 0.05 + 0.90 * q / (double)(cut_grid_n - 1);
      double c = qtype1(s, p);
      if (c == last) continue; last = c;
      if (c > lo + 1e-3 && c < hi - 1e-3) {
        int nl = (int)(std::lower_bound(s.begin(), s.end(), c) - s.begin());  // # points < c
        if (nl >= min_leaf && n - nl >= min_leaf) out[j].push_back(c);
      }
    }
  }
  nd->cache->cuts.swap(out);
  nd->cache->cuts_ready = true;
  return nd->cache->cuts;
}
static inline bool has_valid(const std::vector< std::vector<double> >& vc) {
  for (size_t j = 0; j < vc.size(); ++j) if (!vc[j].empty()) return true; return false;
}
static inline bool msplittable(MNode* nd, const arma::mat& pts, int max_depth, int min_leaf, int cg) {
  return nd->leaf() && nd->depth < max_depth && has_valid(mvalid_cuts(nd, pts, min_leaf, cg));
}
static inline double mSfac(MNode* nd, const arma::mat& pts, int max_depth, int min_leaf, int cg,
                           double alpha, double eta) {
  if (nd->depth < max_depth && has_valid(mvalid_cuts(nd, pts, min_leaf, cg)))
    return 1.0 - mrho(nd->depth, alpha, eta);
  return 1.0;
}
static void make_children(MNode* nd, int j, double c, const arma::mat& pts, MNode*& Lo, MNode*& Ro) {
  arma::vec colj = pts.col(j); arma::vec xj = colj.elem(nd->idx);
  arma::uvec lm = arma::find(xj < c), rm = arma::find(xj >= c);
  arma::mat bl = nd->box; bl(j,1) = c; arma::mat br = nd->box; br(j,0) = c;
  std::unique_ptr<MNode> left(new MNode(bl, nd->idx.elem(lm), nd->depth + 1));
  std::unique_ptr<MNode> right(new MNode(br, nd->idx.elem(rm), nd->depth + 1));
  left->parent = nd;
  right->parent = nd;
  Lo = left.release();
  Ro = right.release();
}
// Draw leaf intensities lambda_l|x ~ Ga(a+n, b+area), evaluate them at the
// requested locations, and retain the exact point-process compensator and
// in-sample log likelihood for the same intensity draw.
static arma::vec predict_grid(MNode* root, const arma::mat& grid,
                              double a, double b, double& integral,
                              double& poisson_loglik) {
  int np = grid.n_rows, d = grid.n_cols; arma::vec out(np, arma::fill::zeros);
  std::vector<MNode*> lv; collect_leaves(root, lv);
  integral = 0.0;
  poisson_loglik = 0.0;
  for (size_t k = 0; k < lv.size(); ++k) {
    MNode* lf = lv[k]; double lam = R::rgamma(a + lf->n(), 1.0 / (b + lf->area()));
    lf->lambda = lam;
    integral += lam * lf->area();
    if (lf->n() > 0)
      poisson_loglik += lf->n() * std::log(std::max(lam, std::numeric_limits<double>::min()));
    arma::uvec inside(np, arma::fill::ones);                  // membership over all axes
    for (int j = 0; j < d; ++j) {
      arma::vec gj = grid.col(j);                             // materialise column (proven pattern)
      inside %= ((gj >= lf->box(j,0)) % (gj <= lf->box(j,1)));
    }
    arma::uvec in_idx = arma::find(inside);
    out.elem(in_idx).fill(lam);
  }
  poisson_loglik -= integral;
  return out;
}
static double tree_loglik(MNode* root, double a, double b) {
  std::vector<MNode*> lv; collect_leaves(root, lv); double s = 0;
  for (size_t k = 0; k < lv.size(); ++k) s += mlogQ0(lv[k]->n(), lv[k]->area(), a, b);
  return s;
}

// ---- Locally balanced informed MH -----------------------------------------
struct MInformedControl {
  const arma::mat& pts;
  int max_depth, min_leaf, cut_grid_n;
  double a, b, alpha, eta;
};

static const MNodeCache& mlocal_scores(MNode* nd, const MInformedControl& ctl) {
  MNodeCache& cache = *nd->cache;
  if (cache.scores_ready) return cache;
  const std::vector< std::vector<double> >& cuts =
    mvalid_cuts(nd, ctl.pts, ctl.min_leaf, ctl.cut_grid_n);
  const bool can_split = nd->depth < ctl.max_depth && has_valid(cuts);
  const double rho = mrho(nd->depth, ctl.alpha, ctl.eta);
  cache.leaf_log_weight = mlogQ0(nd->n(), nd->area(), ctl.a, ctl.b) +
    (can_split ? std::log1p(-rho) : 0.0);
  // Build into a temporary vector: an R interrupt cannot leave a partial cache.
  std::vector<MLocalSplit> scores;
  if (can_split) {
    for (size_t j = 0; j < cuts.size(); ++j) {
      for (size_t k = 0; k < cuts[j].size(); ++k) {
        if ((scores.size() & 63) == 0) Rcpp::checkUserInterrupt();
        MNode *left_raw, *right_raw;
        make_children(nd, j, cuts[j][k], ctl.pts, left_raw, right_raw);
        std::unique_ptr<MNode> left(left_raw), right(right_raw);
        bool gl = msplittable(left.get(), ctl.pts, ctl.max_depth, ctl.min_leaf, ctl.cut_grid_n);
        bool gr = msplittable(right.get(), ctl.pts, ctl.max_depth, ctl.min_leaf, ctl.cut_grid_n);
        double sl = gl ? std::log1p(-mrho(left->depth, ctl.alpha, ctl.eta)) : 0.0;
        double sr = gr ? std::log1p(-mrho(right->depth, ctl.alpha, ctl.eta)) : 0.0;
        MLocalSplit score;
        score.axis = j;
        score.cut = cuts[j][k];
        score.log_weight = std::log(rho) - std::log((double)cuts.size()) -
          std::log((double)cuts[j].size()) +
          mlogQ0(left->n(), left->area(), ctl.a, ctl.b) + sl +
          mlogQ0(right->n(), right->area(), ctl.a, ctl.b) + sr;
        score.left_growable = gl;
        score.right_growable = gr;
        scores.push_back(score);
      }
    }
  }
  cache.splits.swap(scores);
  cache.scores_ready = true;
  return cache;
}

static std::unique_ptr<MNode> mclone_tree(MNode* node, MNode* parent = nullptr) {
  std::unique_ptr<MNode> copy(new MNode(node->box, node->idx, node->depth));
  copy->axis = node->axis;
  copy->cut = node->cut;
  copy->lambda = node->lambda;
  copy->cache = node->cache;
  copy->parent = parent;
  if (!node->leaf()) {
    copy->L = mclone_tree(node->L, copy.get()).release();
    copy->R = mclone_tree(node->R, copy.get()).release();
  }
  return copy;
}

struct MNeighbor {
  int move;                           // 0 grow, 1 prune, 2 change
  int node_index;                     // preorder index in the source tree
  int axis;
  double cut;
  double log_target_ratio, log_q0, log_q0_reverse, log_ratio, log_weight;
  // A neighbor is deterministic given this current tree. Retain only its
  // scalar reverse normalizer, not its tree or candidate table. The vector
  // owning this scalar is discarded whenever the chain accepts a new state.
  double reverse_log_normalizer = std::numeric_limits<double>::quiet_NaN();
};
struct MNeighborhood {
  std::vector<MNeighbor> candidates;
  double log_normalizer = -std::numeric_limits<double>::infinity();
};

static void madd_neighbor(MNeighborhood& out, int move, int node_index,
                          int axis, double cut, double log_target_ratio,
                          double log_q0, double log_q0_reverse) {
  MNeighbor candidate;
  candidate.move = move;
  candidate.node_index = node_index;
  candidate.axis = axis;
  candidate.cut = cut;
  candidate.log_target_ratio = log_target_ratio;
  candidate.log_q0 = log_q0;
  candidate.log_q0_reverse = log_q0_reverse;
  candidate.log_ratio = log_target_ratio + log_q0_reverse - log_q0;
  candidate.log_weight = log_q0 + 0.5 * candidate.log_ratio;
  if (!std::isfinite(candidate.log_weight))
    Rcpp::stop("Non-finite informed proposal weight; check the data and hyperparameters.");
  out.candidates.push_back(candidate);
}

static MNeighborhood mneighborhood(MNode* root, const MInformedControl& ctl) {
  MNeighborhood out;
  std::vector<MNode*> nodes;
  collect_nodes(root, nodes);
  int growable = 0, prunable = 0;
  for (size_t i = 0; i < nodes.size(); ++i) {
    if (msplittable(nodes[i], ctl.pts, ctl.max_depth, ctl.min_leaf, ctl.cut_grid_n)) ++growable;
    if (!nodes[i]->leaf() && nodes[i]->L->leaf() && nodes[i]->R->leaf()) ++prunable;
  }
  const double log3 = std::log(3.0), logd = std::log((double)ctl.pts.n_cols);
  for (size_t i = 0; i < nodes.size(); ++i) {
    if ((i & 31) == 0) Rcpp::checkUserInterrupt();
    MNode* node = nodes[i];
    bool can_grow = msplittable(node, ctl.pts, ctl.max_depth, ctl.min_leaf, ctl.cut_grid_n);
    bool cherry = !node->leaf() && node->L->leaf() && node->R->leaf();
    if (!can_grow && !cherry) continue;
    const MNodeCache& local = mlocal_scores(node, ctl);
    const std::vector< std::vector<double> >& cuts = node->cache->cuts;
    if (can_grow) {
      bool parent_cherry = node->parent && node->parent->L->leaf() && node->parent->R->leaf();
      int reverse_prunable = prunable + 1 - (parent_cherry ? 1 : 0);
      for (size_t k = 0; k < local.splits.size(); ++k) {
        const MLocalSplit& split = local.splits[k];
        madd_neighbor(out, 0, i, split.axis, split.cut,
          split.log_weight - local.leaf_log_weight,
          -log3 - std::log((double)growable) - logd - std::log((double)cuts[split.axis].size()),
          -log3 - std::log((double)reverse_prunable));
      }
    } else {
      const MLocalSplit* previous = nullptr;
      for (size_t k = 0; k < local.splits.size(); ++k)
        if (local.splits[k].axis == node->axis && local.splits[k].cut == node->cut)
          previous = &local.splits[k];
      if (!previous) Rcpp::stop("Tree contains a split outside the RJ-MCMC candidate support.");
      int reverse_growable = growable - (previous->left_growable ? 1 : 0) -
        (previous->right_growable ? 1 : 0) + 1;
      madd_neighbor(out, 1, i, node->axis, node->cut,
        local.leaf_log_weight - previous->log_weight,
        -log3 - std::log((double)prunable),
        -log3 - std::log((double)reverse_growable) - logd - std::log((double)cuts[node->axis].size()));
      for (size_t k = 0; k < local.splits.size(); ++k) {
        const MLocalSplit& split = local.splits[k];
        // The original CHANGE move may redraw its current rule. It is a null
        // transition and is excluded, without renormalizing the remaining q0.
        if (&split == previous) continue;
        madd_neighbor(out, 2, i, split.axis, split.cut,
          split.log_weight - previous->log_weight,
          -log3 - std::log((double)prunable) - logd - std::log((double)cuts[split.axis].size()),
          -log3 - std::log((double)prunable) - logd - std::log((double)cuts[node->axis].size()));
      }
    }
  }
  if (!out.candidates.empty()) {
    double max_weight = out.candidates[0].log_weight, sum = 0.0;
    for (size_t i = 1; i < out.candidates.size(); ++i)
      max_weight = std::max(max_weight, out.candidates[i].log_weight);
    for (size_t i = 0; i < out.candidates.size(); ++i)
      sum += std::exp(out.candidates[i].log_weight - max_weight);
    out.log_normalizer = max_weight + std::log(sum);
  }
  return out;
}

static std::unique_ptr<MNode> mpropose_neighbor(MNode* root, const MNeighbor& candidate,
                                               const MInformedControl& ctl) {
  std::unique_ptr<MNode> proposed = mclone_tree(root);
  std::vector<MNode*> nodes;
  collect_nodes(proposed.get(), nodes);
  MNode* node = nodes[candidate.node_index];
  if (!node->leaf()) {
    delete node->L;
    delete node->R;
    node->L = nullptr;
    node->R = nullptr;
  }
  if (candidate.move == 1) {
    node->axis = -1;
    node->cut = NA_REAL;
  } else {
    make_children(node, candidate.axis, candidate.cut, ctl.pts, node->L, node->R);
    node->axis = candidate.axis;
    node->cut = candidate.cut;
  }
  return proposed;
}

static size_t mdraw_neighbor(const MNeighborhood& neighborhood) {
  const double u = R::unif_rand();
  double cdf = 0.0;
  for (size_t i = 0; i < neighborhood.candidates.size(); ++i) {
    cdf += std::exp(neighborhood.candidates[i].log_weight - neighborhood.log_normalizer);
    if (u < cdf) return i;
  }
  return neighborhood.candidates.size() - 1; // rounding at the last CDF endpoint
}

// ===========================================================================
static Rcpp::List mfit_tree(const arma::mat& pts, const arma::mat& grid, const arma::mat& region,
                        int niter = 4000, int burnin = 1000, int max_depth = 8, int min_leaf_n = 1,
                        int cut_grid_n = 50, double a = 0.5, double b = 0.0,
                        double alpha = 0.95, double eta = 2.0, int n_pred = 300,
                        bool informed = false) {
  ppt_checked_depth(max_depth);
  int N = pts.n_rows, d = pts.n_cols, np = grid.n_rows;
  std::unique_ptr<MNode> root_owner(new MNode(region, arma::regspace<arma::uvec>(0, N - 1), 0));
  MNode* root = root_owner.get();
  arma::ivec nl_trace(niter, arma::fill::zeros);
  arma::ivec md_trace(niter, arma::fill::zeros);
  arma::vec  ll_trace(niter, arma::fill::zeros);
  int pred_every = std::max(1, niter / std::max(1, n_pred));
  std::vector<arma::vec> preds;
  std::vector<double> integrals;
  std::vector<double> poisson_loglik;
  std::vector<Rcpp::List> tree_draws;
  long acc_g=0, acc_p=0, acc_c=0, tot_g=0, tot_p=0, tot_c=0;

  MInformedControl ctl{pts, max_depth, min_leaf_n, cut_grid_n, a, b, alpha, eta};
  MNeighborhood current;
  long reverse_cache_hits = 0, reverse_cache_misses = 0;
  arma::ivec neighborhood_size(niter, arma::fill::zeros);
  arma::vec log_normalizer(niter, arma::fill::zeros);
  if (informed) current = mneighborhood(root, ctl);

  for (int it = 0; it < niter + burnin; ++it) {
    if (informed) {
      if (!current.candidates.empty()) {
        MNeighbor& candidate = current.candidates[mdraw_neighbor(current)];
        if (candidate.move == 0) ++tot_g;
        else if (candidate.move == 1) ++tot_p;
        else ++tot_c;
        std::unique_ptr<MNode> proposed;
        MNeighborhood next;
        if (std::isnan(candidate.reverse_log_normalizer)) {
          ++reverse_cache_misses;
          proposed = mpropose_neighbor(root, candidate, ctl);
          next = mneighborhood(proposed.get(), ctl);
          candidate.reverse_log_normalizer = next.log_normalizer;
        } else ++reverse_cache_hits;
        if (std::log(R::unif_rand()) < current.log_normalizer - candidate.reverse_log_normalizer) {
          // Cached normalizers allow a repeated rejected proposal to avoid
          // constructing and scoring its tree. Materialize an accepted state
          // only after the test; these deterministic builders use no RNG.
          if (!proposed) {
            proposed = mpropose_neighbor(root, candidate, ctl);
            next = mneighborhood(proposed.get(), ctl);
          }
          if (candidate.move == 0) ++acc_g;
          else if (candidate.move == 1) ++acc_p;
          else ++acc_c;
          root_owner = std::move(proposed);
          root = root_owner.get();
          current = std::move(next);
        }
        // A rejection keeps the current tree, candidate table, and known
        // reverse normalizers. Memory is linear in the current neighborhood.
      }
    } else {
    double u = R::unif_rand();

    if (u < 1.0/3.0) {                                   // ---- GROW ----
      tot_g++;
      std::vector<MNode*> lv; collect_leaves(root, lv);
      std::vector<MNode*> B;
      for (size_t k=0;k<lv.size();++k) if (msplittable(lv[k],pts,max_depth,min_leaf_n,cut_grid_n)) B.push_back(lv[k]);
      if (!B.empty()) {
        MNode* ell = B[runif_int(B.size())]; int nB = B.size();
        const std::vector< std::vector<double> >& vc = mvalid_cuts(ell, pts, min_leaf_n, cut_grid_n);
        int j = runif_int(d);
        if (!vc[j].empty()) {
          double c = vc[j][runif_int(vc[j].size())];
          MNode *L, *R; make_children(ell, j, c, pts, L, R);
          double dm = mlogQ0(L->n(),L->area(),a,b) + mlogQ0(R->n(),R->area(),a,b) - mlogQ0(ell->n(),ell->area(),a,b);
          double SL = mSfac(L,pts,max_depth,min_leaf_n,cut_grid_n,alpha,eta);
          double SR = mSfac(R,pts,max_depth,min_leaf_n,cut_grid_n,alpha,eta);
          double rd = mrho(ell->depth, alpha, eta);
          double dp = std::log(rd) + std::log(SL) + std::log(SR) - std::log(1.0 - rd);
          ell->L=L; ell->R=R; ell->axis=j; ell->cut=c;                       // tentatively attach
          std::vector<MNode*> pr; collect_prunable(root, pr); int Wp = pr.size();
          double logA = dm + dp + std::log((double)nB) - std::log((double)Wp);
          if (std::log(R::unif_rand()) < logA) acc_g++;
          else { delete ell->L; delete ell->R; ell->L=nullptr; ell->R=nullptr; ell->axis=-1; }
        }
      }

    } else if (u < 2.0/3.0) {                            // ---- PRUNE ----
      tot_p++;
      std::vector<MNode*> pr; collect_prunable(root, pr);
      if (!pr.empty()) {
        MNode* v = pr[runif_int(pr.size())]; int nW = pr.size();
        double dm = mlogQ0(v->L->n(),v->L->area(),a,b) + mlogQ0(v->R->n(),v->R->area(),a,b) - mlogQ0(v->n(),v->area(),a,b);
        double SL = mSfac(v->L,pts,max_depth,min_leaf_n,cut_grid_n,alpha,eta);
        double SR = mSfac(v->R,pts,max_depth,min_leaf_n,cut_grid_n,alpha,eta);
        double rd = mrho(v->depth, alpha, eta);
        double dp = std::log(rd) + std::log(SL) + std::log(SR) - std::log(1.0 - rd);
        std::vector<MNode*> lv; collect_leaves(root, lv); int Bnow=0;         // splittable leaves now
        for (size_t k=0;k<lv.size();++k) if (msplittable(lv[k],pts,max_depth,min_leaf_n,cut_grid_n)) Bnow++;
        int sL = msplittable(v->L,pts,max_depth,min_leaf_n,cut_grid_n)?1:0;
        int sR = msplittable(v->R,pts,max_depth,min_leaf_n,cut_grid_n)?1:0;
        int sV = (v->depth<max_depth && has_valid(mvalid_cuts(v,pts,min_leaf_n,cut_grid_n)))?1:0;
        int Bp = Bnow - sL - sR + sV;                                         // splittable leaves after prune
        double logA = -dm - dp + std::log((double)nW) - std::log((double)Bp);
        if (std::log(R::unif_rand()) < logA) { delete v->L; delete v->R; v->L=nullptr; v->R=nullptr; v->axis=-1; acc_p++; }
      }

    } else {                                             // ---- CHANGE (terminal split) ----
      tot_c++;
      std::vector<MNode*> pr; collect_prunable(root, pr);
      if (!pr.empty()) {
        MNode* v = pr[runif_int(pr.size())];
        const std::vector< std::vector<double> >& vc = mvalid_cuts(v, pts, min_leaf_n, cut_grid_n);
        int j = runif_int(d);
        if (!vc[j].empty()) {
          double c = vc[j][runif_int(vc[j].size())];
          MNode *nL, *nR; make_children(v, j, c, pts, nL, nR);
          double dm = (mlogQ0(nL->n(),nL->area(),a,b)+mlogQ0(nR->n(),nR->area(),a,b))
                    - (mlogQ0(v->L->n(),v->L->area(),a,b)+mlogQ0(v->R->n(),v->R->area(),a,b));
          double dS = std::log(mSfac(nL,pts,max_depth,min_leaf_n,cut_grid_n,alpha,eta))
                    + std::log(mSfac(nR,pts,max_depth,min_leaf_n,cut_grid_n,alpha,eta))
                    - std::log(mSfac(v->L,pts,max_depth,min_leaf_n,cut_grid_n,alpha,eta))
                    - std::log(mSfac(v->R,pts,max_depth,min_leaf_n,cut_grid_n,alpha,eta));
          if (std::log(R::unif_rand()) < dm + dS) { delete v->L; delete v->R; v->L=nL; v->R=nR; v->axis=j; v->cut=c; acc_c++; }
          else { delete nL; delete nR; }
        }
      }
    }

    }
    if (it >= burnin) {
      int t = it - burnin;
      if (informed) {
        neighborhood_size(t) = current.candidates.size();
        log_normalizer(t) = current.log_normalizer;
      }
      std::vector<MNode*> lv; collect_leaves(root, lv);
      nl_trace(t) = (int)lv.size();
      md_trace(t) = maximum_leaf_depth(root);
      ll_trace(t) = tree_loglik(root, a, b);
      if (t % pred_every == 0) {
        double integral = 0.0, pp_loglik = 0.0;
        preds.push_back(predict_grid(root, grid, a, b, integral, pp_loglik));
        integrals.push_back(integral);
        poisson_loglik.push_back(pp_loglik);
        tree_draws.push_back(mtree_to_R(root));
      }
    }
    if ((it & 1023) == 0) Rcpp::checkUserInterrupt();
  }

  int K = preds.size(); arma::mat draws(np, K);
  for (int k = 0; k < K; ++k) draws.col(k) = preds[k];
  arma::vec integral_draws(K), poisson_loglik_draws(K);
  for (int k = 0; k < K; ++k) {
    integral_draws(k) = integrals[k];
    poisson_loglik_draws(k) = poisson_loglik[k];
  }
  arma::vec lam_mean = arma::mean(draws, 1), lam_med = arma::median(draws, 1);
  arma::vec lam_lo = arma::quantile(draws, arma::vec({0.025}), 1);
  arma::vec lam_hi = arma::quantile(draws, arma::vec({0.975}), 1);
  Rcpp::List trees(K);
  for (int k = 0; k < K; ++k) trees[k] = tree_draws[k];
  Rcpp::List result = Rcpp::List::create(
    Rcpp::_["lambda"] = Rcpp::List::create(Rcpp::_["mean"]=lam_mean, Rcpp::_["median"]=lam_med,
              Rcpp::_["lower95"]=lam_lo, Rcpp::_["upper95"]=lam_hi, Rcpp::_["draws"]=draws),
    Rcpp::_["nleaves"]=nl_trace, Rcpp::_["max_depth"]=md_trace,
    Rcpp::_["loglik"]=ll_trace,
    Rcpp::_["poisson_loglik"]=poisson_loglik_draws,
    Rcpp::_["integrated_intensity"]=integral_draws,
    Rcpp::_["tree_draws"]=trees,
    Rcpp::_["accept"]=Rcpp::NumericVector::create(Rcpp::_["grow"]=(double)acc_g/std::max(1L,tot_g),
              Rcpp::_["prune"]=(double)acc_p/std::max(1L,tot_p), Rcpp::_["change"]=(double)acc_c/std::max(1L,tot_c)),
    Rcpp::_["niter"]=niter, Rcpp::_["burnin"]=burnin);
  if (informed) result["informed"] = Rcpp::List::create(
    Rcpp::_["balancing"] = "sqrt",
    Rcpp::_["neighborhood_size"] = neighborhood_size,
    Rcpp::_["log_normalizer"] = log_normalizer,
    Rcpp::_["reverse_normalizer_cache_hits"] = (double)reverse_cache_hits,
    Rcpp::_["reverse_normalizer_cache_misses"] = (double)reverse_cache_misses);
  return result;
}

// [[Rcpp::export]]
Rcpp::List PPT_fit_MCMC(const arma::mat& pts, const arma::mat& grid, const arma::mat& region,
                        int niter = 4000, int burnin = 1000, double max_depth = 8, int min_leaf_n = 1,
                        int cut_grid_n = 50, double a = 0.5, double b = 0.0,
                        double alpha = 0.95, double eta = 2.0, int n_pred = 300) {
  const int depth = ppt_checked_depth(max_depth);
  return mfit_tree(pts, grid, region, niter, burnin, depth, min_leaf_n,
                   cut_grid_n, a, b, alpha, eta, n_pred, false);
}

// [[Rcpp::export]]
Rcpp::List PPT_fit_IMCMC(const arma::mat& pts, const arma::mat& grid, const arma::mat& region,
                         int niter = 4000, int burnin = 1000, double max_depth = 8, int min_leaf_n = 1,
                         int cut_grid_n = 50, double a = 0.5, double b = 0.0,
                         double alpha = 0.95, double eta = 2.0, int n_pred = 300) {
  const int depth = ppt_checked_depth(max_depth);
  return mfit_tree(pts, grid, region, niter, burnin, depth, min_leaf_n,
                   cut_grid_n, a, b, alpha, eta, n_pred, true);
}

// A deterministic diagnostic entry point for exact finite-state tests. Trees
// are encoded by rows (binary-heap node id, zero-based axis, cut); an empty
// 0-by-3 matrix denotes the root leaf. This is intentionally an internal API.
static void mdecode_splits(MNode* node, double id, const arma::mat& splits,
                           std::vector<bool>& used, const MInformedControl& ctl) {
  for (arma::uword k = 0; k < splits.n_rows; ++k) {
    if (splits(k, 0) != id) continue;
    if (node->depth >= ctl.max_depth) Rcpp::stop("Encoded tree exceeds max_depth.");
    int axis = (int)splits(k, 1);
    double cut = splits(k, 2);
    const std::vector< std::vector<double> >& cuts =
      mvalid_cuts(node, ctl.pts, ctl.min_leaf, ctl.cut_grid_n);
    if (axis < 0 || axis >= (int)cuts.size() ||
        std::find(cuts[axis].begin(), cuts[axis].end(), cut) == cuts[axis].end())
      Rcpp::stop("Encoded tree contains an invalid split.");
    make_children(node, axis, cut, ctl.pts, node->L, node->R);
    node->axis = axis;
    node->cut = cut;
    used[k] = true;
    mdecode_splits(node->L, 2.0 * id, splits, used, ctl);
    mdecode_splits(node->R, 2.0 * id + 1.0, splits, used, ctl);
    return;
  }
}
static void mencoded_rows(MNode* node, double id, std::vector<arma::rowvec>& rows) {
  if (node->leaf()) return;
  if (id > 4503599627370495.0) Rcpp::stop("Tree is too deep for numeric heap-id encoding.");
  rows.push_back(arma::rowvec({id, (double)node->axis, node->cut}));
  mencoded_rows(node->L, 2.0 * id, rows);
  mencoded_rows(node->R, 2.0 * id + 1.0, rows);
}
static arma::mat mencode_splits(MNode* root) {
  std::vector<arma::rowvec> rows;
  mencoded_rows(root, 1.0, rows);
  arma::mat out(rows.size(), 3);
  for (size_t k = 0; k < rows.size(); ++k) out.row(k) = rows[k];
  return out;
}
static double mlog_target(MNode* node, const MInformedControl& ctl) {
  const std::vector< std::vector<double> >& cuts =
    mvalid_cuts(node, ctl.pts, ctl.min_leaf, ctl.cut_grid_n);
  if (node->leaf()) return mlogQ0(node->n(), node->area(), ctl.a, ctl.b) +
    ((node->depth < ctl.max_depth && has_valid(cuts)) ?
       std::log1p(-mrho(node->depth, ctl.alpha, ctl.eta)) : 0.0);
  return std::log(mrho(node->depth, ctl.alpha, ctl.eta)) -
    std::log((double)cuts.size()) - std::log((double)cuts[node->axis].size()) +
    mlog_target(node->L, ctl) + mlog_target(node->R, ctl);
}
static double mheap_id(MNode* node) {
  std::vector<int> path;
  while (node->parent) {
    path.push_back(node == node->parent->R ? 1 : 0);
    node = node->parent;
  }
  double id = 1.0;
  for (std::vector<int>::reverse_iterator it = path.rbegin(); it != path.rend(); ++it)
    id = 2.0 * id + *it;
  return id;
}

// [[Rcpp::export]]
Rcpp::List PPT_IMCMC_transition(const arma::mat& pts, const arma::mat& region,
                                const arma::mat& splits, double max_depth = 8,
                                int min_leaf_n = 1, int cut_grid_n = 50,
                                double a = 0.5, double b = 0.0,
                                double alpha = 0.95, double eta = 2.0) {
  const int depth = ppt_checked_depth(max_depth);
  if (pts.n_rows == 0 || pts.n_cols == 0 || region.n_rows != pts.n_cols || region.n_cols != 2 ||
      !pts.is_finite() || !region.is_finite() || arma::any(region.col(1) <= region.col(0)) ||
      depth < 0 || min_leaf_n < 1 || cut_grid_n < 1 || !std::isfinite(a) || a <= 0 ||
      !std::isfinite(b) || b < 0 || !std::isfinite(alpha) || alpha <= 0 || alpha >= 1 ||
      !std::isfinite(eta) || eta < 0)
    Rcpp::stop("Invalid data, region, or informed-MH controls.");
  if (splits.n_cols != 3 || !splits.is_finite())
    Rcpp::stop("splits must be a finite numeric matrix with three columns.");
  for (arma::uword k = 0; k < splits.n_rows; ++k) {
    if (splits(k, 0) < 1 || splits(k, 0) > 4503599627370495.0 ||
        splits(k, 0) != std::floor(splits(k, 0)) ||
        splits(k, 1) < 0 || splits(k, 1) >= (double)pts.n_cols ||
        splits(k, 1) != std::floor(splits(k, 1)))
      Rcpp::stop("Invalid heap id or axis in splits.");
    for (arma::uword j = 0; j < k; ++j)
      if (splits(k, 0) == splits(j, 0)) Rcpp::stop("Duplicate heap id in splits.");
  }
  MInformedControl ctl{pts, depth, min_leaf_n, cut_grid_n, a, b, alpha, eta};
  std::unique_ptr<MNode> root(new MNode(region, arma::regspace<arma::uvec>(0, pts.n_rows - 1), 0));
  std::vector<bool> used(splits.n_rows, false);
  mdecode_splits(root.get(), 1.0, splits, used, ctl);
  for (size_t k = 0; k < used.size(); ++k)
    if (!used[k]) Rcpp::stop("Encoded split has a missing ancestor.");
  MNeighborhood current = mneighborhood(root.get(), ctl);
  std::vector<MNode*> nodes;
  collect_nodes(root.get(), nodes);
  int count = current.candidates.size();
  Rcpp::CharacterVector move(count);
  Rcpp::NumericVector node_id(count), cut(count), log_q0(count), log_q0_reverse(count),
    log_ratio(count), log_weight(count), proposal_probability(count), acceptance_probability(count),
    transition_probability(count), log_neighbor_normalizer(count), log_target_ratio(count);
  Rcpp::IntegerVector axis(count);
  Rcpp::List trees(count);
  double offdiagonal = 0.0, baseline_offdiagonal = 0.0;
  for (int k = 0; k < count; ++k) {
    const MNeighbor& candidate = current.candidates[k];
    std::unique_ptr<MNode> proposed = mpropose_neighbor(root.get(), candidate, ctl);
    MNeighborhood next = mneighborhood(proposed.get(), ctl);
    move[k] = candidate.move == 0 ? "grow" : (candidate.move == 1 ? "prune" : "change");
    node_id[k] = mheap_id(nodes[candidate.node_index]);
    axis[k] = candidate.axis;
    cut[k] = candidate.cut;
    log_q0[k] = candidate.log_q0;
    log_q0_reverse[k] = candidate.log_q0_reverse;
    log_ratio[k] = candidate.log_ratio;
    log_weight[k] = candidate.log_weight;
    proposal_probability[k] = std::exp(candidate.log_weight - current.log_normalizer);
    acceptance_probability[k] = std::exp(std::min(0.0, current.log_normalizer - next.log_normalizer));
    transition_probability[k] = proposal_probability[k] * acceptance_probability[k];
    log_neighbor_normalizer[k] = next.log_normalizer;
    log_target_ratio[k] = candidate.log_target_ratio;
    offdiagonal += transition_probability[k];
    baseline_offdiagonal += std::exp(candidate.log_q0);
    trees[k] = mencode_splits(proposed.get());
  }
  return Rcpp::List::create(
    Rcpp::_["log_target"] = mlog_target(root.get(), ctl),
    Rcpp::_["log_normalizer"] = current.log_normalizer,
    Rcpp::_["neighbors"] = Rcpp::DataFrame::create(
      Rcpp::_["move"] = move, Rcpp::_["node_id"] = node_id,
      Rcpp::_["axis"] = axis, Rcpp::_["cut"] = cut,
      Rcpp::_["log_q0"] = log_q0, Rcpp::_["log_q0_reverse"] = log_q0_reverse,
      Rcpp::_["log_ratio"] = log_ratio, Rcpp::_["log_weight"] = log_weight,
      Rcpp::_["proposal_probability"] = proposal_probability,
      Rcpp::_["acceptance_probability"] = acceptance_probability,
      Rcpp::_["transition_probability"] = transition_probability,
      Rcpp::_["log_neighbor_normalizer"] = log_neighbor_normalizer,
      Rcpp::_["log_target_ratio"] = log_target_ratio),
    Rcpp::_["trees"] = trees,
    Rcpp::_["self_probability"] = std::max(0.0, 1.0 - offdiagonal),
    Rcpp::_["baseline_self_probability"] = std::max(0.0, 1.0 - baseline_offdiagonal));
}
